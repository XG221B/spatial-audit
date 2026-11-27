#!/usr/bin/env Rscript
# One independent R process per cohort and method. All coordinates are in um.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_sota_benchmark.R METHOD INPUT_DIR OUTPUT_JSON")
method <- args[1L]
directory <- args[2L]
output <- args[3L]
suppressPackageStartupMessages(library(jsonlite))
packages <- c(edgeR="edgeR", DESeq2="DESeq2", DESpace2="DESpace",
              nnSVG="nnSVG", "SPARK-X"="SPARK", Cardinal="Cardinal")
pkg <- unname(packages[[method]])
version <- NULL
notes <- character(0)
set.seed(17)

# Keep the prespecified panel as the multiple-testing family even when a
# package filters a feature. Missing test statistics stay missing.
family_result <- function(feature, p, universe, scope) {
  if (length(feature) != length(p) || anyDuplicated(feature) || !all(feature %in% universe))
    stop("Unexpected feature mapping in package output")
  if (any(is.finite(p) & (p < 0 | p > 1))) stop("Package returned p outside [0,1]")
  raw <- rep(NA_real_, length(universe))
  raw[match(feature, universe)] <- p
  raw[!is.finite(raw)] <- NA_real_
  q <- p.adjust(raw, "BH", n=length(universe))
  list(scope=scope, n_family=length(universe), n_tests=sum(is.finite(q)),
       n_called=sum(q < 0.05, na.rm=TRUE),
       p=I(unname(raw)), q=I(unname(q)))
}

run_method <- function() {
  # loadNamespace errors retain their diagnostic (e.g. a missing DLL).
  if (!nzchar(system.file(package=pkg)))
    return(list(status="unavailable", reason=paste("R package unavailable:", pkg)))
  loadNamespace(pkg)
  version <<- as.character(packageVersion(pkg))
  if (method == "DESpace2" && packageVersion(pkg) < "2.0.0")
    return(list(status="unavailable", reason="Requires DESpace >= 2.0.0 with dsp_test"))
  manifest <- read_json(file.path(directory, "input.json"), simplifyVector=TRUE)
  features <- as.character(manifest$features)
  patients <- read.csv(file.path(directory, "patients.csv"), stringsAsFactors=FALSE,
                       check.names=FALSE, colClasses="character", na.strings=character(),
                       fileEncoding="UTF-8")
  patients$group <- factor(patients$group)
  rownames(patients) <- patients$patient_id
  sections <- lapply(seq_len(manifest$n_sections), function(i) {
    mat <- as(Matrix::readMM(file.path(directory, sprintf("s%05d.mtx", i))), "CsparseMatrix")
    obs <- read.csv(file.path(directory, sprintf("s%05d.csv", i)), check.names=FALSE,
                    colClasses=c("numeric","numeric","character","numeric"), na.strings=character(),
                    fileEncoding="UTF-8")
    stopifnot(nrow(mat) == length(features), ncol(mat) == nrow(obs))
    rownames(mat) <- features
    colnames(mat) <- paste0("s",i,"_",seq_len(ncol(mat)))
    list(X=mat, obs=obs, patient=manifest$patients[i], group=manifest$groups[i])
  })

  if (method %in% c("edgeR","DESeq2")) {
    pb <- matrix(0, length(features), nrow(patients),
                 dimnames=list(features,patients$patient_id))
    for (s in sections) pb[,s$patient] <- pb[,s$patient] + Matrix::rowSums(s$X)
    if (any(colSums(pb) == 0)) stop("Patient library is empty in the supplied feature panel")
    keep <- rowSums(pb) > 0
    counts <- pb[keep,,drop=FALSE]
    if (method == "edgeR") {
      y <- edgeR::DGEList(counts=counts, group=patients$group)
      y <- edgeR::calcNormFactors(y)
      design <- model.matrix(~group, patients)
      y <- edgeR::estimateDisp(y,design)
      fit <- edgeR::glmQLFit(y,design)
      tt <- edgeR::glmQLFTest(fit,coef=2L)$table
      p <- tt$PValue
      tested <- rownames(tt)
    } else {
      if (max(counts) > .Machine$integer.max) stop("Counts exceed DESeq2 integer range")
      storage.mode(counts) <- "integer"
      dds <- DESeq2::DESeqDataSetFromMatrix(counts,patients,~group)
      dds <- DESeq2::DESeq(dds,quiet=TRUE,sfType="poscounts")
      tt <- DESeq2::results(dds,independentFiltering=FALSE)
      tested <- rownames(tt)
      p <- tt$pvalue
    }
    panels <- list(family_result(tested,p,features,"cohort"))
  } else if (method == "DESpace2") {
    # DESpace 2.4 evaluates SCE helper names in the calling search path.
    suppressPackageStartupMessages(library(SingleCellExperiment))
    sections <- lapply(sections,function(s) {
      keep <- nzchar(s$obs$domain)
      s$X <- s$X[,keep,drop=FALSE]
      s$obs <- s$obs[keep,,drop=FALSE]
      s
    })
    mat <- do.call(cbind,lapply(sections, function(s) s$X))
    sample <- unlist(lapply(sections,function(s) rep(s$patient,ncol(s$X))))
    groups <- unlist(lapply(sections,function(s) rep(s$group,ncol(s$X))))
    obs <- do.call(rbind,lapply(sections,function(s) s$obs))
    spe <- SpatialExperiment::SpatialExperiment(
      assays=list(counts=mat), spatialCoords=as.matrix(obs[,c("x","y")]),
      colData=S4Vectors::DataFrame(patient=factor(sample), group=factor(groups),
        domain=factor(obs$domain), row.names=colnames(mat)))
    res <- DESpace::dsp_test(spe, cluster_col="domain", sample_col="patient",
      condition_col="group",test="QLF",min_counts=1,min_non_zero_spots=1,
      filter_gene=FALSE,filter_cluster=TRUE,verbose=FALSE)
    tab <- if (is.data.frame(res)) res else as.data.frame(res$gene_results)
    fc <- intersect(c("gene_id","gene"), names(tab))
    pc <- intersect(c("PValue","pvalue","P.Value","p_val"), names(tab))
    if (!length(pc)) stop(paste("Unsupported DESpace columns:",paste(names(tab), collapse=", ")))
    names_tested <- if (length(fc)) as.character(tab[[fc[1L]]]) else rownames(tab)
    panels <- list(family_result(names_tested,tab[[pc[1L]]],features,"cohort"))
  } else if (method == "Cardinal") {
    # meansTest summarizes spectra at the declared biological replicate.
    # Supply a feature-major matrix and request row-wise analyses explicitly.
    mat <- do.call(cbind,lapply(sections,function(s) as.matrix(s$X)))
    sample <- unlist(lapply(sections,function(s) rep(s$patient,ncol(s$X))))
    groups <- unlist(lapply(sections,function(s) rep(s$group,ncol(s$X))))
    obs <- data.frame(condition=factor(groups), patient=factor(sample))
    fit <- Cardinal::meansTest(mat,data=obs,fixed=~condition,samples=obs$patient,
             byrow=TRUE,verbose=FALSE,BPPARAM=BiocParallel::SerialParam())
    p <- S4Vectors::mcols(fit)$pvalue
    panels <- list(family_result(features,p,features,"cohort"))
  } else {
    panels <- lapply(seq_along(sections), function(i) {
      s <- sections[[i]]
      mat <- s$X
      normalized <- method=="nnSVG" && isTRUE(manifest$log_normalized)
      observed <- if (normalized) rep(TRUE,ncol(mat)) else if (method=="nnSVG") s$obs$depth>0 else Matrix::colSums(mat)>0
      mat <- mat[,observed,drop=FALSE]
      pos <- as.matrix(s$obs[observed,c("x","y"),drop=FALSE])
      keep <- Matrix::rowSums(mat)>0
      mat <- mat[keep,,drop=FALSE]
      scope <- paste0("section:",i)
      if (ncol(mat)<3L || !nrow(mat)) return(family_result(character(),numeric(),features,scope))
      if (method=="nnSVG") {
        if (normalized) {
          lognorm <- as.matrix(mat)
        } else {
          norm <- mat %*% Matrix::Diagonal(x=1e4 / s$obs$depth[observed])
          lognorm <- as.matrix(log1p(norm))
        }
        rownames(lognorm) <- rownames(mat)
        variable <- apply(lognorm, 1L, var) > 0
        if (!any(variable)) return(family_result(character(),numeric(),features,scope))
        # The matrix API expects genes in rows; return values include raw p-values.
        run_nn <- function(panel) nnSVG::nnSVG(panel,spatial_coords=pos,
                n_neighbors=min(10L,ncol(mat)-1L), n_threads=1L,verbose=FALSE,
                order=if (nrow(pos)<65L) "Sum_coords" else "AMMD",
                BPPARAM=BiocParallel::SerialParam())
        selected <- lognorm[variable,,drop=FALSE]
        tt <- tryCatch(as.data.frame(run_nn(selected)), error=function(e) NULL)
        if (is.null(tt)) {
          # Retain estimable genes when one covariance factorization fails.
          tested <- rownames(selected)
          p <- vapply(seq_len(nrow(selected)), function(j) {
            tryCatch({
              fit <- as.data.frame(run_nn(selected[j,,drop=FALSE]))
              if (!("pval" %in% names(fit)) || nrow(fit)!=1L) stop("unexpected nnSVG output")
              as.numeric(fit$pval)
            },error=function(e) {
              notes <<- unique(c(notes,paste(scope,tested[j],conditionMessage(e))))
              NA_real_
            })
          }, numeric(1L))
        } else {
          if (!("pval" %in% names(tt))) stop("nnSVG output lacks pval")
          tested <- rownames(tt)
          p <- tt$pval
        }
      } else {
        res <- SPARK::sparkx(mat,pos,numCores=1L,option="mixture",verbose=FALSE)
        tt <- as.data.frame(res$res_mtest)
        if (!("combinedPval" %in% names(tt))) stop("SPARK-X output lacks combinedPval")
        tested <- rownames(tt)
        p <- tt$combinedPval
      }
      family_result(tested,p,features,scope)
    })
  }
  list(status="ok", panels=panels)
}

result <- tryCatch(withCallingHandlers(run_method(),warning=function(w) {
  notes <<- unique(c(notes,conditionMessage(w)))
  invokeRestart("muffleWarning")
}),error=function(e) list(status="error",reason=conditionMessage(e)))
result$version <- version
result$r_version <- R.version.string
result$warnings <- I(head(notes,12L))
write_json(result,output,auto_unbox=TRUE,na="null",null="null",digits=NA)

