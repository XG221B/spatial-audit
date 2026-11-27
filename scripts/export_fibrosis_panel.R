#!/usr/bin/env Rscript
# Export a bounded marker panel from the corrected archived Seurat object.
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% c(2L,3L)) stop("Expected the RDS, export directory and optional validated matrix projection")
dir.create(args[2L],recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratObject))
obj <- readRDS(args[1L])
meta <- obj@meta.data
message("Assays: ", paste(names(obj@assays),collapse=", "),
        ". Cell metadata: ", paste(colnames(meta),collapse=", "))
if (!"RNA" %in% names(obj@assays)) stop("RNA assay absent in corrected RDS")
projected <- length(args)==3L
if (projected) {
  if (names(obj@assays)[1L]!="RNA") stop("First serialized assay must be RNA")
  attrs <- obj@assays[["RNA"]]@counts[[1L]]
  dim <- attrs[["dim"]]
  names <- attrs[["dimnames"]]
  if (!identical(dim,c(343L,1630319L)) || !identical(lengths(names),dim))
    stop("Original serialized RNA dimnames disagree with checked projection")
  rows <- head(names[[1L]],32L)
  mat <- matrix(readBin(file.path(args[3L],"counts.i32"),"integer",n=32L*dim[2],size=4L,endian="little"),
                nrow=32L,dimnames=list(rows,names[[2L]]))
  library_size <- readBin(file.path(args[3L],"depth.f64"),"double",n=dim[2],size=8L,endian="little")
} else {
  mat <- SeuratObject::LayerData(obj, assay="RNA",layer="counts")
  library_size <- colSums(mat)
}
if (!identical(rownames(meta),colnames(mat)) || anyDuplicated(rownames(mat)) ||
    anyDuplicated(rownames(meta))) stop("Inconsistent or repeated feature or cell IDs")
fields <- c("sample","patient","disease_status","x_centroid","y_centroid","tma","run",
            "final_CT","final_lineage")
if (!all(fields %in% names(meta))) stop("Required sample, annotation or native centroid fields absent")
if (anyNA(meta[fields]) || any(!is.finite(as.matrix(meta[c("x_centroid","y_centroid")]))))
  stop("Missing specimen assignment, expression-derived cell class or native position")
if (!all(library_size==meta$nCount_RNA)) stop("Archived RNA count totals disagree with independently decoded all-target depths")
optional <- intersect(c("donor","donor_id","Donor_ID","status","sample_type","lineage"),names(meta))
meta <- meta[,c(fields,optional),drop=FALSE]
meta$library_size <- library_size
genes <- head(rownames(mat)[!grepl("^(NegControl|Unassigned|Blank|antisense|Deprecated)",
                                  rownames(mat),ignore.case=TRUE)],32L)
if (length(genes) != 32L) stop("Insufficient RNA targets")
writeLines(genes,file.path(args[2L],"features.txt"),useBytes=TRUE)
samples <- sort(unique(as.character(meta$sample)))
manifest <- data.frame(sample=samples, cells=integer(length(samples)), file=character(length(samples)))
for (i in seq_along(samples)) {
  ix <- which(as.character(meta$sample)==samples[i])
  data <- as.matrix(mat[genes,ix,drop=FALSE])
  if (any(!is.finite(data)) || any(data<0) || any(data!=floor(data))) stop("Invalid counts in registered panel")
  file <- sprintf("sample-%03d.csv.gz",i)
  con <- gzfile(file.path(args[2L],file),"wt")
  tryCatch(write.csv(cbind(cell_id=rownames(meta)[ix],meta[ix,,drop=FALSE],t(data)),con,row.names=FALSE),
           finally=close(con))
  manifest$cells[i] <- length(ix)
  manifest$file[i] <- file
}
write.csv(manifest,file.path(args[2L],"samples.csv"),row.names=FALSE)
