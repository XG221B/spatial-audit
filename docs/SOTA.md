# Reproducing the R comparison

The benchmark dispatches one fresh R process per method and cohort. Each process reads the same declared feature panel. Every output stores the patient assignments, section sizes, input file hashes, raw p-values, harmonized BH values, package version, and scope of each test. Input manifests are the same ones used by the Python audit.

The released comparison includes a table of measured calls by cohort. Each cell gives the registered family and the number of estimable hypotheses. The companion JSON includes the feature-level probabilities and sampling records. [Measured package results](../results/public/sota.md) · [Machine-readable results](../results/public/sota.json)

## Environment

Use a supported Bioconductor release and its matching R version. Install the following in an R session. SPARK-X is distributed in the authors' **SPARK** package. Building it from source requires the compiler suite for your R installation. DESpace2 is provided by **DESpace** from version 2.0. An installed library directory can be selected using the standard R_LIBS environment variable.

```r
install.packages(c("BiocManager", "remotes", "jsonlite"))
BiocManager::install(c("edgeR", "DESeq2", "DESpace", "nnSVG", "Cardinal",
                       "SpatialExperiment", "SingleCellExperiment", "BiocParallel"))
remotes::install_github("xzhoulab/SPARK@a8b4bf27b804604dfda53da42992f100b8e4e727",
                        upgrade = "never")
```

Pass multiple manifests to run multiple cohorts. Each output contains a single provenance kind. Declare the RNA count layer in the manifest, or set `"measurement":"counts"` for an already-standardized NPZ count matrix. DESpace2 also needs domain annotations supplied independently of the group contrast. Set `--rscript` to the installed executable when it is outside PATH.

## Reproduce the cross-platform suite

The public registry gives the accession, adapter, analysis ID, expected manifest path relative to a staging directory, and within-section SVG sampling cap. Prepare each manifest under the registered path using the referenced adapter. GSE220442 and GSE348499 write manifests inside their output directories. For Stereo-seq, post-acute pulmonary fibrosis, Mayr and the BioStudies lung cohort, select the adapter's count manifest. The Kuppe comparison uses the deposited log-normalized layer. The dashboard separately includes full-panel and full-cohort analyses of studies with a registered smaller cross-method comparison. [Machine-readable suite registry](data/sota_registry.json) · [Archive preparation and evidence](COHORTS.md)

```sh
python scripts/prepare_gse220442.py --download --raw-dir cache/ad --out manifests/GSE220442
python scripts/prepare_gse348499.py --download --verify-geo --raw-dir cache/xenium --out manifests/GSE348499
python scripts/run_sota_benchmark.py --release-inputs manifests --out results/public/sota.json --markdown results/public/sota.md
```

Run the last command after staging the other registered manifests. For a single prepared cohort, pass `--manifest manifests/cosmx.json --cap-for GSE327879=100` in place of `--release-inputs manifests`. Spot limits sample nuclei, spots or bins without replacement using deterministic independent random streams for each section. Both within-section methods receive the same sampled positions. Patient and domain contrasts retain the full registered observation set. A fresh R process isolates each method; the Python driver opens one input manifest at a time.

The R implementation lives in `scripts/run_sota_benchmark.R`. The Python driver prepares a private Matrix Market directory and passes its path to a snapshot of the R script. This protects running workers from edits to the checkout. Its DESpace adapter attaches the SingleCellExperiment namespace for the package's unqualified helper. The Windows driver maps an inherited POSIX C UTF-8 locale to the system UTF-8 locale for Unicode library paths. If a package is missing, inspect `.libPaths()` in the same R installation as `Rscript`, then install the listed package there or set `R_LIBS`. For SPARK installation on Windows use the Rtools compiler distribution matching R. A successful package check is `stopifnot("sparkx" %in% getNamespaceExports("SPARK"), "dsp_test" %in% getNamespaceExports("DESpace"))`.

## Interpreting each row

| Method | Hypothesis and sampling unit | Output family |
|---|---|---|
| edgeR | TMM-normalized patient pseudobulk counts differ between groups, quasi-likelihood F test | Features in cohort |
| DESeq2 | Patient pseudobulk counts differ between groups, Wald test with positive-count size factors | Features in cohort |
| DESpace2 | Group differences across annotated spatial domains, domain–patient pseudobulk, QLF | Features in cohort |
| nnSVG | A spatial Gaussian-process component explains variability within one section of log-normalized data | Features within each section |
| SPARK-X | Kernel-based association of transcript counts with spatial positions within one section | Features within each section |
| Cardinal | Patient means of the registered ion intensity differ between clinical groups, likelihood-ratio test | Ions in cohort |

Within-section SVG rows do not test a group label. The SVG summary counts feature–section opportunities, applies BH separately in each section and aggregates the number of calls. No significance test pools physical coordinates across slides. SPARK-X drops positions with zero registered panel counts. nnSVG retains an observation whenever its normalization library has positive depth. Features with no counts in a section cannot contribute a spatial test. The registered panel remains the denominator. Null entries in the per-feature vectors and the tested fraction show how many hypotheses the model could evaluate.

For nnSVG on count assays, the driver uses `log1p(10000 * counts / depth)`. If the standardized input supplies archived depth in an NPZ array, a hashed sidecar or a mapped H5AD observation field, the driver uses that full library. A backed H5AD supplies the sum over all columns of the declared count layer when no depth field was registered. For a projected input with a truncated feature panel and no depth source, it records panel-closed normalization. A log-normalized measurement can be passed directly to the nnSVG likelihood; count models require integer counts. The exact depth source and transformation appear in each method row. DESpace excludes unlabeled or unassigned positions (including the GSE220442 noise category), records their count, and requires at least two shared domains. Numerical factorization failures within an nnSVG batch are retried one gene at a time; remaining non-estimable feature-section pairs retain null p-values in the original family.

The spot gap is defined for marginal contrasts only. The RNA reference ranks supplied spot counts, while edgeR and DESeq2 normalize patient libraries using the registered panel. The resulting call-rate difference reflects changes in both sampling unit and normalization. It does not estimate a false discovery rate from unknown ground truth. When the panel is small, library normalization itself can be sensitive to panel composition. DESpace2, nnSVG and SPARK-X have different null hypotheses; the corresponding spot-gap field is unavailable.

Resource measurements include R process startup, package imports, reads of exported matrices and model fitting. Python loading, Wilcoxon and reference patient-level tests, and Matrix Market export take place before timing. RSS is the maximum observed sum of resident sizes across the process tree at 30 ms sampling intervals. Short peaks can be missed and shared pages can be counted in multiple processes. Timings support comparisons under the recorded workload, host load and environment. Missing packages, incompatible assay types, model failures and timeouts have distinct statuses. The CLI returns a nonzero exit for errors and timeouts after writing the diagnostic report.

## Primary implementations

- Patient count models. [edgeR](https://bioconductor.org/packages/edgeR/) · [DESeq2](https://bioconductor.org/packages/DESeq2/)
- Differential spatial patterns. [DESpace](https://bioconductor.org/packages/DESpace/)
- Within-section spatial variability. [nnSVG](https://bioconductor.org/packages/nnSVG/) · [SPARK-X source](https://github.com/xzhoulab/SPARK/tree/a8b4bf27b804604dfda53da42992f100b8e4e727)
- Imaging-mass-spectrometry inference. [Cardinal](https://bioconductor.org/packages/Cardinal/)

