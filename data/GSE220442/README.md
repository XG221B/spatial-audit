# GSE220442 reconstruction

The adapter joins the two GEO supplementary files at donor and spot level. The archived `patientID` and `category` columns supply donor identities and group labels. The barcode suffix is checked against the donor's archive directory. Cortical layer and white matter indicators represent histological compartments. [GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442) · [Original study](https://doi.org/10.1186/s40478-022-01494-6)

From the repository root:

```sh
python scripts/prepare_gse220442.py --raw-dir data-cache/raw --out data-cache/visium --download
spatial-audit run --manifest data-cache/visium/manifest.json --out results
spatial-audit run --manifest data-cache/visium/manifest-full-panel.json --out results
```

The first command downloads missing assets at URLs recorded in the catalog and checks both files against the retained archive fingerprints. To reuse a local copy, give its directory through `--raw-dir`. To check a reconstruction against existing `GSE220442_ready` files, add `--reference <directory>`. The resulting `dataset_manifest.json` records observed group assignments, calibration and file fingerprints. Archive fingerprints were computed independently from the retained bytes.

`feature_panel.tsv` freezes Ensembl IDs, labels and order from the prepared benchmark. The upstream selection used average donor detection of at least 5%, the top 2,500 genes among eligible genes ranked by mean CP10K, and curated AD pathway genes. The panel file fully specifies both analysis families. The comparator manifest takes its first 32 genes. The full-panel manifest takes all frozen genes and requests an in-memory matrix for efficient repeated column access on this cohort. Both analyses use raw counts. Library normalization for the reconstructed H5AD X uses the complete original library before selecting the panel. The audit maps the archived staggered array row and column to a hexagonal lattice at the nominal Visium 100 µm center spacing. The reconstruction receipt compares these distances with the original image-based transformation for each section.

The portable reconstruction evidence records exact, elementwise comparison with the local reference across all selected features and all sections. It covers donor and group labels, spot IDs, raw counts, normalized X, the legacy image-coordinate transform and audit covariates. The output physical coordinates use the independently calibrated array addresses. [Reconstruction evidence](../../results/public/GSE220442_reconstruction.json)


Space Ranger documents the image-diameter convention separately from the capture oligo diameter. The physical coordinate field in this release is calibrated from the lattice. [10x Genomics spatial outputs](https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs)
