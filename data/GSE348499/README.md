# GSE348499 reconstruction

The source snapshot records the author-supplied patient title and treatment outcome for each GEO sample and the indexed URLs for its Xenium cell count matrix and cell centroid table. Both the original outcome string and the binary contrast label are retained. `feature_panel.txt` is the ordered gene family used by the shipped public benchmark. [GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348499) · [GEO sample-level SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348499&targ=gsm&form=text&view=full)

From the repository root:

```sh
python -m pip install -e ".[archives]"
python scripts/prepare_gse348499.py --raw-dir data-cache/xenium/raw --out data-cache/xenium/prepared --verify-geo --download
spatial-audit run --manifest data-cache/xenium/prepared/manifest-subset.json --out results
```

Use `manifest.json` for all donors, or add `--subset-only` at preparation time to download only the first three GSM accessions per outcome in lexical order. Missing files are stored in donor subdirectories. Both downloaded and retained files must match their pinned byte lengths and SHA-256 digests. `--verify-geo` independently re-reads current sample titles, clinical outcome characteristics and both URLs from GEO. Missing, duplicate or changed design information stops preparation.

The adapter filters for Gene Expression features, uses exact symbol lookup, and joins HDF5 count barcodes to the Parquet `cell_id` field. The output X is a raw integer-count matrix. The original x/y cell centroids are already in micrometres. Audit manifests include SHA-256 for each prepared H5AD. Add `--reference <directory>` to compare the selected counts and coordinates element by element with earlier per-GSM H5ADs and check donor labels and cell ordering. [Reconstruction evidence](../../results/public/GSE348499_reconstruction.json)

