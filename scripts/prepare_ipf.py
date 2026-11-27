"""Audit all Visium sections in the Mayr et al. human IPF H5AD with donor-level replication."""

from __future__ import annotations

import argparse
import hashlib
import os
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse

from spatial_audit.schema import sha256, strict_read, write_json
from spatial_audit.visium_geometry import PITCH, PITCH_TOLERANCE, V1_LAYOUT, checked_visium_lattice

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "10.5281/zenodo.10012934"
V5_LAYOUT: dict = {
    "name": "CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz",
    "url": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz",
    "bytes": 34400239,
    "sha256": "1125e5cb0372c9e14d0d3023d50518df87201f3a2a32d987f4c8a157534f3421",
}
REPORT_ID = "zenodo-10012934"
GENES = (ROOT / "data/lung_marker_panel.txt").read_text(encoding="utf-8").splitlines()
GEOMETRY_SOURCE = "https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf"
SUPPLEMENT = "https://pmc-oa-opendata.s3.amazonaws.com/PMC11313858.1/sciadv.adl5473_sm.pdf"
SELECTION = "All Visium and CytAssist sections and every QC-retained spot in the deposited H5AD. Each biological donor has unit weight. Predefined lung marker panel shared with the post-acute lung fibrosis adapter; ambiguous or absent symbols excluded."
META = [
    "sampleID",
    "patient",
    "treatment",
    "VisiumSlides",
    "technology",
    "array_row",
    "array_col",
    "Niche_NMF",
]


def check_design(obs):
    if obs[META].isna().any().any():
        raise ValueError("Incomplete deposited sample, location or domain annotation")
    sample = obs[
        ["sampleID", "patient", "treatment", "VisiumSlides", "technology"]
    ].drop_duplicates()
    if sample["sampleID"].duplicated().any():
        raise ValueError("Inconsistent donor or acquisition metadata within a section")
    # Supplementary Figure S1A reports which donors contribute adjacent sections.
    dual = {"Donor_Control_1", "Donor_Control_2", "Donor_IPF_1", "Donor_IPF_2"}
    standard = {"Donor_Control_3", "Donor_Control_4", "Donor_IPF_3"}
    if set(sample["patient"]) != dual | standard:
        raise ValueError("The deposited donor set differs from Figure S1")
    for patient, group in sample.groupby("patient", observed=True):
        tech = ["CytAssist", "Visium"] if patient in dual else ["Visium"]
        treatment = "IPF" if str(patient).startswith("Donor_IPF_") else "Healthy"
        if sorted(group["technology"]) != tech or set(group["treatment"]) != {treatment}:
            raise ValueError("Deposited case or duplicate-section mapping conflicts with Figure S1")
    return sample


def counts_and_depth(node, columns, block=512):
    if node.attrs["encoding-type"] != "csr_matrix":
        raise ValueError("Expected the deposited CSR count layer")
    n, p = map(int, node.attrs["shape"])
    ptr = node["indptr"][:]
    if len(ptr) != n + 1 or ptr[0] != 0 or np.any(np.diff(ptr) < 0) or ptr[-1] != len(node["data"]):
        raise ValueError("Invalid CSR count-layer pointers")
    counts = np.empty((n, len(columns)), dtype=np.int32)
    depth = np.empty(n, dtype=np.int64)
    for start in range(0, n, block):
        stop = min(start + block, n)
        left, right = int(ptr[start]), int(ptr[stop])
        vals = node["data"][left:right]
        idx = node["indices"][left:right]
        if len(idx) != len(vals) or (idx < 0).any() or (idx >= p).any():
            raise ValueError("Invalid stored feature index")
        if (
            not np.isfinite(vals).all()
            or (vals < 0).any()
            or not np.equal(vals, np.rint(vals)).all()
            or (vals > np.iinfo(np.int32).max).any()
        ):
            raise ValueError("Noninteger or unsupported deposited count entry")
        m = sparse.csr_matrix(
            (vals.astype(np.int64), idx, ptr[start : stop + 1] - left), shape=(stop - start, p)
        )
        counts[start:stop] = m[:, columns].toarray()
        depth[start:stop] = np.asarray(m.sum(axis=1)).ravel()
    if np.any(depth <= 0):
        raise ValueError("Zero-depth spot after archive QC")
    return counts, depth


def prepare(
    cache: Path, manifest: Path, evidence: Path, count_manifest: Path | None, fetch_inputs=False
):
    item = next(
        c for c in strict_read(ROOT / "docs/data/catalog.json")["cohorts"] if c["id"] == ACCESSION
    )
    asset = next(a for a in item["assets"] if a["name"].endswith(".h5ad"))
    source = cache / asset["name"]
    if fetch_inputs:
        from fetch_and_verify import fetch

        fetch(
            {
                "files": [
                    {**asset, "accession": ACCESSION},
                    {**V1_LAYOUT, "accession": "10x-Visium-v1-layout"},
                    {**V5_LAYOUT, "accession": "10x-Visium-V5-layout"},
                ]
            },
            cache,
        )
    vendor_layouts = {}
    for tech, ref in (("Visium", V1_LAYOUT), ("CytAssist", V5_LAYOUT)):
        layout_file = cache / ref["name"]
        if sha256(layout_file) != ref["sha256"]:
            raise ValueError("Vendor array layout fingerprint mismatch")
        with tarfile.open(layout_file) as archive:
            member = (
                "spatial/tissue_positions_list.csv"
                if tech == "Visium"
                else "spatial/tissue_positions.csv"
            )
            handle = archive.extractfile(member)
            if handle is None:
                raise ValueError("Missing tissue positions in vendor reference")
            if tech == "Visium":
                table = pd.read_csv(handle, header=None).set_index(0)[[2, 3]]
            else:
                table = pd.read_csv(handle).set_index("barcode")[["array_row", "array_col"]]
        if not table.index.is_unique or table.duplicated().any():
            raise ValueError("Ambiguous manufacturer capture addresses")
        vendor_layouts[tech] = table
    if source.stat().st_size != asset["bytes"]:
        raise ValueError("Zenodo H5AD length mismatch")
    with source.open("rb") as f:
        if hashlib.file_digest(f, "md5").hexdigest() != asset["archive_checksum"].removeprefix(
            "md5:"
        ):
            raise ValueError("Zenodo H5AD MD5 mismatch")
    receipt = {**asset, "accession": ACCESSION, "sha256": sha256(source)}
    with h5py.File(source) as h:
        obs_node = h["obs"]
        obs_names = pd.Index(read_elem(obs_node[obs_node.attrs["_index"]])).astype(str)
        obs = pd.DataFrame({k: np.asarray(read_elem(obs_node[k])) for k in META}, index=obs_names)
        design = check_design(obs)
        symbols = pd.Index(read_elem(h["var/SYMBOL"])).astype(str)
        unique = ~symbols.duplicated(keep=False)
        positions = {sym: j for j, sym in enumerate(symbols) if unique[j] and sym in GENES}
        genes = [g for g in GENES if g in positions]
        if not genes:
            raise ValueError("No assay features match the prespecified lung panel")
        columns = [positions[g] for g in genes]
        counts, filtered_depth = counts_and_depth(h["layers/counts"], columns)
        ensembl = pd.Index(read_elem(h["var/ENSEMBL"])).astype(str)
        raw_ids = pd.Index(read_elem(h["raw/var/ENSEMBL"])).astype(str)
        if not raw_ids.is_unique:
            raise ValueError("Ambiguous raw gene IDs")
        raw_columns = raw_ids.get_indexer(ensembl[columns])
        if (raw_columns < 0).any():
            raise ValueError("One or more panel genes missing from the raw assay")
        raw_counts, depth = counts_and_depth(h["raw/X"], raw_columns)
        if not np.array_equal(raw_counts, counts) or (filtered_depth > depth).any():
            raise ValueError("Filtered counts disagree with the complete raw assay")
        prior_depth = np.asarray(read_elem(obs_node["total_counts"]), dtype=float)
        if (
            not np.isfinite(prior_depth).all()
            or not np.equal(prior_depth, np.rint(prior_depth)).all()
            or np.any(prior_depth < depth)
        ):
            raise ValueError("Invalid archive prefilter QC library size")
        stored = sparse_dataset(h["X"])
        archived = np.empty(counts.shape, dtype=np.float32)
        for start in range(0, len(obs), 512):
            stop = min(start + 512, len(obs))
            block = stored[start:stop, :]
            if not sparse.issparse(block):
                raise TypeError("Expected a sparse block from archive X")
            archived[start:stop] = sparse.csr_matrix(block)[:, columns].toarray()
        reconstructed = np.log1p(10000 * counts.astype(float) / depth[:, None])
        if not np.isfinite(archived).all() or (archived < 0).any():
            raise ValueError("Invalid archived log-normalized expression")
        close = np.isclose(archived, reconstructed, rtol=2e-6, atol=1e-6)
        abundance = read_elem(h["obsm/q05_cell_abundance_w_sf"])
        if not isinstance(abundance, pd.DataFrame):
            raise TypeError("Expected a labeled cell2location factor table")
        if not abundance.index.equals(obs.index) or not abundance.columns.is_unique:
            raise ValueError(
                "cell2location abundance rows or factors disagree with archive observation axis"
            )
        prefix = "q05cell_abundance_w_sf_"
        if not all(k.startswith(prefix) for k in abundance.columns):
            raise ValueError("Unknown cell abundance factors")
        c = abundance.to_numpy(dtype=float, copy=True)
        total = c.sum(axis=1)
        if not np.isfinite(c).all() or (c < 0).any() or (total <= 0).any():
            raise ValueError("Invalid estimated abundance for one or more spots")
        c /= total[:, None]
        fractions = [k.removeprefix(prefix) for k in abundance.columns]
        pix = h["obsm/spatial"][:]
        if pix.shape != (len(obs), 2) or len(counts) != len(obs):
            raise ValueError("Coordinates or count matrix differ from the observation axis")
        entries = []
        sections = []
        out = cache / "derived"
        out.mkdir(exist_ok=True)
        for record in design.to_dict(orient="records"):
            name = str(record["sampleID"])
            sel = obs["sampleID"].to_numpy() == name
            if not np.any(sel) or not obs.index[sel].is_unique:
                raise ValueError("Duplicate or missing spots within section")
            diam = h[f"uns/spatial/{name}/scalefactors/spot_diameter_fullres"][()]
            reference = vendor_layouts[record["technology"]]
            barcodes = obs.index[sel].str.removeprefix(name + "_")
            if not barcodes.is_unique or len(barcodes.difference(reference.index)):
                raise ValueError("Section has absent or ambiguous vendor barcodes")
            xy, calibration = checked_visium_lattice(
                obs.loc[sel, ["array_row", "array_col"]],
                pix[sel],
                diam,
                reference.loc[barcodes].to_numpy(),
            )
            local = counts[sel]
            norm = archived[sel]
            metadata = pd.DataFrame(
                {
                    "domain": obs.loc[sel, "Niche_NMF"].astype(str),
                    "library_size": depth[sel],
                    "archive_prefilter_total_counts": prior_depth[sel].astype(np.int64),
                },
                index=obs.index[sel],
            )
            sec = ad.AnnData(
                X=norm,
                obs=metadata,
                var=pd.DataFrame(index=genes),
                obsm={"spatial": xy, "composition": c[sel]},
                layers={"counts": local},
            )
            dest = out / (name + ".h5ad")
            if dest.parent != out:
                raise ValueError("Nonportable section name")
            with tempfile.NamedTemporaryFile(dir=out, suffix=".h5ad", delete=False) as tmp:
                scratch = Path(tmp.name)
            try:
                sec.write_h5ad(scratch, compression="gzip")
                os.replace(scratch, dest)
            finally:
                scratch.unlink(missing_ok=True)
            ident = {
                "patient_id": str(record["patient"]),
                "group": str(record["treatment"]),
                "section_id": name,
                "slide_id": str(record["VisiumSlides"]),
                "batch": str(record["technology"]),
            }
            digest = sha256(dest)
            entries.append((dest, {**ident, "sha256": digest}))
            sections.append(
                {
                    **ident,
                    **calibration,
                    "spots": int(sel.sum()),
                    "section_bytes": dest.stat().st_size,
                    "section_sha256": digest,
                    "retained_assay_count_sum": int(depth[sel].sum()),
                    "raw_feature_axis": len(raw_ids),
                    "raw_normalization_checked_entries": int(close[sel].size),
                    "raw_normalization_mismatched_entries": int((~close[sel]).sum()),
                    "raw_normalization_max_absolute_difference": float(
                        np.max(np.abs(norm - reconstructed[sel]))
                    ),
                }
            )
    spec = {
        "id": REPORT_ID,
        "accession": ACCESSION,
        "kind": "public",
        "source": item["source"],
        "platform": "Visium",
        "units": "um",
        "preprocessing_evidence": evidence.name,
        "features": genes,
        "covariate_key": "composition",
        "domain": "domain",
        "matrix_mode": "memory",
        "library_size_column": "library_size",
        "selection": SELECTION,
        "batch_label": "assay technology",
    }
    for path, extra in [
        (manifest, {"measurement": "log_normalized"}),
        (count_manifest, {"measurement": "counts", "layer": "counts"}),
    ]:
        if path is not None:
            files = [
                {**entry, "path": os.path.relpath(f, path.parent).replace(os.sep, "/")}
                for f, entry in entries
            ]
            write_json(path, {"cohorts": [{**spec, **extra, "files": files}]})
    write_json(
        evidence,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": item["source"],
            "paper": item["paper"],
            "supplement": SUPPLEMENT,
            "files": [receipt],
            "sections": sections,
            "preprocessing": {
                "selection": SELECTION,
                "candidate_genes": GENES,
                "target_genes": genes,
                "missing_or_ambiguous_targets": sorted(set(GENES) - set(genes)),
                "cell2location_factors": fractions,
                "composition": "Per-spot relative weights computed by dividing every posterior 5th-percentile abundance by their sum. These are normalized marginal quantiles, and are estimated using the same spatial transcriptome.",
                "expression": "Original H5AD X retained for every panel entry. Counts from layers/counts match the same Ensembl IDs in raw/X. Library size is summed across all deposited raw/X genes before feature selection. Every archived X entry is compared against log1p(counts times 10000 divided by the raw/X library sum). Sections report discrepancies with absolute tolerance 1e-6 and relative tolerance 2e-6. The earlier deposited obs.total_counts is preserved as separate QC metadata.",
                "coordinates": "Space Ranger array_row and array_col mapped to the nominal hexagonal pitch. Pixel distances are independently checked with archived spot_diameter_fullres and Space Ranger's theoretical display diameter.",
                "geometry_source": GEOMETRY_SOURCE,
                "vendor_layouts": {"Visium": V1_LAYOUT, "CytAssist": V5_LAYOUT},
                "registered_images": "Image-scale disagreements are retained as sensitivity diagnostics. Nominal physical centers are determined independently by the validated manufacturer barcode-address map.",
                "spot_pitch_um": PITCH,
                "max_relative_pitch_discrepancy": PITCH_TOLERANCE,
                "replication": "Original patient, sampleID and treatment mapping agrees with Supplementary Fig. S1A and Table S2. Technical sections are nested within patient.",
                "acquisition": "VisiumSlides is the deposited slide label. The batch diagnostic here uses deposited technology (direct-placement Visium or CytAssist). Absolute instrument stage coordinates unavailable.",
            },
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--counts-manifest", type=Path)
    p.add_argument("--provenance", type=Path, required=True)
    p.add_argument("--fetch", action="store_true")
    a = p.parse_args()
    prepare(
        a.cache.resolve(),
        a.manifest.resolve(),
        a.provenance.resolve(),
        a.counts_manifest.resolve() if a.counts_manifest else None,
        a.fetch,
    )


if __name__ == "__main__":
    main()
