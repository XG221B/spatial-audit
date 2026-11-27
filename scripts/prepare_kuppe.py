"""Rebuild the pure fibrotic-zone versus donor-control Visium panel from Kuppe et al."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.schema import sha256, strict_read, write_json
from spatial_audit.visium_geometry import V1_LAYOUT as LAYOUT

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "10.5281/zenodo.6578047"
REPORT_ID = "zenodo-6578047-FZ"
GEOMETRY_SOURCE = "https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf"
PITCH = 100.0
GENES = [
    "ACTB",
    "GAPDH",
    "B2M",
    "RPLP0",
    "RPL13A",
    "RPS18",
    "RPS3",
    "RPL10",
    "EEF1A1",
    "PPIA",
    "HPRT1",
    "CALM1",
    "TTN",
    "MYH7",
    "MYH6",
    "TNNT2",
    "TNNI3",
    "MYL2",
    "MYL3",
    "ATP2A2",
    "PLN",
    "NPPA",
    "NPPB",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "PECAM1",
    "VWF",
    "ACTA2",
    "TAGLN",
]
COVARIATES = [
    "Adipocyte",
    "Cardiomyocyte",
    "Cycling.cells",
    "Endothelial",
    "Fibroblast",
    "Lymphoid",
    "Mast",
    "Myeloid",
    "Neuronal",
    "Pericyte",
    "vSMCs",
]
SAMPLE = re.compile(r"Visium_(FZ|control)_(P\d+)\.h5ad")
SELECT = "All single-zone FZ files and all donor controls in Zenodo record 6578047. Mixed-zone and other MI sections excluded. Fixed housekeeping, cardiomyocyte, matrix and vascular marker panel."


def physical(layout, barcodes, deposited):
    """Check each embedded coordinate against its vendor capture address."""
    if not layout.index.is_unique or not barcodes.is_unique:
        raise ValueError("Duplicate Visium barcode")
    ix = layout.index.get_indexer(barcodes)
    if (ix < 0).any():
        raise ValueError("Visium barcode not found in the manufacturer layout")
    rc = layout.iloc[ix].to_numpy(dtype=float)
    if not np.equal(rc, np.rint(rc)).all() or np.any((rc[:, 0] + rc[:, 1]) % 2):
        raise ValueError("Unexpected Visium grid parity")
    spread = np.ptp(rc, axis=0)
    if (spread == 0).any():
        raise ValueError("Cannot verify a one-dimensional spatial embedding")
    scaled = (rc - rc.min(axis=0)) / spread
    if deposited.shape != rc.shape:
        raise ValueError("Invalid archive display coordinate shape")
    transforms = []
    for axes in ([0, 1], [1, 0]):
        for reflected in ([False, False], [True, False], [False, True], [True, True]):
            display = np.where(reflected, 1 - scaled[:, axes], scaled[:, axes])
            if np.allclose(display, deposited, atol=1e-6, rtol=0):
                transforms.append({"array_axes": axes, "reflected_display_axes": reflected})
    if len(transforms) != 1:
        raise ValueError("Archive display coordinates disagree with vendor barcode map")
    # Doubled column addresses in a staggered hexagonal lattice.
    return np.column_stack([rc[:, 1] * PITCH / 2, rc[:, 0] * PITCH * np.sqrt(3) / 2]), transforms[0]


def receipt(path, a):
    if path.stat().st_size != a["bytes"]:
        raise ValueError(f"Incorrect archive size for {path.name}")
    with path.open("rb") as file:
        if hashlib.file_digest(file, "md5").hexdigest() != a["archive_checksum"].removeprefix(
            "md5:"
        ):
            raise ValueError(f"Incorrect official Zenodo MD5 for {path.name}")
    return {**a, "accession": ACCESSION, "sha256": sha256(path)}


def ordered_columns(matrix, indices):
    ordered = np.sort(indices)
    subset = matrix[:, ordered]
    dense = subset.toarray() if sparse.issparse(subset) else np.asarray(subset)
    return dense[:, np.searchsorted(ordered, indices)]


def shared_panel(assets, cache):
    common = set(GENES)
    missing = {}
    for asset in assets:
        a = ad.read_h5ad(cache / asset["name"], backed="r")
        try:
            if not a.var_names.is_unique:
                raise ValueError("Nonunique H5AD feature names")
            missing[asset["name"]] = sorted(set(GENES).difference(a.var_names))
            common.intersection_update(a.var_names)
        finally:
            a.file.close()
    if not common:
        raise ValueError("No candidate feature measured in all selected tissues")
    return [gene for gene in GENES if gene in common], missing


def project(path, layout, group_token, patient, genes):
    a = ad.read_h5ad(path, backed="r")
    try:
        title = (
            f"Visium-{group_token}_{patient}-Spatial multi-omic map of human myocardial infarction"
        )
        if a.uns["title"] != title or a.uns["X_normalization"] != "LogNormalization":
            raise ValueError(f"Archive specimen identity or normalization mismatch for {path.name}")
        expected = "MONDO:0005068" if group_token == "FZ" else "PATO:0000461"
        if set(a.obs["disease_ontology_term_id"]) != {expected}:
            raise ValueError(f"Inconsistent pathology ontology in {path.name}")
        coords, orientation = physical(layout, a.obs_names, a.obsm["X_spatial"])
        c = np.asarray(a.obs[COVARIATES], dtype=float)
        if (
            not np.isfinite(c).all()
            or (c < 0).any()
            or not np.allclose(c.sum(axis=1), 1, atol=1e-5)
        ):
            raise ValueError("Invalid deposited cell-type weights")
        ix = a.var_names.get_indexer(genes)
        if (ix < 0).any():
            raise ValueError("Missing common H5AD feature")
        # Sparse original matrix remains backed until this section is processed.
        x = ordered_columns(a.X, ix)
        if not np.isfinite(x).all() or (x < 0).any():
            raise ValueError("Invalid archived log-normalized expression")
        return coords, x, c, title, orientation
    finally:
        a.file.close()


def prepare(cache: Path, manifest: Path, evidence: Path, fetch_inputs: bool = False):
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    item = next(r for r in catalog["cohorts"] if r["id"] == ACCESSION)
    assets = [a for a in item["assets"] if SAMPLE.fullmatch(a["name"])]
    layout_file = cache / LAYOUT["name"]
    if fetch_inputs:
        from fetch_and_verify import fetch

        fetch(
            {
                "files": [{**a, "accession": ACCESSION} for a in assets]
                + [{**LAYOUT, "accession": "10x-Visium-v1-layout"}]
            },
            cache,
        )
    if sha256(layout_file) != LAYOUT["sha256"]:
        raise ValueError("Vendor barcode reference fingerprint mismatch")
    with tarfile.open(layout_file) as archive:
        handle = archive.extractfile("spatial/tissue_positions_list.csv")
        if handle is None:
            raise ValueError("Missing barcode grid in vendor archive")
        table = pd.read_csv(handle, header=None).set_index(0)[[2, 3]]
    receipts = [receipt(cache / asset["name"], asset) for asset in assets]
    genes, excluded = shared_panel(assets, cache)
    sections = []
    measurements = []
    for asset in assets:
        path = cache / asset["name"]
        match = SAMPLE.fullmatch(path.name)
        if match is None:
            raise ValueError("Unexpected section name")
        token, pid = match.groups()
        coords, matrix, c, title, orientation = project(path, table, token, pid, genes)
        output_dir = cache / "derived"
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / (path.stem + ".npz")
        with tempfile.NamedTemporaryFile(dir=output_dir, suffix=".npz", delete=False) as f:
            scratch = Path(f.name)
        try:
            np.savez_compressed(scratch, X=matrix, coords=coords, C=c, features=np.asarray(genes))
            os.replace(scratch, out)
        finally:
            scratch.unlink(missing_ok=True)
        group = "MI pure fibrotic zone" if token == "FZ" else "control donor"
        digest = sha256(out)
        common = {"patient_id": pid, "group": group, "section_id": path.stem}
        sections.append(
            {
                **common,
                "path": os.path.relpath(out, manifest.parent).replace(os.sep, "/"),
                "sha256": digest,
            }
        )
        measurements.append(
            {
                **common,
                "source_file": path.name,
                "archived_title": title,
                "archive_display_orientation": orientation,
                "validated_barcodes": len(coords),
                "section_bytes": out.stat().st_size,
                "section_sha256": digest,
            }
        )
    if len({row["patient_id"] for row in sections}) != len(sections):
        raise ValueError("Pure-region contrast requires unique heart donors")
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": REPORT_ID,
                    "accession": ACCESSION,
                    "kind": "public",
                    "platform": "Visium",
                    "source": item["source"],
                    "units": "um",
                    "preprocessing_evidence": evidence.name,
                    "covariates": COVARIATES,
                    "layer": "archived H5AD X; LogNormalization",
                    "measurement": "log_normalized",
                    "selection": SELECT,
                    "features": genes,
                    "files": sections,
                }
            ]
        },
    )
    write_json(
        evidence,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": item["source"],
            "paper": item["paper"],
            "files": receipts,
            "sections": measurements,
            "preprocessing": {
                "selection": SELECT,
                "candidate_genes": GENES,
                "target_genes": genes,
                "missing_candidates_by_asset": excluded,
                "covariates": COVARIATES,
                "composition": "The archive's spot-level cell-type weights sum to one. Estimated fractions share measurements with the outcomes.",
                "grid_layout_reference": LAYOUT,
                "grid_design_source": GEOMETRY_SOURCE,
                "spot_pitch_um": PITCH,
                "coordinates": "Match every barcode to the vendor capture-array layout; independently min-max rescale each array axis and compare to deposited X_spatial. Use only array-row/column addresses for nominal micrometer coordinates. Vendor image pixels and tissue identity are unused. No across-slide stage coordinates are available.",
            },
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve(), args.fetch)


if __name__ == "__main__":
    main()
