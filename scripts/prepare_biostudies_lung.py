"""Reconstruct the public human lung fibrosis contrast in BioStudies S-BSST1410."""

from __future__ import annotations

import argparse
import io
import json
import os
import tempfile
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.schema import sha256, strict_read, write_json
from spatial_audit.visium_geometry import PITCH, checked_visium_lattice

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "S-BSST1410"
DATA = "hs_visium_spaceranger_output.zip"
FRACTIONS = "cell2location_habermann2020.zip"
METADATA = "hs_visium_metadata.tsv"
GENES = (ROOT / "data/lung_marker_panel.txt").read_text(encoding="utf-8").splitlines()
SELECTION = "All human lung sections listed in the archive sample table. Each spot must have a positive filtered library, an in-tissue Space Ranger barcode and a deposited cell2location estimate. Fixed lung marker panel. The batch diagnostic represents procurement category from the primary article. Sequencing-run batch is unavailable."
PROCUREMENT = {
    **{f"HC_{i}": "procurement: postmortem" for i in range(1, 5)},
    **{f"IPF_{i}": "procurement: transplant surgery" for i in range(1, 5)},
}
COLUMNS = [
    "sample_id",
    "subject_alias",
    "condition",
    "slide_id",
    "slide_ca",
    "tissue_alias",
    "replicate",
    "species",
    "tissue",
]


def samples(path):
    tab = pd.read_csv(path, sep="\t", dtype=str)
    if tab[COLUMNS].isna().any().any() or tab["sample_id"].duplicated().any():
        raise ValueError("Missing or duplicated section metadata")
    if not tab["sample_id"].str.fullmatch(r"[A-Za-z0-9_-]+").all():
        raise ValueError("Nonportable section identifier")
    if set(tab["species"]) != {"human"} or set(tab["tissue"]) != {"lung"}:
        raise ValueError("Unexpected species or tissue")
    if (
        set(tab["condition"]) != {"IPF", "control"}
        or tab.groupby("subject_alias")["condition"].nunique().max() != 1
    ):
        raise ValueError("Inconsistent subject-to-group assignment")
    if not (tab["slide_id"] + "-" + tab["slide_ca"]).equals(tab["sample_id"]):
        raise ValueError(
            "Archive slide ID and capture-area fields do not reconstruct the section ID"
        )
    return tab


def abundance_table(archive, sid):
    src = f"cell2location_habermann2020/{sid}_spot_cell_abundances_5pc.csv"
    with archive.open(src) as f:
        frame = pd.read_csv(f).set_index("spot_id")
    if not frame.index.str.startswith(sid + "_").all() or not frame.index.is_unique:
        raise ValueError("Unexpected abundance table barcode")
    frame.index = frame.index.str.removeprefix(sid + "_")
    prefix = "q05cell_abundance_w_sf_"
    if not all(k.startswith(prefix) for k in frame.columns):
        raise ValueError("Unknown abundance table column")
    frame.columns = frame.columns.str.removeprefix(prefix)
    if not frame.columns.is_unique:
        raise ValueError("Ambiguous archived cell2location factors")
    return frame


def extract(archive, abundances, sid):
    matrix_member = f"{sid}/filtered_feature_bc_matrix.h5"
    with h5py.File(io.BytesIO(archive.read(matrix_member))) as file:
        h = file["matrix"]
        m = sparse.csc_matrix(
            (h["data"][:], h["indices"][:], h["indptr"][:]), shape=tuple(h["shape"][:])
        )
        if m.data.min(initial=0) < 0 or not np.issubdtype(m.dtype, np.integer):
            raise ValueError("Invalid stored count")
        symbols = pd.Index(h["features/name"].asstr()[:])
        barcodes = pd.Index(h["barcodes"].asstr()[:])
        if (
            not barcodes.is_unique
            or m.shape != (len(symbols), len(barcodes))
            or set(h["features/feature_type"].asstr()[:]) != {"Gene Expression"}
        ):
            raise ValueError("Invalid gene or observation axis")
        unique = ~symbols.duplicated(keep=False)
        lookup = {g: i for i, g in enumerate(symbols) if unique[i] and g in GENES}
        genes = [g for g in GENES if g in lookup]
        panel = m[[lookup[g] for g in genes], :].T.toarray()
        depth = np.asarray(m.sum(axis=0), dtype=np.int64).ravel()
        chemistry = str(file.attrs["chemistry_description"])
        software = str(file.attrs["software_version"])
    with archive.open(f"{sid}/spatial/tissue_positions_list.csv") as f:
        sites = pd.read_csv(
            f, header=None, names=["barcode", "in_tissue", "row", "col", "py", "px"]
        ).set_index("barcode")
    if not sites.index.is_unique or len(barcodes.difference(sites.index)):
        raise ValueError("Invalid position table or missing expression barcode")
    sites = sites.loc[barcodes]
    abundance = abundance_table(abundances, sid)
    has = barcodes.isin(abundance.index)
    positions = sites["in_tissue"].to_numpy() == 1
    valid = has & positions & (depth > 0)
    if not valid.any():
        raise ValueError("No aligned spots with counts and deconvolution")
    selected = barcodes[valid]
    a = abundance.loc[selected].to_numpy(dtype=float, copy=True)
    total = a.sum(axis=1)
    if not np.isfinite(a).all() or (a < 0).any() or (total <= 0).any():
        raise ValueError("Invalid or absent estimated cell abundance")
    a /= total[:, None]
    section = sites.loc[selected]
    scale = json.loads(archive.read(f"{sid}/spatial/scalefactors_json.json"))
    xy, proof = checked_visium_lattice(
        section[["row", "col"]], section[["px", "py"]], scale["spot_diameter_fullres"]
    )
    return (
        selected,
        xy,
        genes,
        panel[valid],
        depth[valid],
        abundance.columns.tolist(),
        a,
        {
            **proof,
            "matrix_member": matrix_member,
            "chemistry": chemistry,
            "software": software,
            "filtered_matrix_spots": len(barcodes),
            "missing_abundance_spots": int((~has).sum()),
            "abundance_only_spots": len(abundance.index.difference(barcodes)),
            "excluded_off_tissue_spots": int((~positions).sum()),
            "zero_depth_spots": int((depth == 0).sum()),
            "retained_spots": len(selected),
            "missing_or_ambiguous_targets": sorted(set(GENES) - set(genes)),
        },
    )


def prepare(
    cache: Path, manifest: Path, evidence: Path, count_manifest: Path | None, fetch_inputs=False
):
    item = next(
        c for c in strict_read(ROOT / "docs/data/catalog.json")["cohorts"] if c["id"] == ACCESSION
    )
    selected = [
        {**a, "accession": ACCESSION}
        for a in item["assets"]
        if a["name"] in {DATA, FRACTIONS, METADATA}
    ]
    if len(selected) != 3:
        raise ValueError("Expected all three archived products in the evidence catalog")
    if fetch_inputs:
        from fetch_and_verify import fetch

        fetch({"files": selected}, cache)
    receipts = []
    for a in selected:
        p = cache / a["name"]
        if p.stat().st_size != a["bytes"]:
            raise ValueError("BioStudies asset size mismatch")
        receipts.append({**a, "sha256": sha256(p)})
    tab = samples(cache / METADATA)
    if set(tab["subject_alias"]) != set(PROCUREMENT):
        raise ValueError("Archive subject aliases differ from the original procurement crosswalk")
    sections = []
    parts = []
    with zipfile.ZipFile(cache / DATA) as arc, zipfile.ZipFile(cache / FRACTIONS) as abund:
        expected = set(tab["sample_id"])
        members = [
            n.split("/")[0]
            for n in arc.namelist()
            if n.count("/") == 1 and n.endswith("/filtered_feature_bc_matrix.h5")
        ]
        if len(members) != len(set(members)) or set(members) != expected:
            raise ValueError(
                "Space Ranger archive and sample table do not describe the same sections"
            )
        for record in tab.to_dict(orient="records"):
            sid = record["sample_id"]
            barcodes, xy, genes, x, depth, names, c, proof = extract(arc, abund, sid)
            meta = {
                "patient_id": record["subject_alias"],
                "section_id": sid,
                "group": record["condition"],
                "slide_id": record["slide_id"],
                "batch": PROCUREMENT[record["subject_alias"]],
            }
            sections.append(
                {
                    **meta,
                    **proof,
                    "tissue_alias": record["tissue_alias"],
                    "replicate": record["replicate"],
                }
            )
            parts.append((meta, barcodes, xy, genes, x, depth, names, c))
    genes = [g for g in GENES if all(g in row[3] for row in parts)]
    names = parts[0][6]
    if not genes or any(set(z[6]) != set(names) for z in parts):
        raise ValueError("No shared marker panel or inconsistent cell abundance annotations")
    dest = cache / "derived"
    dest.mkdir(exist_ok=True)
    entries = []
    for (meta, barcodes, xy, targets, x, depth, labels, c), proof in zip(
        parts, sections, strict=True
    ):
        raw = x[:, [targets.index(g) for g in genes]]
        fractions = c[:, [labels.index(label) for label in names]]
        norm = np.log1p(10000 * raw.astype(float) / depth[:, None]).astype(np.float32)
        adata = ad.AnnData(
            X=norm,
            obs=pd.DataFrame({"library_size": depth}, index=barcodes),
            var=pd.DataFrame(index=genes),
            obsm={"spatial": xy, "composition": fractions},
            layers={"counts": raw},
        )
        out = dest / (meta["section_id"] + ".h5ad")
        with tempfile.NamedTemporaryFile(suffix=".h5ad", dir=dest, delete=False) as f:
            tmp = Path(f.name)
        try:
            adata.write_h5ad(tmp, compression="gzip")
            os.replace(tmp, out)
        finally:
            tmp.unlink(missing_ok=True)
        digest = sha256(out)
        proof.update(section_sha256=digest, section_bytes=out.stat().st_size)
        entries.append((out, {**meta, "sha256": digest}))
    base = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "kind": "public",
        "source": item["source"],
        "platform": "Visium",
        "units": "um",
        "preprocessing_evidence": evidence.name,
        "matrix_mode": "memory",
        "library_size_column": "library_size",
        "features": genes,
        "covariate_key": "composition",
        "selection": SELECTION,
        "batch_label": "procurement category (paper-derived)",
    }
    for p, extra in [
        (manifest, {"measurement": "log_normalized"}),
        (count_manifest, {"measurement": "counts", "layer": "counts"}),
    ]:
        if p is not None:
            rows = [
                {**row, "path": os.path.relpath(path, p.parent).replace(os.sep, "/")}
                for path, row in entries
            ]
            write_json(p, {"cohorts": [{**base, **extra, "files": rows}]})
    write_json(
        evidence,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": item["source"],
            "paper": item["paper"],
            "files": receipts,
            "sections": sections,
            "preprocessing": {
                "selection": SELECTION,
                "candidate_genes": GENES,
                "target_genes": genes,
                "cell2location_factors": names,
                "spot_pitch_um": PITCH,
                "expression": "Depth-normalized log1p using the complete filtered gene matrix per spot. Raw integer counts exported in a separate layer.",
                "composition": "Normalized marginal 5th-percentile abundances from the published Habermann et al. cell2location reference; the resulting relative weights are derived from the same spatial transcriptome.",
                "coordinate_validation": "Deposited Space Ranger spot array addresses cross-checked with adjacent pixels and the archived spot diameter. Coordinates remain local to each capture area.",
                "replication": "Patients identified by subject_alias; distinct slide/capture areas identified by sample_id. Tissue blocks and repeated technical sections are recorded under each donor. The engine weights sections equally within donors, so multiply sampled blocks may have greater within-donor influence.",
                "procurement_by_donor": PROCUREMENT,
                "acquisition": "slide_id from the archive. The batch field represents procurement category. The primary article's Results and Fig. 1 describe IPF 1–4 specimens collected during transplantation and HC 1–4 controls collected postmortem. Sequencing-run batches and cross-slide stage coordinates are unavailable.",
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
