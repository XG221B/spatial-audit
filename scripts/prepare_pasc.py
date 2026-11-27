"""Project the human GSE267226 count matrices onto the verified Visium V4 capture grid."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE267226"
ARCHIVE = ACCESSION + "_RAW.tar"
ARCHIVE_BYTES = 152811520
ARCHIVE_SHA256 = "777fdc8223da43898c4aeea9fe6b617b6bf2d73f2da3ce307e2469e89a93cf5a"
LAYOUT_REFERENCES = [
    (
        "2.0.0",
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",
        32680802,
        "cd636c52be5493031d3f58e56162906acdec16ab711f6d8cd05d18d315cc76e0",
    ),
    (
        "2.1.0",
        "CytAssist_FFPE_Protein_Expression_Human_Tonsil",
        34480354,
        "12027dcef2531604897ccfcaa7d64f2e8855ad9854d048a20bd6e0a77c933003",
    ),
]
GLOSSARY = "https://www.10xgenomics.com/support/software/space-ranger/latest/getting-started/space-ranger-glossary"
GEOMETRY_SOURCE = "https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf"
SPOT_PITCH = 100.0
GSM = ["GSM8262701", "GSM8262702", "GSM8262703", "GSM8262704", "GSM8262705"]
GENES = (ROOT / "data/lung_marker_panel.txt").read_text(encoding="utf-8").splitlines()
SELECTION = "All deposited human Visium samples in GSE267226. Predefined housekeeping, immune, epithelial and stromal panel restricted to unambiguous targeted gene symbols shared by every H5. Full-library counts per ten thousand followed by log1p are used for the main audit."


def layout_assets():
    refs = []
    for ver, name, size, digest in LAYOUT_REFERENCES:
        base = f"https://cf.10xgenomics.com/samples/spatial-exp/{ver}/{name}/{name}_"
        refs.append(
            {
                "name": name + "_spatial.tar.gz",
                "url": base + "spatial.tar.gz",
                "accession": "10x-Visium-V4-layout",
                "bytes": size,
                "sha256": digest,
                "manufacturer_assay_summary": base + "web_summary.html",
            }
        )
    return refs


def sample_assets():
    return [
        {
            "name": gsm + ".soft",
            "accession": ACCESSION,
            "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=brief",
        }
        for gsm in GSM
    ]


def vendor_layout(cache):
    tables = []
    for a in layout_assets():
        path = cache / a["name"]
        if sha256(path) != a["sha256"]:
            raise ValueError("Vendor Visium V4 layout fingerprint mismatch")
        with tarfile.open(path) as archive:
            member = archive.extractfile("spatial/tissue_positions.csv")
            if member is None:
                raise ValueError("Missing vendor Visium V4 tissue_positions.csv")
            frame = (
                pd.read_csv(member).set_index("barcode")[["array_row", "array_col"]].sort_index()
            )
        if not frame.index.is_unique or frame.duplicated().any() or frame.isna().any().any():
            raise ValueError("Ambiguous Visium V4 capture address")
        tables.append(frame)
    if not tables[0].equals(tables[1]) or len(tables[0]) != 4992:
        raise ValueError("Independent manufacturer Visium V4 designs differ")
    return tables[0]


def read_sample(path, gsm):
    soft = path.read_text(encoding="utf-8")

    def field(expr):
        matches = re.findall(expr, soft, re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f"Expected one matching GEO field in {gsm}")
        return matches[0].strip()

    title = field(r"^!Sample_title = (.+)$")
    disease = field(r"^!Sample_characteristics_ch1 = disease state: (.+)$")
    accession = field(r"^!Sample_geo_accession = (.+)$")
    if accession != gsm or field(r"^!Sample_organism_ch1 = (.+)$") != "Homo sapiens":
        raise ValueError("Unexpected human spatial sample identity")
    if disease == "control" and re.fullmatch(r"Control\d+", title):
        return title, "control donor"
    if disease == "COVID" and re.fullmatch(r"COVID\d+", title):
        return title, "PASC pulmonary fibrosis"
    raise ValueError(f"Unrecognized GEO lung group in {gsm}")


def assay(path, layout):
    with h5py.File(path) as handle:
        if handle.attrs["chemistry_description"] != "Visium V4 Slide":
            raise ValueError("A V4 barcode map requires the V4 chemistry declaration")
        m = handle["matrix"]
        obs = pd.Index(m["barcodes"].asstr()[:])
        idx = layout.index.get_indexer(obs)
        if not obs.is_unique or (idx < 0).any():
            raise ValueError("Duplicate or non-V4 assay barcode")
        rc = layout.iloc[idx].to_numpy(dtype=float)
        if np.any((rc[:, 0] + rc[:, 1]) % 2):
            raise ValueError("Invalid staggered V4 spot")
        xy = np.column_stack([rc[:, 1] * SPOT_PITCH / 2, rc[:, 0] * SPOT_PITCH * np.sqrt(3) / 2])
        fg = m["features"]
        if set(fg["feature_type"].asstr()[:]) != {"Gene Expression"} or set(
            fg["genome"].asstr()[:]
        ) != {"GRCh38"}:
            raise ValueError("Unexpected assay feature set")
        features = pd.Index(fg["name"].asstr()[:])
        unique = ~features.duplicated(keep=False)
        positions = {gene: i for i, gene in enumerate(features) if unique[i] and gene in GENES}
        raw = sparse.csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])
        )
        if (
            raw.shape != (len(features), len(obs))
            or not np.issubdtype(raw.data.dtype, np.integer)
            or (raw.data < 0).any()
        ):
            raise ValueError("Invalid counts or matrix axes")
        depth = np.asarray(raw.sum(axis=0), dtype=float).ravel()
        valid = depth > 0
        if not valid.any():
            raise ValueError("No positive-depth spot in filtered matrix")
        kept = [gene for gene in GENES if gene in positions]
        values = raw[[positions[g] for g in kept], :].T.toarray()
        info = {
            "chemistry": handle.attrs["chemistry_description"],
            "software": handle.attrs["software_version"],
            "library_ids": [x.decode() for x in handle.attrs["library_ids"]],
            "covered_array_barcodes": len(obs),
            "excluded_empty_libraries": int((~valid).sum()),
            "retained_spots": int(valid.sum()),
            "full_library_feature_count": len(features),
            "total_library_counts": int(depth.sum()),
            "missing_or_ambiguous_targets": sorted(set(GENES) - set(kept)),
        }
    return xy[valid], obs[valid], kept, values[valid], depth[valid], info


def prepare(
    cache: Path, manifest: Path, evidence: Path, count_manifest: Path | None, download=False
):
    study = next(
        c for c in strict_read(ROOT / "docs/data/catalog.json")["cohorts"] if c["id"] == ACCESSION
    )
    source_archive = next(a for a in study["assets"] if a["name"] == ARCHIVE)
    pin = {
        **source_archive,
        "accession": ACCESSION,
        "bytes": ARCHIVE_BYTES,
        "sha256": ARCHIVE_SHA256,
    }
    if download:
        from fetch_and_verify import fetch

        fetch({"files": [pin, *sample_assets(), *layout_assets()]}, cache)
    archive_path = cache / ARCHIVE
    if archive_path.stat().st_size != ARCHIVE_BYTES or sha256(archive_path) != ARCHIVE_SHA256:
        raise ValueError("GEO tar fingerprint mismatch")
    layout = vendor_layout(cache)
    samples = []
    measured = []
    sources = [pin]
    with tempfile.TemporaryDirectory(dir=cache) as tmp_name, tarfile.open(archive_path) as tar:
        members = [
            m
            for m in tar.getmembers()
            if m.isfile() and m.name.endswith("_filtered_feature_bc_matrix.h5")
        ]
        if {m.name.split("_", 1)[0] for m in members} != set(GSM) or len(members) != len(GSM):
            raise ValueError("Unexpected GEO set of human H5 files")
        for member in sorted(members, key=lambda m: m.name):
            gsm = member.name.split("_", 1)[0]
            sample_asset = next(a for a in sample_assets() if a["name"] == gsm + ".soft")
            soft = cache / sample_asset["name"]
            pid, group = read_sample(soft, gsm)
            sources.append({**sample_asset, "bytes": soft.stat().st_size, "sha256": sha256(soft)})
            extracted = Path(tmp_name) / (gsm + ".h5")
            reader = tar.extractfile(member)
            if reader is None:
                raise ValueError("Cannot read matrix inside tar")
            with reader, extracted.open("wb") as f:
                shutil.copyfileobj(reader, f)
            xy, obs, genes, x, depth, info = assay(extracted, layout)
            common = {"patient_id": pid, "section_id": gsm, "group": group}
            info.update(
                common,
                matrix_member=member.name,
                member_sha256=sha256(extracted),
                member_bytes=member.size,
            )
            measured.append(info)
            samples.append((common, xy, obs, genes, x, depth))
            extracted.unlink()
    shared = [g for g in GENES if all(g in sample[3] for sample in samples)]
    if not shared or len({s[0]["patient_id"] for s in samples}) != len(samples):
        raise ValueError("No common feature panel or repeated human lung identifiers")
    output = cache / "derived"
    output.mkdir(exist_ok=True)
    entries = []
    for (common, xy, obs, genes, x, depth), stats in zip(samples, measured, strict=True):
        counts = x[:, [genes.index(g) for g in shared]]
        normalized = np.log1p(10000 * counts.astype(float) / depth[:, None]).astype(np.float32)
        sec = ad.AnnData(
            X=normalized,
            obs=pd.DataFrame({"library_size": depth.astype(np.int64)}, index=obs),
            var=pd.DataFrame(index=shared),
            obsm={"spatial": xy},
            layers={"counts": counts},
        )
        dest = output / (common["section_id"] + ".h5ad")
        with tempfile.NamedTemporaryFile(dir=output, suffix=".h5ad", delete=False) as f:
            temp_path = Path(f.name)
        try:
            sec.write_h5ad(temp_path, compression="gzip")
            os.replace(temp_path, dest)
        finally:
            temp_path.unlink(missing_ok=True)
        digest = sha256(dest)
        stats.update(section_sha256=digest, section_bytes=dest.stat().st_size)
        entries.append((dest, {**common, "sha256": digest}))
    spec = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "kind": "public",
        "source": study["source"],
        "platform": "Visium",
        "units": "um",
        "preprocessing_evidence": evidence.name,
        "matrix_mode": "memory",
        "library_size_column": "library_size",
        "features": shared,
        "selection": SELECTION,
        "covariates": [],
    }
    for p, extra in [
        (manifest, {"measurement": "log_normalized"}),
        (count_manifest, {"layer": "counts", "measurement": "counts"}),
    ]:
        if p is not None:
            rows = [
                {**entry, "path": os.path.relpath(path, p.parent).replace(os.sep, "/")}
                for path, entry in entries
            ]
            write_json(p, {"cohorts": [{**spec, **extra, "files": rows}]})
    write_json(
        evidence,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": study["source"],
            "paper": study["paper"],
            "files": sources,
            "sections": measured,
            "preprocessing": {
                "selection": SELECTION,
                "candidate_genes": GENES,
                "genes": shared,
                "design_source": GLOSSARY,
                "geometry_source": GEOMETRY_SOURCE,
                "spot_pitch_um": SPOT_PITCH,
                "manufacturer_layouts": layout_assets(),
                "layout_size": len(layout),
                "layout_validation": "Every V4 assay barcode is unique and found on two independent, mutually identical official V4 reference layouts. Rows and columns of the manufactured hexagonal capture array supply relative positions, independent of which spots contain tissue. No tissue image or image-derived alignment was transferred across datasets.",
                "sample_mapping": "GEO sample titles and disease states are unique. Extended Data Fig. 1c–g explicitly gives the independent control and PASC-PF human lung counts for spatial profiling.",
                "coordinate_frame": "Manufacturer capture addresses are local to each array. Per-slide image orientation, histological registration and between-slide absolute acquisition coordinates are unavailable.",
                "composition": "No deposited cell-type fractions provided. The engine reports coverage of this dimension.",
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
