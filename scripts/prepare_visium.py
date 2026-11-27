"""Convert Visium v1/v2 6.5-mm filtered counts and archived spot positions to the audit format.

Use --cache for the registered GSE202011 lesion vs healthy contrast. For additional
studies supply --map (TSV with section_id, patient_id, group, matrix, positions and
optional batch/slide_id), --id and --source. Paths in the TSV are relative to the TSV.
The matrix must be a 10x filtered feature-barcode HDF5 output. Every barcode must match.
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.schema import sha256, write_json

ROOT = Path(__file__).resolve().parents[1]
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
GEOMETRY = "https://www.10xgenomics.com/support/software/space-ranger/latest/getting-started/space-ranger-glossary"
PANEL = [
    "ACTB",
    "GAPDH",
    "RPLP0",
    "PPIA",
    "B2M",
    "EEF1A1",
    "RPS18",
    "RPL13A",
    "EPCAM",
    "KRT5",
    "KRT14",
    "KRT1",
    "KRT10",
    "KRT16",
    "KRT17",
    "S100A7",
    "S100A8",
    "S100A9",
    "LOR",
    "FLG",
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM",
    "VIM",
    "PECAM1",
    "VWF",
    "PTPRC",
    "CD3D",
    "CD68",
    "LYZ",
    "HLA-DRA",
]


def parse_spot_map(raw: bytes) -> pd.DataFrame:
    """Validate the complete barcode-to-address map and tissue mask."""
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    frame = pd.read_csv(io.BytesIO(raw), header=None)
    if str(frame.iloc[0, 0]).lower() == "barcode":
        frame = frame.iloc[1:]
    if frame.shape[1] != 6 or len(frame) != 4992:
        raise ValueError("Expected the full 6.5-mm Visium spot map with six columns")
    frame.columns = ["barcode", "in_tissue", "row", "col", "image_y", "image_x"]
    frame = frame.set_index("barcode")
    if not frame.index.is_unique:
        raise ValueError("Visium barcodes must be unique per capture area")
    values = frame[["row", "col", "in_tissue"]].to_numpy(dtype=float)
    row, col, tissue = values.T
    if (
        not np.isfinite(values).all()
        or np.any(values != np.floor(values))
        or np.any((row < 0) | (row >= 78) | (col < 0) | (col >= 128) | ((row + col) % 2 != 0))
        or len(np.unique(values[:, :2], axis=0)) != len(frame)
        or not np.isin(tissue, [0, 1]).all()
    ):
        raise ValueError("Invalid Visium 6.5-mm hexagonal array geometry")
    return frame[["row", "col", "in_tissue"]].astype(int)


def positions_um(raw: bytes | pd.DataFrame, barcodes: np.ndarray) -> np.ndarray:
    """Use a section map or a verified chemistry-specific barcode address table."""
    frame = parse_spot_map(raw) if isinstance(raw, bytes) else raw
    if not pd.Index(barcodes).is_unique:
        raise ValueError("Duplicate barcodes in matrix")
    indices = frame.index.get_indexer(barcodes)
    if (indices < 0).any() or (
        "in_tissue" in frame and (frame.in_tissue.to_numpy()[indices] != 1).any()
    ):
        raise ValueError("Missing or off-tissue barcodes in filtered matrix")
    rc = frame.iloc[indices]
    return np.column_stack((50 * rc.col, 50 * np.sqrt(3) * rc.row))


def read_panel(
    matrix, positions: bytes | pd.DataFrame, panel: list[str], chemistry: str | None = None,
    *, with_depth=False
):
    with h5py.File(matrix, "r") as file:
        if chemistry is not None and file.attrs.get("chemistry_description") != chemistry:
            raise ValueError("Matrix chemistry does not match the vendor address reference")
        h = file["matrix"]
        feature = pd.Index(h["features/name"].asstr()[:])
        types = h["features/feature_type"].asstr()[:]
        selected = []
        for gene in panel:
            loc = np.flatnonzero((feature == gene) & (types == "Gene Expression"))
            if len(loc) != 1:
                raise ValueError(f"Expected exactly one RNA feature named {gene}")
            selected.append(int(loc[0]))
        cols = h["barcodes"].asstr()[:]
        coords = positions_um(positions, cols)
        counts = sparse.csc_matrix(
            (h["data"][:], h["indices"][:], h["indptr"][:]), shape=tuple(h["shape"][:])
        )
        x = counts[selected, :].T.toarray()
        if not np.isfinite(x).all() or (x < 0).any() or np.any(x != np.floor(x)):
            raise ValueError("Non-count measurement in HDF5 RNA panel")
        if with_depth:
            depth = np.asarray(counts[types == "Gene Expression"].sum(axis=0)).ravel()
            return x, coords, depth
        return x, coords


def geo_records(raw: str) -> list[dict]:
    rows = []
    for part in re.split(r"(?m)^\^SAMPLE = ", raw)[1:]:
        gsm = part.splitlines()[0].strip()
        characteristics = dict(
            re.findall(r"(?m)^!Sample_characteristics_ch1 = ([^:\n]+): (.+)$", part)
        )
        # A single filtered matrix and one donor are required per GEO record.
        files = re.findall(r"(?m)^!Sample_supplementary_file(?:_\d+)? = (.+\.h5)$", part)
        if len(files) != 1 or not re.fullmatch(r"GSM\d+", gsm):
            raise ValueError(f"Invalid filtered-matrix mapping for {gsm}")
        disease = characteristics["disease state"]
        status = characteristics["disease status"]
        donor = characteristics["individual"]
        selected = (disease == "Psoriasis (PsO)" and status == "Lesional") or (
            disease == "control" and status == "Healthy"
        )
        rows.append(
            {
                "section_id": gsm,
                "patient_id": donor,
                "group": disease,
                "batch": characteristics["sequencing run / batch"],
                "matrix_name": PurePosixPath(files[0]).name,
                "selected": selected,
                "status": status,
            }
        )
    # PsA donor numbers overlap with PsO in GEO. Exclude PsA before interpreting IDs.
    seen = set()
    for row in sorted(rows, key=lambda v: v["section_id"]):
        if row["selected"]:
            key = (row["patient_id"], row["group"])
            row["selected"] = key not in seen
            seen.add(key)
    return rows


def package_section(out, manifest, metadata, x, coords, panel, library_size=None):
    name = metadata["section_id"]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
        raise ValueError("Use a portable section ID")
    dest = out / (name + ".npz")
    with tempfile.NamedTemporaryFile(dir=out, suffix=".npz", delete=False) as handle:
        tmp = Path(handle.name)
    try:
        depth = {} if library_size is None else {"library_size": library_size}
        np.savez_compressed(tmp, coords=coords, X=x, features=np.asarray(panel), **depth)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    return {
        k: metadata[k]
        for k in ["section_id", "patient_id", "group", "batch", "slide_id"]
        if metadata.get(k)
    } | {"path": dest.relative_to(manifest.parent).as_posix(), "sha256": sha256(dest)}


def prepare(args) -> None:
    manifest = args.manifest.resolve()
    out = manifest.parent / (manifest.stem + ".sections")
    out.mkdir(parents=True, exist_ok=True)
    entries = []
    evidence: list[dict] = []
    selection = []
    if args.map:
        accession = args.id
        source = args.source
        if not accession or not source or not source.startswith("https://"):
            raise ValueError("--map requires --id and an official --source URL")
        frame = pd.read_csv(args.map, sep="\t", dtype=str, keep_default_na=False)
        if frame.section_id.duplicated().any():
            raise ValueError("Duplicate section identifiers in mapping")
        for row in frame.to_dict(orient="records"):
            expr = args.map.parent / row["matrix"]
            pos = args.map.parent / row["positions"]
            x, xy, depth = read_panel(expr, pos.read_bytes(), PANEL, with_depth=True)
            entries.append(package_section(out, manifest, row, x, xy, PANEL, depth))
            evidence.extend(
                {"name": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
                for p in [expr, pos]
            )
        description = "Mapped independent patient contrast; declared input is RNA counts."
    else:
        accession = "GSE202011"
        source = GEO + accession
        cache = args.cache
        names = ["GSE202011_samples.soft", "GSE202011_RAW.tar", "GSE202011_st_images.tar.gz"]
        urls = [source + "&targ=gsm&form=text&view=full"] + [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202011/suppl/" + n
            for n in names[1:]
        ]
        from spatial_audit.schema import strict_read

        cat = strict_read(ROOT / "docs/data/catalog.json")
        study = next(c for c in cat["cohorts"] if c["id"] == accession)
        for name, url in zip(names, urls, strict=True):
            asset = next(
                a for a in study["assets"] if a["name"] == name and a["role"] == "audit_input"
            )
            p = cache / name
            digest = sha256(p)
            if digest != asset["sha256"] or p.stat().st_size != asset["bytes"]:
                raise ValueError(f"Archived release asset changed: {name}")
            evidence.append({"name": name, "url": url, "bytes": p.stat().st_size, "sha256": digest})
        selection = geo_records((cache / names[0]).read_text(encoding="utf-8"))
        chosen = [r for r in selection if r["selected"]]
        if {r["group"] for r in chosen} != {"control", "Psoriasis (PsO)"}:
            raise ValueError("Both prespecified groups must be present")
        position_bytes = {}
        with tarfile.open(cache / names[2], mode="r|gz") as images:
            for member in images:
                path = PurePosixPath(member.name)
                if (
                    path.name
                    in {
                        "tissue_positions_list.csv.gz",
                        "tissue_positions_list.csv",
                        "tissue_positions.csv.gz",
                        "tissue_positions.csv",
                    }
                    and member.isfile()
                ):
                    data = images.extractfile(member)
                    if data is None or path.parent.name in position_bytes:
                        raise ValueError("Repeated or unreadable spot mapping")
                    position_bytes[path.parent.name] = data.read()
        with tarfile.open(cache / names[1], mode="r:") as archive:
            for row in chosen:
                name = row["matrix_name"]
                section_key = name.removeprefix(row["section_id"] + "_").removesuffix(
                    "_filtered_feature_bc_matrix.h5"
                )
                positions = position_bytes[section_key + "_spatial"]
                stream = archive.extractfile(name)
                if stream is None:
                    raise ValueError(f"Missing count matrix for {name}")
                x, xy, depth = read_panel(io.BytesIO(stream.read()), positions, PANEL, with_depth=True)
                entries.append(package_section(out, manifest, row, x, xy, PANEL, depth))
        description = (
            "First deposited section per donor within PsO lesional and healthy skin. "
            "Exclude PsA, non-lesional tissue and additional sections from the same donor."
        )
    if len({e["section_id"] for e in entries}) != len(entries) or not entries:
        raise ValueError("Empty or duplicate section set")
    spec = {
        "id": accession,
        "kind": "public",
        "source": source,
        "platform": "Visium",
        "measurement": "counts",
        "features": PANEL,
        "units": "um",
        "files": entries,
        "selection": description
        + " Predefined skin, immune, stromal and housekeeping marker panel. Within-array coordinates "
        "use the Visium v1/v2 nominal physical pitch. Acquisition stage and composition unavailable.",
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        args.provenance,
        {
            "kind": "public",
            "accession": accession,
            "files": evidence,
            "source": source,
            "physical_calibration": GEOMETRY,
            "panel": PANEL,
            "selection": description,
            "screening": selection,
            "sections": len(entries),
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--cache", type=Path)
    modes.add_argument("--map", type=Path)
    parser.add_argument("--id")
    parser.add_argument("--source")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    prepare(parser.parse_args())
