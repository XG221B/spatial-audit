"""Reconstruct the frozen Visium analysis panel from two official GEO supplementary files."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import tarfile
from contextlib import nullcontext
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.archive import verified_file
from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE220442"
METADATA = ACCESSION + "_metadata.csv.gz"
ARCHIVE = ACCESSION + "_counts_and_images.tar.gz"
PANEL = ROOT / "data" / ACCESSION / "feature_panel.tsv"
COVARIATES = [f"ct_Layer_{i}" for i in range(2, 7)] + ["ct_White_Matter"]
MEMBER = re.compile(
    r"counts_and_images/([A-Za-z0-9_-]+)/(filtered_feature_bc_matrix\.h5|spatial/tissue_positions_list\.csv|spatial/scalefactors_json\.json)(\.gz)?"
)


def verify_assets(raw, download):
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    row = next(r for r in catalog["cohorts"] if r["id"] == ACCESSION)
    urls = {a["name"]: a["url"] for a in row["assets"] if a["role"] == "data"}
    raw.mkdir(parents=True, exist_ok=True)
    records = []
    for pin in strict_read(PANEL.parent / "archive_fingerprints.json")["files"]:
        name = pin["name"]
        file = raw / name
        verified_file(file, urls[name], pin, download)
        records.append({**pin, "url": urls[name], "reference_digest_matched": True})
    return records


def selected_members(archive, directory):
    samples = {}
    with tarfile.open(archive, "r|gz") as tar:
        for item in tar:
            match = MEMBER.fullmatch(item.name)
            if not match or not item.isfile():
                continue
            sample, member, compressed = match.groups()
            file = directory / sample / Path(member).name
            if file.exists():
                raise ValueError(f"Duplicate archive member {item.name}")
            file.parent.mkdir(exist_ok=True)
            stream = tar.extractfile(item)
            if stream is None:
                raise ValueError(f"Unreadable archive member {item.name}")
            with (
                stream,
                gzip.GzipFile(fileobj=stream) if compressed else nullcontext(stream) as source,
                file.open("wb") as target,
            ):
                shutil.copyfileobj(source, target, 1 << 20)
            samples[sample] = file.parent
    return samples


def make_section(folder, patient, metadata, panel):
    """Match by donor, spot barcode and Ensembl ID before normalization or geometry."""
    m = metadata.loc[metadata["patientID"] == patient].copy()
    labels = m["category"].unique().tolist()
    if len(labels) != 1 or labels[0] not in {"AD", "control"}:
        raise ValueError(f"Ambiguous archived disease category for donor {patient}")
    barcode = m.iloc[:, 0].astype(str).str.rsplit("-", n=1, expand=True)
    if not barcode[1].eq(patient).all() or barcode[0].duplicated().any():
        raise ValueError(f"Nonunique donor/barcode map for {patient}")
    m.index = barcode[0].to_numpy()
    with h5py.File(folder / "filtered_feature_bc_matrix.h5", "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(map(int, matrix["shape"][()]))
        counts = sparse.csr_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]), shape=shape[::-1]
        ).astype(np.int32)
        var = pd.Index(matrix["features/id"].asstr()[:])
        bc = pd.Index(matrix["barcodes"].asstr()[:], name="barcode")
    if (
        not var.is_unique
        or not bc.is_unique
        or not panel.index.is_unique
        or not panel.ensembl_id.is_unique
    ):
        raise ValueError("Nonunique expression or panel identifiers")
    fids = var.get_indexer(panel.ensembl_id)
    if (fids < 0).any():
        raise ValueError(f"Selected Ensembl IDs missing in {patient}")
    pos = pd.read_csv(
        folder / "tissue_positions_list.csv",
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "py", "px"],
    ).set_index("barcode")
    if not pos.index.is_unique:
        raise ValueError("Duplicate position barcode")
    pos = pos.loc[bc]
    keep = pos["in_tissue"].eq(1).to_numpy()
    pos = pos.iloc[keep]
    bc = bc[keep]
    names = bc.str.rsplit("-", n=1).str[0]
    if set(names) != set(m.index):
        raise ValueError(f"Filtered barcodes disagree with source metadata for {patient}")
    annotation = m.loc[names, "Layer"].to_numpy()
    if pd.isna(annotation).any():
        raise ValueError("Missing tissue layer annotation")
    raw = counts[keep][:, fids].tocsr()
    lib = np.asarray(counts[keep].sum(1)).ravel().astype(float)
    normalized = (sparse.diags(1e4 / np.maximum(lib, 1)) @ raw).tocsr().astype(np.float32)
    normalized.data = np.log1p(normalized.data)
    diameter = float(
        json.loads((folder / "scalefactors_json.json").read_text(encoding="utf-8"))[
            "spot_diameter_fullres"
        ]
    )
    if not np.isfinite(diameter) or diameter <= 0:
        raise ValueError("Missing positive spot diameter from Space Ranger scale factors")
    pixel = pos[["px", "py"]].to_numpy(dtype=float)
    rc = pos[["array_row", "array_col"]].to_numpy(dtype=float)
    if (
        not np.isfinite(rc).all() or not np.isfinite(pixel).all()
        or np.any(rc != np.floor(rc)) or np.any(rc < 0)
        or np.any(rc[:, 0] >= 78) or np.any(rc[:, 1] >= 128)
        or np.any((rc.sum(axis=1) % 2) != 0) or len(np.unique(rc, axis=0)) != len(rc)
    ):
        raise ValueError("Invalid 6.5-mm Visium array addresses or image coordinates")
    sites = {tuple(address): point for address, point in zip(rc, pixel, strict=True)}
    distances = [
        np.linalg.norm(sites[other] - point)
        for (row, col), point in sites.items()
        for other in [(row, col + 2), (row + 1, col + 1), (row + 1, col - 1)]
        if other in sites
    ]
    if not distances or min(distances) <= 0:
        raise ValueError("No valid image distance between adjacent array addresses")
    pixel_pitch = float(np.median(distances))
    # Array address spacing supplies physical coordinates independent of H&E magnification.
    array_xy = np.column_stack((50 * rc[:, 1], 50 * np.sqrt(3.0) * rc[:, 0]))
    section = ACCESSION + "_" + patient
    obs = pd.DataFrame(
        {
            "patient_id": patient,
            "section_id": section,
            "condition": labels[0],
            "layer_clean": annotation,
            "total_counts": lib,
        },
        index=bc,
    )
    for name in [
        "Layer 1",
        "Layer 2",
        "Layer 3",
        "Layer 4",
        "Layer 5",
        "Layer 6",
        "White Matter",
        "Noise",
    ]:
        obs["ct_" + name.replace(" ", "_")] = (annotation == name).astype(np.float32)
    obj = ad.AnnData(X=normalized, obs=obs, var=panel.copy())
    obj.layers["counts"] = raw
    obj.obsm["spatial"] = array_xy
    obj.obsm["spatial_pixels"] = pixel
    obj.uns["unit"] = "um"
    obj.uns["spot_diameter_fullres"] = diameter
    return obj, {
        "patient_id": patient,
        "section_id": section,
        "archive_sample": folder.name,
        "group": labels[0],
        "spots": len(bc),
        "nominal_array_pitch_um": 100.0,
        "median_adjacent_spacing_pixels": pixel_pitch,
        "legacy_image_derived_pitch_um": pixel_pitch * (55.0 / diameter),
        "library_denominator": "all filtered-matrix genes before panel subset",
    }


def check_reference(obj, path):
    ref = ad.read_h5ad(path)
    compare = {
        "barcodes_equal": obj.obs_names.equals(ref.obs_names),
        "features_equal": obj.var_names.equals(ref.var_names),
        "counts_equal": (obj.layers["counts"] != ref.layers["counts"]).nnz == 0,
        "normalized_X_equal": (obj.X != ref.X).nnz == 0,
        "legacy_image_transform_reproduced": np.array_equal(
            obj.obsm["spatial_pixels"] * (55.0 / obj.uns["spot_diameter_fullres"]),
            ref.obsm["spatial"],
        ),
        "patient_group_layer_equal": all(
            np.array_equal(obj.obs[k].astype(str), ref.obs[k].astype(str))
            for k in ("patient_id", "condition", "layer_clean")
        ),
        "audit_covariates_equal": np.array_equal(
            np.asarray(obj.obs[COVARIATES]), np.asarray(ref.obs[COVARIATES])
        ),
    }
    if not all(compare.values()):
        raise ValueError(f"Scientific matrices differ for {path.name}: {compare}")
    return compare


def prepare(raw, out, download=False, reference=None, evidence=None):
    out.mkdir(parents=True, exist_ok=True)
    files = verify_assets(raw, download)
    meta = pd.read_csv(raw / METADATA, dtype={"patientID": str})
    panel = pd.read_csv(PANEL, sep="\t").set_index("symbol")
    sections = []
    inputs = []
    with TemporaryDirectory(prefix="spatial-audit-visium-", dir=out) as temp:
        groups = selected_members(raw / ARCHIVE, Path(temp))
        names = {k.replace("-", "_"): v for k, v in groups.items()}
        if len(groups) != len(names) or set(names) != set(meta["patientID"]):
            raise ValueError("Archive directories and deposited donor IDs disagree")
        for donor, folder in sorted(names.items(), key=lambda pair: pair[0]):
            obj, row = make_section(folder, donor, meta, panel)
            name = row["section_id"] + ".h5ad"
            if reference:
                row["reference_comparison"] = check_reference(obj, reference / name)
            pending = Path(temp) / name
            obj.write_h5ad(pending, compression="gzip")
            digest = sha256(pending)
            product = f"{row['section_id']}-{digest[:16]}.h5ad"
            pending.replace(out / product)
            inputs.append({"path": product, "sha256": digest})
            row["prepared_file"] = product
            row["prepared_sha256"] = digest
            sections.append(row)
    # Preserve the analysis feature budget used by the shipped public report.
    schema = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "measurement": "counts",
        "preprocessing_evidence": evidence.name if evidence else "dataset_manifest.json",
        "kind": "public",
        "platform": "Visium",
        "units": "um",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442",
        "layer": "counts",
        "feature_limit": 32,
        "selection": "First 32 genes in the frozen, ordered Visium benchmark panel (six donors).",
        "domain": "layer_clean",
        "library_size_column": "total_counts",
        "columns": {"group": "condition"},
        "covariates": COVARIATES,
        "files": inputs,
    }
    write_json(out / "manifest.json", {"cohorts": [schema]})
    full = {
        **schema,
        "id": ACCESSION + "-full-panel",
        "feature_limit": len(panel),
        "matrix_mode": "memory",
        "parent_accession": ACCESSION,
        "selection": "Every feature in the frozen Visium analysis panel, all six archived donors.",
    }
    write_json(out / "manifest-full-panel.json", {"cohorts": [full]})
    record = {
        "kind": "public",
        "record_type": "archive_reconstruction",
        "accession": ACCESSION,
        "coordinate_calibration": "GEO archive barcode-matched Visium array addresses at nominal 100 um nearest-neighbor spacing. Legacy H5AD image scaling is checked separately from the physical lattice.",
        "input_files": files,
        "feature_panel_sha256": sha256(PANEL),
        "panel_size": len(panel),
        "audit_feature_limit": schema["feature_limit"],
        "n_patients": len(sections),
        "n_observations": sum(r["spots"] for r in sections),
        "sections": sections,
        "mapping": "Direct patientID and category columns of the archived metadata CSV; "
        "per-spot donor barcode suffix checked against each archive directory.",
        "covariates": "Cortical layer and white-matter indicator columns from archived "
        "Layer labels. They describe anatomy; no cell-type fractions were estimated.",
    }
    write_json(out / "dataset_manifest.json", record)
    if evidence:
        write_json(evidence, record)
    print(f"Prepared {len(sections)} Visium sections in {out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--raw-dir",
        type=Path,
        required=True,
        help="Directory for the two pinned GEO supplementary files",
    )
    p.add_argument("--out", type=Path, required=True, help="Output H5ADs and relative manifest")
    p.add_argument(
        "--download",
        action="store_true",
        help="Fetch any missing pinned asset from its official URL",
    )
    p.add_argument(
        "--reference",
        type=Path,
        help="Optional existing GSE220442_ready directory for array comparison",
    )
    p.add_argument(
        "--evidence", type=Path, help="Optional additional copy of portable reconstruction evidence"
    )
    args = p.parse_args()
    prepare(args.raw_dir, args.out, args.download, args.reference, args.evidence)
