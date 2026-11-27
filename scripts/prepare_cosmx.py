"""Reconstruct the public GSE327879 CosMx patient, FOV and physical-coordinate mappings.

Run fetch_and_verify.py --accession GSE327879 --role audit_input before using this adapter.
The expression matrices and FOV origins are from GEO. Donor and cell annotations are
from the exact Git commit cited in the accession catalog.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE327879"
SLIDES = (("GSM9667707", "Ctrl"), ("GSM9667708", "Covid"))
PANEL_SIZE = 32


def physical_pitch(metadata: pd.DataFrame, positions: pd.DataFrame) -> tuple[np.ndarray, float]:
    """Fit each axis to the GEO FOV origins; reject inconsistencies within and across FOVs."""
    df = pd.DataFrame(
        {
            "fov": metadata["fov"],
            "x": metadata["CenterX_global_px"] - metadata["CenterX_local_px"],
            "y": metadata["CenterY_global_px"] + metadata["CenterY_local_px"],
        }
    )
    fov = df.groupby("fov")[["x", "y"]]
    if not np.isfinite(df).all(axis=None) or (fov.max() - fov.min()).to_numpy().max() > 1e-5:
        raise ValueError("Cell offsets are inconsistent within a field of view")
    origins = fov.median()
    stage = positions.set_index("FOV")
    if not stage.index.is_unique:
        raise ValueError("GEO FOV identifier is repeated")
    if set(origins.index) != set(stage.index):
        raise ValueError("FOV identifiers differ between the coordinate tables")
    mm = stage.loc[origins.index, ["X_mm", "Y_mm"]].to_numpy(dtype=float)
    pixel = origins.to_numpy(dtype=float)
    if not np.isfinite(mm).all() or min(np.ptp(pixel, axis=0)) <= 0:
        raise ValueError("Need finite and distinct FOV origins in both axes")
    slopes = np.array([np.polyfit(pixel[:, axis], mm[:, axis], 1) for axis in (0, 1)])
    error = float(np.max(abs(pixel * slopes[:, 0] + slopes[:, 1] - mm)))
    if error > 1e-5 or np.any(slopes[:, 0] <= 0):
        raise ValueError("Archive origins do not support a common positive pixel scale")
    return slopes[:, 0] * 1000.0, error


def match_cells(metadata: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    keys = ["fov", "cell_ID"]
    joined = metadata.merge(counts, on=keys, how="outer", validate="one_to_one", indicator=True)
    if not (joined["_merge"] == "both").all():
        raise ValueError("Donor metadata and expression matrices refer to different cell sets")
    return joined.drop(columns="_merge")


def prepare(cache: Path, manifest: Path, provenance: Path) -> None:
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    study = next(c for c in catalog["cohorts"] if c["id"] == ACCESSION)
    receipts = []
    for a in study["assets"]:
        if a["role"] != "audit_input":
            continue
        path = cache / a["name"]
        if path.stat().st_size != a["bytes"] or sha256(path) != a["sha256"]:
            raise ValueError(f"Published benchmark input has changed: {path.name}")
        receipts.append(a)
    if len(receipts) != 5:
        raise ValueError("Expected the author metadata and two GEO input pairs")

    metadata = pd.read_csv(cache / "Metadata_CosMx.csv")
    keys = [
        "Condition",
        "slide_ID",
        "fov",
        "cell_ID",
        "patient_id",
        "cell_id",
        "annotation",
        "CenterX_global_px",
        "CenterY_global_px",
        "CenterX_local_px",
        "CenterY_local_px",
    ]
    if metadata[keys].isna().any(axis=None) or metadata["cell_id"].duplicated().any():
        raise ValueError("Missing or nonunique archived biological mappings")
    if set(metadata["Condition"]) != {condition for _, condition in SLIDES}:
        raise ValueError("Unexpected donor conditions")
    if metadata.groupby("patient_id").Condition.nunique().max() != 1:
        raise ValueError("Donor is assigned to multiple conditions")
    typenames = sorted(metadata["annotation"].unique())
    files = []
    audit = []
    common: list[str] = []
    out = manifest.parent / (manifest.stem + ".sections")
    out.mkdir(parents=True, exist_ok=True)

    for accession, group in SLIDES:
        s = metadata.loc[metadata["Condition"] == group].copy()
        if s["slide_ID"].nunique() != 1 or s["patient_id"].nunique() != 3:
            raise ValueError("Slide or patient mapping disagrees with the GEO sample record")
        expected = (
            "c_"
            + s["slide_ID"].astype(str)
            + "_"
            + s["fov"].astype(str)
            + "_"
            + s["cell_ID"].astype(str)
        )
        if not (expected == s["cell_id"]).all():
            raise ValueError("Archived cell identifiers disagree with slide, FOV or cell keys")
        position = pd.read_csv(cache / f"{accession}_{group}_fov_positions_file.csv.gz")
        if set(position["Slide"]) != set(s["slide_ID"]):
            raise ValueError("GEO FOV archive and author metadata disagree on slide identity")
        pitch, residual = physical_pitch(s, position)

        expr = cache / f"{accession}_{group}_exprMat_file.csv.gz"
        header = pd.read_csv(expr, nrows=0).columns
        all_targets = [
            f
            for f in header
            if f not in {"X", "fov", "cell_ID"} and not f.startswith(("Unnamed:", "Negative"))
        ]
        features = all_targets[:PANEL_SIZE]
        if len(features) != PANEL_SIZE or (common and features != common):
            raise ValueError("Target panel differs across slides")
        common = features
        panels = []
        for chunk in pd.read_csv(expr, usecols=["cell_ID", "fov", *all_targets], chunksize=512):
            values = chunk[all_targets].to_numpy(dtype=float)
            if (
                not np.isfinite(values).all()
                or np.any(values < 0)
                or np.any(values != np.floor(values))
            ):
                raise ValueError("Non-count target in the complete registered CosMx library")
            part = chunk[["cell_ID", "fov", *features]].copy()
            part["library_size"] = values.sum(axis=1)
            panels.append(part)
        joined = match_cells(s, pd.concat(panels, ignore_index=True))
        x = joined[features].to_numpy()
        if not np.isfinite(x).all() or np.any(x < 0) or np.any(x != np.floor(x)):
            raise ValueError("Non-count values in archived RNA panel")

        for field, cells in joined.groupby("fov", sort=True):
            if cells.patient_id.nunique() != 1:
                raise ValueError("An FOV has ambiguous donor labels")
            donor = str(cells["patient_id"].iloc[0])
            sid = f"{accession}_fov{field}"
            coords = cells[["CenterX_local_px", "CenterY_local_px"]].to_numpy(float) * pitch
            labels = pd.Categorical(cells["annotation"], categories=typenames).codes
            composition = np.eye(len(typenames), dtype=np.uint8)[labels]
            dest = out / (sid + ".npz")
            with tempfile.NamedTemporaryFile(dir=out, suffix=".npz", delete=False) as handle:
                tmp = Path(handle.name)
            try:
                np.savez_compressed(
                    tmp,
                    X=cells[features].to_numpy(dtype=np.int32),
                    C=composition,
                    coords=coords,
                    features=np.asarray(features),
                )
                os.replace(tmp, dest)
            finally:
                tmp.unlink(missing_ok=True)
            depth_path = out / (sid + ".depth.npy")
            np.save(depth_path, cells["library_size"].to_numpy(float), allow_pickle=False)
            files.append(
                {
                    "path": dest.relative_to(manifest.parent).as_posix(),
                    "library_size_path": depth_path.relative_to(manifest.parent).as_posix(),
                    "library_size_sha256": sha256(depth_path),
                    "group": group,
                    "patient_id": donor,
                    "section_id": sid,
                    "slide_id": accession,
                    "sha256": sha256(dest),
                }
            )
        audit.append(
            {
                "slide": accession,
                "group": group,
                "patients": sorted(s.patient_id.unique()),
                "cells": len(s),
                "fields": int(s.fov.nunique()),
                "um_per_pixel_xy": pitch,
                "max_affine_residual_mm": residual,
            }
        )

    selection = (
        "All author-mapped cells, one graph per FOV. First 32 archived gene targets "
        "in original column order. Clinical contrast and slide are completely aliased. "
        "C contains the original study's expression-derived cell-type indicators. "
        "Shares estimate in-sample annotation sensitivity. FOV widths limit coarse wavelengths."
    )
    spec = {
        "id": ACCESSION,
        "kind": "public",
        "platform": "CosMx SMI",
        "source": study["source"],
        "units": "um",
        "features": common,
        "measurement": "counts",
        "selection": selection,
        "files": files,
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        provenance,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": study["source"],
            "paper": study["paper"],
            "files": receipts,
            "slides": audit,
            "features": common,
            "covariate_levels": typenames,
            "selection": selection,
            "coordinate_frame": "Within-FOV physical centroids. No between-slide "
            "stage registration assumed. GEO X_mm/Y_mm establish each slide's "
            "length calibration from global-minus-local x and global-plus-local y.",
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
