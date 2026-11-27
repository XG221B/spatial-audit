"""Prepare the baseline Xenium UC/healthy contrast from Dataset 1 of GSE282123."""

import argparse
import os
import re
import tempfile
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from fetch_and_verify import fetch
from scipy import sparse

from spatial_audit.schema import sha256, write_json

ACCESSION = "GSE282123"
SOURCE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + ACCESSION
FIGSHARE = "https://doi.org/10.6084/m9.figshare.27327813"
PAPER = "https://doi.org/10.1172/JCI202488"
TABLE = "https://dm5migu4zj3pb.cloudfront.net/manuscripts/202000/202488/JCI202488.sdt1-11.xlsx"
CODE = (
    "https://github.com/UCSF-DSCOLAB/spatial_transcriptomics_colitis_analysis_2024/blob/"
    "ef610b37ed40f37a6549465e33cb95ac706b34d1/Preprocessing/Preprocessing_XeniumRep1and2.ipynb"
)
RAW = "25_11_12_Xenium_Dataset1_290_IntReps1and2_Raw.h5ad"
ANNOTATED = "25_11_22_Xenium_Dataset1_290_IntReps1and2_Annotated.h5ad"
CLINICAL = "jci-136-202488-s339.xlsx"
PLAN = [
    {
        "name": RAW,
        "url": "https://ndownloader.figshare.com/files/59547476",
        "bytes": 442568032,
        "archive_checksum": "md5:bc1a7f5bd66be6e1350de89a7081744a",
    },
    {
        "name": ANNOTATED,
        "url": "https://ndownloader.figshare.com/files/59799746",
        "bytes": 1523664984,
        "archive_checksum": "md5:419cab98a3be238d310db6e3fc92bc42",
    },
    {
        "name": CLINICAL,
        "url": TABLE,
        "bytes": 440862,
        "sha256": "fb95db639a496ebcde7507fa5173466d7d835039cdf94565452970bdfe6a1767",
    },
]


def clinical_groups(path):
    """Read covariate-blind donor assignments from published Supplemental Table 1."""
    with ZipFile(path) as book:
        workbook = ElementTree.fromstring(book.read("xl/workbook.xml"))
        first_sheet = workbook.find(".//{*}sheet")
        if first_sheet is None or first_sheet.get("name") != "ST1":
            raise ValueError("Unexpected order in the pinned supplementary workbook")
        strings = [
            "".join(s.itertext()) for s in ElementTree.fromstring(book.read("xl/sharedStrings.xml"))
        ]
        sheet = ElementTree.fromstring(book.read("xl/worksheets/sheet1.xml"))
    result = {}
    for row in sheet.findall(".//{*}sheetData/{*}row"):
        values = {}
        for c in row.findall("{*}c"):
            value = c.findtext("{*}v", "")
            if c.get("t") == "s":
                value = strings[int(value)]
            values[re.sub(r"\d+", "", c.attrib["r"])] = value
        donor = values.get("A", "")
        if re.fullmatch(r"HS\d+", donor):
            group = values["D"]
            if group not in {"HC", "UC"} or donor in result:
                raise ValueError("Ambiguous donor or disease status in Supplemental Table 1")
            result[donor] = group
    if not result:
        raise ValueError("No Dataset 1 donors found")
    return result


def baseline_patient(label, clinical):
    if label == "unassigned" or "_POST_" in label:
        return None
    m = re.fullmatch(r"(HS\d+)_+(HC|PRE_VDZ_(?:NR|R))", label)
    if not m:
        raise ValueError(f"Unrecognized annotated biopsy label: {label}")
    group = "HC" if m[2] == "HC" else "UC"
    if clinical.get(m[1]) != group:
        raise ValueError(f"Disease assignment disagrees with Supplemental Table 1: {label}")
    return m[1], group


def acquisition_coords(coords, batch):
    if batch not in {"Sept", "Nov"}:
        raise ValueError("Unrecognized acquisition label")
    original = np.array(coords, dtype=float, copy=True)
    if batch == "Nov":
        original[:, 0] -= 15_000
    if not np.isfinite(original).all():
        raise ValueError("Nonfinite original Xenium centroid")
    return original


def relative_file(path, base):
    try:
        return Path(os.path.relpath(path, base)).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def write_section(target, **arrays):
    with tempfile.NamedTemporaryFile(dir=target.parent, suffix=".npz", delete=False) as tmp:
        name = Path(tmp.name)
    try:
        np.savez_compressed(name, **arrays)
        name.replace(target)
    finally:
        name.unlink(missing_ok=True)


def matrix_rows(node, start, end):
    if node.attrs["encoding-type"] != "csr_matrix":
        raise ValueError("Pinned count matrix must use CSR encoding")
    ptr = node["indptr"][start : end + 1]
    first, last = int(ptr[0]), int(ptr[-1])
    return sparse.csr_matrix(
        (node["data"][first:last], node["indices"][first:last], ptr - first),
        shape=(end - start, int(node.attrs["shape"][1])),
    )


def project(raw_path, ann_path, clinical, out, manifest_dir):
    with h5py.File(raw_path) as raw, h5py.File(ann_path) as ann:
        obs = ann["obs"]
        sample_col = "24_01_17_HS"
        core_col = "Patient_ID_cores_combined"
        type_col = "25_06_11_Common_Coarse_Xenium_Combinedv3"
        frame = pd.DataFrame(
            {
                k: read_elem(obs[k])
                for k in (sample_col, core_col, type_col, "24_01_17_Condition", "batch")
            }
        )
        batch_labels = frame["batch"].astype(str)
        barcode = pd.Series(read_elem(obs["cell_id"])).astype(str)
        # The deposited Nov annotated identifiers already carry the acquisition
        # prefix, while Sept identifiers carry the unprefixed Xenium barcode.
        prefix = barcode.str.extract(r"^(Sept|Nov)-", expand=False)
        if (prefix.notna() & (prefix != batch_labels)).any():
            raise ValueError("Cell identifier acquisition disagrees with batch")
        cell_key = batch_labels + "-" + barcode.str.replace(r"^(Sept|Nov)-", "", regex=True)
        raw_ids = pd.Index(read_elem(raw["obs/_index"]))
        if not raw_ids.is_unique or not cell_key.is_unique:
            raise ValueError("Barcode + acquisition must uniquely identify each cell")
        raw_index = raw_ids.get_indexer(cell_key)
        if (raw_index < 0).any():
            raise ValueError("Annotated barcode missing in raw file")
        raw_genes = pd.Index(read_elem(raw["var/_index"])).astype(str).to_numpy(dtype=str)
        if len(set(raw_genes)) != len(raw_genes) or not np.array_equal(
            raw_genes, pd.Index(read_elem(ann["var/_index"])).astype(str).to_numpy(dtype=str)
        ):
            raise ValueError("Features disagree between count and annotation releases")
        genes = raw_genes[:32]
        xy = np.column_stack(
            [raw["obs/" + field][:][raw_index] for field in ("x_centroid", "y_centroid")]
        )
        ann_xy = np.column_stack([obs[field][:] for field in ("x_centroid", "y_centroid")])
        if not np.array_equal(xy, ann_xy):
            raise ValueError("Centroids disagree between the two author releases")
        counts = np.empty((len(frame), len(genes)), np.int32)
        depths = np.empty(len(frame), np.int64)
        # Decode raw counts once in observation chunks; validate the matched
        # annotated layer and retain full targeted-library sizes independently.
        raw_n = len(raw_ids)
        lookup = np.full(raw_n, -1, np.int64)
        lookup[raw_index] = np.arange(len(frame))
        for start in range(0, raw_n, 4096):
            end = min(start + 4096, raw_n)
            ids = lookup[start:end]
            keep = ids >= 0
            m = matrix_rows(raw["X"], start, end)[keep]
            if (
                not np.isfinite(m.data).all()
                or np.any(m.data < 0)
                or np.any(m.data != np.floor(m.data))
            ):
                raise ValueError("Raw Xenium counts must be nonnegative integers")
            counts[ids[keep]] = m[:, : len(genes)].toarray()
            depths[ids[keep]] = np.asarray(m.sum(axis=1, dtype=np.int64)).ravel()
        for start in range(0, len(frame), 4096):
            end = min(start + 4096, len(frame))
            m = matrix_rows(ann["layers/raw_counts"], start, end)
            if not np.array_equal(
                m[:, : len(genes)].toarray(), counts[start:end]
            ) or not np.array_equal(
                np.asarray(m.sum(axis=1, dtype=np.float64)).ravel(), depths[start:end]
            ):
                raise ValueError("Archived panel or targeted library size disagrees between files")
        distinct = frame[[sample_col, "24_01_17_Condition"]].drop_duplicates()
        assignments = {}
        for label, condition in distinct.itertuples(index=False, name=None):
            if label == "unassigned":
                if condition != label:
                    raise ValueError("Inconsistent unassigned biopsy")
            elif re.sub(r"^HS\d+_+", "", label) != condition:
                raise ValueError("Within-biopsy condition is ambiguous")
            assignments[label] = baseline_patient(label, clinical)
        levels = list(frame[type_col].cat.categories)
        type_codes = frame[type_col].cat.codes.to_numpy()
        if np.any(type_codes < 0):
            raise ValueError("Missing coarse cell class")
        included = frame[sample_col].isin([k for k, v in assignments.items() if v])
        sections = []
        evidence = []
        for (batch, core), cells in frame.loc[included].groupby(
            ["batch", core_col], observed=True, sort=True
        ):
            labels = cells[sample_col].unique()
            if len(labels) != 1:
                raise ValueError("A TMA core maps to multiple biopsies")
            donor, group = assignments[labels[0]]
            if not str(core).startswith(donor + "_") or not re.fullmatch(
                "[A-Za-z0-9_-]+", str(core)
            ):
                raise ValueError("Core and biopsy patient identifiers disagree")
            sid = f"{batch}_{core}"
            file = out / (sid + ".npz")
            rows = cells.index.to_numpy()
            write_section(
                file,
                X=counts[rows],
                coords=acquisition_coords(xy[rows], batch),
                C=np.eye(len(levels), dtype=np.uint8)[type_codes[rows]],
                features=genes,
                library_size=depths[rows],
            )
            sections.append(
                {
                    "path": relative_file(file, manifest_dir),
                    "patient_id": donor,
                    "group": group,
                    "slide_id": "Dataset1-" + batch,
                    "batch": batch,
                    "section_id": sid,
                    "sha256": sha256(file),
                }
            )
            evidence.append(
                {
                    "patient_id": donor,
                    "group": group,
                    "core": str(core),
                    "run": str(batch),
                    "biopsy": str(labels[0]),
                    "cells": len(rows),
                }
            )
    details = {
        "raw_cells": raw_n,
        "annotated_cells": len(frame),
        "included_cells": int(included.sum()),
        "included_donors": len({s["patient_id"] for s in sections}),
        "excluded_conditions": frame.loc[~included, "24_01_17_Condition"]
        .value_counts()
        .loc[lambda x: x > 0]
        .to_dict(),
        "clinical_table_donors_without_baseline_cells": sorted(
            set(clinical) - {s["patient_id"] for s in sections}
        ),
        "composition_classes": levels,
        "crosscheck": "All annotated cell keys and centroids matched to the raw archive, with exact agreement of selected counts and full targeted-library sums.",
        "sections": evidence,
    }
    return genes.tolist(), sections, details


def prepare(cache, manifest, receipt):
    cache.mkdir(parents=True, exist_ok=True)
    records = []
    for row in PLAN:
        local = cache / row["name"]
        records.extend(
            fetch({"files": [{"accession": ACCESSION, **row}]}, cache, local=local.is_file())
        )
    clinical = clinical_groups(cache / CLINICAL)
    out = cache / "baseline-sections"
    out.mkdir(exist_ok=True)
    genes, sections, details = project(
        cache / RAW, cache / ANNOTATED, clinical, out, manifest.parent
    )
    spec = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "kind": "public",
        "preprocessing_evidence": receipt.name,
        "source": SOURCE,
        "platform": "Xenium",
        "measurement": "counts",
        "units": "um",
        "features": genes,
        "files": sections,
        "covariates": "One-hot archived coarse cell classes, averaged within spatial volumes; expression-derived annotations can share information with the tested panel.",
        "selection": "Dataset 1, HC vs UC pre-vedolizumab; annotated QC-passing cells with assigned biopsies; both technical runs and all eligible cores, grouped by donor; first 32 deposited genes in archive order with full 290-target depth.",
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        receipt,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "accession": ACCESSION,
            "paper": PAPER,
            "figshare": FIGSHARE,
            "clinical": TABLE,
            "coordinate_source": CODE,
            "coordinate_frame": "Xenium centroids in micrometres. The author-added 15,000 um x-offset is removed for Nov. Patient group comparisons use core/run local bounding boxes as acquisition proxies; no registered instrument-stage measurements.",
            "scope": spec["selection"],
            "covariates": spec["covariates"],
            "files": records,
            **details,
        },
    )
    print(
        f"{details['included_cells']:,} baseline cells, {details['included_donors']} donors, {len(sections)} core/run sections"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
