"""Build the GSE250346 patient-level Xenium contrast from corrected GEO data and Supplementary Table 1."""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
import pandas as pd
from project_fibrosis_rds import project

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE250346"
FIELDS = ["sample", "patient", "disease_status", "tma", "run"]
ANNOTATION_CODE = "https://github.com/Banovich-Lab/Spatial_PF/blob/99da8d1954e0c2680f0b87a672ac5eda9f438178/preprocessing_and_celltype_annotation.R"
UNITS = "https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/outputs/xoa-output-understanding-outputs"


def donor_table(path):
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as archive:
        strings = [
            "".join(si.itertext())
            for si in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("x:si", ns)
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet.findall(".//x:row", ns):
            cells = {}
            for cell in row.findall("x:c", ns):
                addr = re.sub("[0-9]", "", cell.attrib["r"])
                value = cell.findtext("x:v", "", ns)
                cells[addr] = strings[int(value)] if cell.get("t") == "s" else value
            if cells.get("C") in {"Control", "Disease"}:
                rows.append(
                    {
                        "patient": cells["A"],
                        "disease_status": cells["C"],
                        "tmas": cells["I"],
                        "runs": cells["J"],
                    }
                )
    t = pd.DataFrame(rows)
    if t.patient.duplicated().any() or t.disease_status.value_counts().to_dict() != {
        "Control": 9,
        "Disease": 26,
    }:
        raise ValueError("Unexpected donor inventory in publication Supplementary Table 1")
    return t.set_index("patient")


def validated_study(cache):
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    study = next(c for c in catalog["cohorts"] if c["id"] == ACCESSION)
    assets = [a for a in study["assets"] if a["role"] == "audit_input"]
    for a in assets:
        p = cache / a["name"]
        if p.stat().st_size != a["bytes"] or sha256(p) != a["sha256"]:
            raise ValueError(f"Source release changed: {p.name}")
    return study, assets


def select_samples(samples: pd.DataFrame, donors: pd.DataFrame) -> pd.DataFrame:
    table = samples.merge(
        donors, left_on="patient", right_index=True, validate="many_to_one", suffixes=("", "_table")
    )
    if (
        set(table.patient) != set(donors.index)
        or len(table) != len(samples)
        or not table.disease_status.equals(table.disease_status_table)
    ):
        raise ValueError("Seurat clinical IDs or status differ from the published donor table")
    for field, summary, prefix in [("tma", "tmas", "TMA"), ("run", "runs", "Run")]:
        match = []
        for code, record in zip(table[field], table[summary], strict=True):
            if re.fullmatch(prefix + r"\d+", code) is None:
                raise ValueError("Malformed core-level acquisition code")
            match.append(code.removeprefix(prefix) in re.findall(r"(?<![.\d])\d+", record))
        table[field + "_in_donor_summary"] = match
    # Supplementary Table 1 summarizes each donor once. Use the recorded core's
    # acquisition fields and retain the per-core agreement flags for review.
    table = table.sort_values(["patient", "sample"]).reset_index(drop=True)
    table["retained"] = ~table.patient.duplicated()
    return table


def package(cache, export, manifest, provenance):
    study, assets = validated_study(cache)
    donors = donor_table(cache / "41588_2025_2080_MOESM5_ESM.xlsx")
    panel = (export / "features.txt").read_text(encoding="utf8").splitlines()
    sections = pd.read_csv(export / "samples.csv")
    if sections["sample"].duplicated().any() or len(panel) != len(set(panel)):
        raise ValueError("Nonunique sample or feature")
    rows = []
    for row in sections.itertuples():
        frame = pd.read_csv(export / row.file, usecols=FIELDS)
        labels = frame.drop_duplicates()
        if len(frame) != row.cells or len(labels) != 1 or labels.iloc[0]["sample"] != row.sample:
            raise ValueError("A core must resolve to one mapped patient, status, TMA and run")
        rows.append({**labels.iloc[0].to_dict(), "file": row.file, "nuclei": row.cells})
    screening = select_samples(pd.DataFrame(rows), donors)
    folder = manifest.parent / (manifest.stem + ".sections")
    folder.mkdir(parents=True, exist_ok=True)
    files = []
    for row in screening.loc[screening.retained].itertuples():
        frame = pd.read_csv(export / row.file)
        xy = frame[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
        x = frame[panel].to_numpy(dtype=np.int32)
        depth = frame.library_size.to_numpy(float)
        if (
            not np.isfinite(xy).all()
            or frame.cell_id.duplicated().any()
            or len(np.unique(xy, axis=0)) != len(xy)
            or np.any(depth < x.sum(axis=1))
        ):
            raise ValueError("Cell identity, raw centroid or depth fails validation")
        c = pd.get_dummies(frame.final_lineage, dtype=np.float32)
        dest = folder / (row.sample + ".npz")
        np.savez_compressed(
            dest, coords=xy, X=x, library_size=depth, features=np.array(panel), C=c.to_numpy()
        )
        files.append(
            {
                "path": dest.relative_to(manifest.parent).as_posix(),
                "sha256": sha256(dest),
                "section_id": row.sample,
                "patient_id": row.patient,
                "group": row.disease_status,
                "slide_id": row.tma,
                "batch": row.run,
            }
        )
    selection = (
        "All mapped publication donors; lexicographically first core in the archived sample field "
        "per donor. First 32 deposited RNA targets in assay order. Native Xenium centroids, within-core "
        "lineage indicators. Original all-target RNA depth verified against the deposited nCount_RNA."
    )
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": ACCESSION,
                    "kind": "public",
                    "platform": "Xenium",
                    "source": study["source"],
                    "units": "um",
                    "measurement": "counts",
                    "selection": selection,
                    "features": panel,
                    "files": files,
                }
            ]
        },
    )
    write_json(
        provenance,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": study["source"],
            "paper": study["paper"],
            "files": assets,
            "screening": screening.to_dict("records"),
            "panel": panel,
            "selection": selection,
            "units_source": UNITS,
            "study_calibration": "Methods, Annotation of histological features. The study applied a 0.2125 micrometre-per-pixel conversion to its image annotations to harmonize them with native cell centroid coordinates.",
            "author_pipeline": ANNOTATION_CODE,
            "covariates": "Archived expression-derived final_lineage indicators.",
            "coordinates": "Native x_centroid and y_centroid within a tissue core, in micrometres. No inter-slide registration.",
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--provenance", type=Path, required=True)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument(
        "--native-rds",
        action="store_true",
        help="Deserialize the full corrected RDS with a high-memory R installation",
    )
    args = p.parse_args()
    cache = args.cache.resolve()
    _, assets = validated_study(cache)
    rds = cache / next(a["name"] for a in assets if a["name"].endswith(".rds"))
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        work = Path(tmp)
        export = work / "export"
        cmd = [args.rscript, "--vanilla", str(ROOT / "scripts/export_fibrosis_panel.R")]
        if args.native_rds:
            cmd.extend([str(rds), str(export)])
        else:
            project(rds, work)
            cmd.extend([str(work / "metadata.rds"), str(export), str(work)])
        subprocess.run(cmd, check=True)
        package(cache, export, args.manifest.resolve(), args.provenance.resolve())


if __name__ == "__main__":
    main()
