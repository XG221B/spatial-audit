"""Prepare the prespecified prolonged versus pneumonia-control arm of GSE190732."""

from __future__ import annotations

import argparse
import re
import tarfile
from pathlib import Path, PurePosixPath

import pandas as pd
from prepare_visium import GEO, ROOT, package_section, parse_spot_map, read_panel

from spatial_audit.schema import sha256, strict_read, write_json

ACCESSION = "GSE190732"
GEOMETRY = "https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf"
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
    "KRT8",
    "KRT18",
    "KRT19",
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "VIM",
    "PECAM1",
    "VWF",
    "ACTA2",
    "PTPRC",
    "CD3D",
    "CD68",
    "LYZ",
    "HLA-DRA",
]


def shared_reference(paths: list[Path]) -> pd.DataFrame:
    maps = []
    if len(paths) != 2 or len(set(paths)) != 2:
        raise ValueError("Two independent Visium v1 reference datasets are required")
    for path in paths:
        with tarfile.open(path, "r:gz") as archive:
            f = archive.extractfile("spatial/tissue_positions_list.csv")
            if f is None:
                raise ValueError("Missing reference capture map")
            maps.append(parse_spot_map(f.read())[["row", "col"]].sort_index())
    if not maps[0].equals(maps[1]):
        raise ValueError("Manufacturer datasets disagree on barcode-to-address mapping")
    return maps[0]


def screen(soft: str) -> list[dict]:
    rows = []
    for part in re.split(r"(?m)^\^SAMPLE = ", soft)[1:]:
        gsm = part.splitlines()[0].strip()
        title = " ".join(re.findall(r"(?m)^!Sample_title = (.+)$", part)[0].split())
        group = re.findall(r"(?m)^!Sample_characteristics_ch1 = disease group: (.+)$", part)[0]
        filenames = re.findall(r"(?m)^!Sample_supplementary_file(?:_\d+)? = (ftp://.+\.h5)$", part)
        if (
            len(filenames) != 1
            or not re.fullmatch(r"(acute|chronic|prolonged|control) case \d+", title)
            or "Visium Spatial Gene Expression 3" not in part
        ):
            raise ValueError(f"Unexpected Visium v1 donor or assay record for {gsm}")
        arm = title.split()[0]
        expected = "non-COVID-19-related pneumonia" if arm == "control" else arm + " COVID-19"
        if group != expected:
            raise ValueError("GEO case and disease labels disagree")
        rows.append(
            {
                "section_id": gsm,
                "patient_id": title,
                "group": group,
                "selected": arm in {"control", "prolonged"},
                "file": PurePosixPath(filenames[0]).name,
            }
        )
    if len({r["patient_id"] for r in rows}) != len(rows):
        raise ValueError("Repeated clinical case in the deposited sample sheet")
    chosen = pd.Series([r["group"] for r in rows if r["selected"]]).value_counts().to_dict()
    if chosen != {"non-COVID-19-related pneumonia": 3, "prolonged COVID-19": 3}:
        raise ValueError("Case inventory differs from the documented spatial cohort")
    return rows


def prepare(cache: Path, manifest: Path, provenance: Path):
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    study = next(s for s in catalog["cohorts"] if s["id"] == ACCESSION)
    assets = [a for a in study["assets"] if a["role"] == "audit_input"]
    for a in assets:
        path = cache / a["name"]
        if path.stat().st_size != a["bytes"] or sha256(path) != a["sha256"]:
            raise ValueError(f"Archived input changed: {path.name}")
    lookup = {a["name"] for a in assets}
    reference = shared_reference(
        [cache / a["name"] for a in assets if a["name"].endswith("_spatial.tar.gz")]
    )
    rows = screen((cache / "GSE190732_samples.soft").read_text(encoding="utf8"))
    out = manifest.parent / (manifest.stem + ".sections")
    out.mkdir(parents=True, exist_ok=True)
    sections = []
    for row in rows:
        if row["selected"]:
            if row["file"] not in lookup:
                raise ValueError("GEO matrix has no verified source receipt")
            x, xy, depth = read_panel(cache / row["file"], reference, PANEL,
                                      chemistry="Spatial 3' v1", with_depth=True)
            sections.append(package_section(out, manifest, row, x, xy, PANEL, depth))
            row["matched_barcodes"] = len(x)
    description = (
        "Prolonged COVID-19 versus non-COVID pneumonia control cases in the deposited "
        "autopsy lung study. Predefined lung epithelial, vascular, immune, matrix and "
        "housekeeping panel. Filtered matrix barcodes joined to canonical Visium v1 "
        "addresses validated against two manufacturer datasets. Acquisition stage and "
        "cell composition unavailable."
    )
    spec = {
        "id": ACCESSION,
        "kind": "public",
        "platform": "Visium",
        "units": "um",
        "measurement": "counts",
        "source": GEO + ACCESSION,
        "features": PANEL,
        "selection": description,
        "files": sections,
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        provenance,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": GEO + ACCESSION,
            "paper": study["paper"],
            "donor_source": "Methods, Spatial transcriptomics; Fig. 1, each column represents paired consecutive sections of one donor.",
            "selection": description,
            "screening": rows,
            "physical_calibration": GEOMETRY,
            "reference_addresses": len(reference),
            "chemistry_verified": "Spatial 3' v1",
            "reference_fields_used": ["barcode", "array_row", "array_col"],
            "panel": PANEL,
            "files": assets,
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--provenance", type=Path, required=True)
    args = p.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
