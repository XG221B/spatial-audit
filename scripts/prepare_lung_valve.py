"""Prepare released Visium lung and aortic-valve donor contrasts from sample-level GEO files."""

import argparse
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from fetch_and_verify import fetch
from prepare_visium import GEOMETRY, package_section, parse_spot_map, positions_um, read_panel
from scipy.io import mmread

from spatial_audit.schema import write_json

PANELS = {
    "GSE318039": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "SFTPA1",
        "SFTPB",
        "SFTPC",
        "AGER",
        "SCGB1A1",
        "MUC5AC",
        "MUC5B",
        "CD3D",
        "CD3E",
        "CD8A",
        "NKG7",
        "GZMB",
        "PRF1",
        "CTSW",
        "PDCD1",
        "CXCL9",
        "CXCL10",
        "LYZ",
        "CD68",
        "PECAM1",
        "COL1A1",
        "DCN",
        "ACTB",
        "B2M",
        "GAPDH",
    ],
    "GSE306288": [
        "COMP",
        "SERPINE1",
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "VIM",
        "ACTA2",
        "RUNX2",
        "SPP1",
        "FN1",
        "VCAN",
        "TAGLN",
        "PECAM1",
        "VWF",
        "CD68",
        "CD163",
        "PTPRC",
        "LYZ",
        "JUN",
        "FOS",
        "GAPDH",
        "ACTB",
        "B2M",
    ],
}


def read_market(matrix_file, feature_file, barcode_file, positions, panel):
    genes = pd.read_csv(feature_file, sep="\t", header=None, dtype=str)
    bc = pd.read_csv(barcode_file, sep="\t", header=None, dtype=str).iloc[:, 0].to_numpy()
    counts = mmread(matrix_file).tocsc()
    if (
        counts.shape != (len(genes), len(bc))
        or np.any(counts.data < 0)
        or np.any(counts.data != np.rint(counts.data))
    ):
        raise ValueError("Expected aligned nonnegative RNA counts")
    chosen = []
    for symbol in panel:
        ix = np.flatnonzero((genes[1] == symbol) & (genes[2] == "Gene Expression"))
        if len(ix) != 1:
            raise ValueError("Ambiguous panel feature " + symbol)
        chosen.append(ix[0])
    frame = parse_spot_map(positions)
    loc = frame.index.get_indexer(bc)
    if np.any(loc < 0) or len(set(bc)) != len(bc):
        raise ValueError("Unmatched or repeated barcodes")
    keep = frame.in_tissue.to_numpy()[loc] == 1
    xy = positions_um(frame, bc[keep])
    lib = np.asarray(counts[genes[2].to_numpy() == "Gene Expression"].sum(axis=0)).ravel()
    return counts[chosen][:, keep].T.toarray(), xy, lib[keep]


def prepare(accession, cache, manifest, evidence):
    source = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + accession
    attr = {
        "url": source + "&targ=gsm&form=text&view=full",
        "name": "gsm.soft",
        "accession": accession,
    }

    def get(asset):
        return fetch({"files": [asset]}, cache, local=(cache / asset["name"]).is_file())[0]

    receipts = [get(attr)]
    rows = []
    for section in re.split(r"(?m)^\^SAMPLE = ", (cache / "gsm.soft").read_text(encoding="utf8"))[
        1:
    ]:
        urls = [
            u.replace("ftp:", "https:")
            for u in re.findall(r"!Sample_supplementary_file(?:_\d+)? = (\S+)", section)
        ]
        position = [
            u
            for u in urls
            if u.endswith(("tissue_positions_list.csv.gz", "tissue_positions.csv.gz"))
        ]
        if not position:
            continue
        attrs = dict(re.findall(r"(?m)^!Sample_characteristics_ch1 = ([^:\n]+): (.+)$", section))
        gsm = section.splitlines()[0]
        title = re.search(r"!Sample_title = (.+)", section)
        if title is None:
            raise ValueError("Missing sample title")
        meta = {"section_id": gsm}
        if accession == "GSE318039":
            donor = re.fullmatch(r"lung, (control\d+|COPD\d+)", title[1])
            if not donor:
                raise ValueError("Missing participant ID")
            meta.update(patient_id=donor[1], group=attrs["phenotype"], batch=attrs["batch"])
            suffixes = ["filtered_feature_bc_matrix.h5"]
        else:
            donor = re.fullmatch(
                r"patient (\d+)[ABV] human aortic valve leaflet tissue, (male|female), diseased",
                title[1],
            )
            if not donor or donor[2] != attrs["Sex"] or attrs["disease state"] != "diseased":
                raise ValueError("Unresolved CAVD participant or disease category")
            meta.update(patient_id=donor[1], group=donor[2])
            suffixes = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
        chosen = []
        for suffix in suffixes:
            file = [u for u in urls if u.endswith(suffix)]
            if len(file) != 1:
                raise ValueError("Missing or repeated archive member: " + suffix)
            chosen += file
        rows.append((meta, chosen + position))
    downloads = [
        {"url": u, "name": u.split("/")[-1], "accession": accession}
        for _, urls in rows
        for u in urls
    ]
    with ThreadPoolExecutor(3) as pool:
        receipts.extend(pool.map(get, downloads))
    sections = []
    selection = []
    out = manifest.parent / (accession + "-sections")
    out.mkdir(parents=True, exist_ok=True)
    for meta, urls in rows:
        paths = [cache / u.split("/")[-1] for u in urls]
        panel = PANELS[accession]
        pos = paths[-1].read_bytes()
        if len(paths) == 2:
            x, xy, library = read_panel(paths[0], pos, panel, with_depth=True)
        else:
            x, xy, library = read_market(paths[0], paths[1], paths[2], pos, panel)
        sections.append(package_section(out, manifest, meta, x, xy, panel, library))
        selection.append({**meta, "observations": len(x), "members": [p.name for p in paths]})
    if len({r["patient_id"] for r in sections}) != ({"GSE318039": 8, "GSE306288": 6}[accession]):
        raise ValueError("Discordant biological-unit inventory")
    note = "All archived Visium capture areas with physical position tables. RNA counts, full RNA library depth and a mechanism-informed fixed feature panel. Repeated sections are nested within donor."
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": accession,
                    "accession": accession,
                    "kind": "public",
                    "platform": "Visium",
                    "units": "um",
                    "source": source,
                    "measurement": "counts",
                    "selection": note,
                    "files": sections,
                }
            ]
        },
    )
    write_json(
        evidence,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "source": source,
            "accession": accession,
            "paper": "https://pubmed.ncbi.nlm.nih.gov/"
            + {"GSE318039": "42351239", "GSE306288": "42206364"}[accession]
            + "/",
            "selection": note,
            "features": PANELS[accession],
            "files": receipts,
            "samples": selection,
            "physical_calibration": GEOMETRY,
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--accession", required=True, choices=sorted(PANELS))
    for name in ("cache", "manifest", "provenance"):
        p.add_argument("--" + name, type=Path, required=True)
    args = p.parse_args()
    prepare(
        args.accession, args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve()
    )
