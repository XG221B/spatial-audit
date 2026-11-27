"""Rebuild the cortical AIN versus ATI Xenium biopsy contrast in GSE336890."""

import argparse
import gzip
import re
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from fetch_and_verify import fetch
from scipy import sparse

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE336890"
TARGETS = [
    "AQP1",
    "LRP2",
    "CUBN",
    "SLC5A2",
    "SLC34A1",
    "UMOD",
    "SLC12A1",
    "SLC12A3",
    "AQP2",
    "HAVCR1",
    "VCAM1",
    "VIM",
    "COL1A1",
    "COL3A1",
    "DCN",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD68",
    "HLA-DRA",
    "LYZ",
    "PECAM1",
    "VWF",
    "EPCAM",
    "ACTB",
    "GAPDH",
    "B2M",
    "PPIA",
    "RPLP0",
]


def prepare(cache, manifest, provenance):
    row = next(
        c for c in strict_read(ROOT / "docs/data/catalog.json")["cohorts"] if c["id"] == ACCESSION
    )

    def get(f):
        return fetch(
            {"files": [dict(f, accession=ACCESSION)]}, cache, local=(cache / f["name"]).exists()
        )[0]

    receipts = [
        get({"name": "samples.soft", "url": row["source"] + "&targ=gsm&form=text&view=full"})
    ]
    design = next(a for a in row["assets"] if "cells_stats" in a["name"])
    receipts.append(get(design))
    files = {}
    text = (cache / "samples.soft").read_text(encoding="utf8")
    for sample in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        r = re.search(r"!Sample_title = Xenium human kidney biopsy, Region(\d+)", sample)
        if not r:
            raise ValueError("Unexpected GEO region title")
        urls = [
            u.replace("ftp:", "https:")
            for u in re.findall(r"!Sample_supplementary_file(?:_\d+)? = (\S+)", sample)
        ]
        files[int(r[1])] = (sample.splitlines()[0], sample, urls)
    sections = []
    evidence = []
    panel = None
    donors: dict[str, str] = {}
    dest = manifest.parent / "cortical-biopsies"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(cache / design["name"]) as archive:
        regions: dict[int, list[tuple[str,str,str,tarfile.TarInfo]]] = {}
        for member in archive.getmembers():
            m = re.fullmatch(
                r"cells_stats_dir/Region_(\d+)/((AIN|ATI|REF)_([CM])_(\d+)_([A-Za-z0-9]+))_cells_stats\.csv",
                member.name,
            )
            if not m:
                continue
            region, label, group, zone, _number, pid = m.groups()
            if donors.setdefault(pid, group) != group:
                raise ValueError("Biopsy identifier crosses diagnostic groups")
            if zone == "C" and group in ("AIN", "ATI"):
                regions.setdefault(int(region), []).append((label, pid, group, member))
        if sorted(list(donors.values()).count(g) for g in set(donors.values())) != [5, 7, 8]:
            raise ValueError(
                "Xenium tissue inventory differs from the published AIN/ATI/reference design"
            )
        for region, specimens in sorted(regions.items()):
            gsm, description, urls = files[region]
            paths = []
            for suffix in ("_cell_feature_matrix.h5", "_cells.parquet.gz"):
                matched = [u for u in urls if u.endswith(suffix)]
                if len(matched) != 1:
                    raise ValueError("Expected one count matrix and one centroid file")
                f = get({"name": matched[0].split("/")[-1], "url": matched[0]})
                receipts.append(f)
                paths.append(cache / f["name"])
            with gzip.open(paths[1], "rb") as buf:
                cells = (
                    pq.read_table(
                        buf, columns=["cell_id", "x_centroid", "y_centroid", "transcript_counts", "total_counts"]
                    )
                    .to_pandas()
                    .set_index("cell_id")
                )
            with h5py.File(paths[0]) as h:
                matrix = h["matrix"]
                barcode = pd.Index(matrix["barcodes"].asstr()[:])
                if (
                    not barcode.is_unique
                    or not cells.index.is_unique
                    or len(barcode.difference(cells.index))
                ):
                    raise ValueError(
                        "Cannot match unique cellular barcodes across archived products"
                    )
                cells = cells.loc[barcode]
                types = matrix["features/feature_type"].asstr()[:]
                names = matrix["features/name"].asstr()[:]
                ids = []
                for gene in TARGETS:
                    match = np.flatnonzero((types == "Gene Expression") & (names == gene))
                    if len(match) == 1:
                        ids.append(match[0])
                this_panel = names[ids].tolist()
                if panel is None:
                    panel = this_panel
                if panel != this_panel or not panel:
                    raise ValueError("Target families disagree between regions")
                counts = sparse.csc_matrix(
                    (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                    shape=tuple(matrix["shape"][:]),
                )
                total = np.asarray(counts.sum(axis=0)).ravel()
                library = cells.transcript_counts.to_numpy()
                if not np.array_equal(total, cells.total_counts) or np.any(library > total):
                    raise ValueError("Complete Xenium matrices and centroid-file totals disagree")
                X = counts[ids, :].T.toarray()
            used: set[str] = set()
            for label, pid, group, member in specimens:
                if label not in description:
                    raise ValueError("Biopsy label missing from GEO region description")
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError("Cell membership table is missing")
                with stream:
                    barcodes = pd.Index(pd.read_csv(stream)["CellID"])
                if not barcodes.is_unique or used.intersection(barcodes):
                    raise ValueError("Overlapping cell assignments")
                used.update(barcodes)
                ix = barcode.get_indexer(barcodes)
                if np.any(ix < 0):
                    raise ValueError("Unmatched assigned cells")
                xy = cells.iloc[ix][["x_centroid", "y_centroid"]].to_numpy()
                path = dest / (label + ".npz")
                np.savez_compressed(
                    path, coords=xy, X=X[ix], features=np.asarray(panel), library_size=library[ix]
                )
                sections.append(
                    {
                        "path": path.relative_to(manifest.parent).as_posix(),
                        "sha256": sha256(path),
                        "patient_id": pid,
                        "group": group,
                        "section_id": label,
                        "batch": gsm,
                    }
                )
                evidence.append(
                    {
                        "gsm": gsm,
                        "biopsy": pid,
                        "section": label,
                        "group": group,
                        "cells": len(ix),
                        "cell_map_member": member.name,
                    }
                )
    selection = "All cortical biopsies assigned to AIN or ATI in the published Xenium tissue subset. Per-cell membership comes from the submitter's region/section CSV; raw feature counts and micrometre Xenium centroids are joined by cell ID. The batch field identifies the acquisition region."
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": ACCESSION,
                    "accession": ACCESSION,
                    "kind": "public",
                    "platform": "Xenium",
                    "source": row["source"],
                    "units": "um",
                    "measurement": "counts",
                    "batch_label": "GEO Xenium acquisition region",
                    "selection": selection,
                    "files": sections,
                }
            ]
        },
    )
    write_json(
        provenance,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "accession": ACCESSION,
            "source": row["source"],
            "paper": row["paper"],
            "selection": selection,
            "features": panel,
            "files": receipts,
            "published_unit_check": "Paper reports eight AIN, seven ATI and five reference Xenium tissues. Unique biopsy codes reproduce those group counts before restriction to cortex.",
            "sections": evidence,
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("cache", "manifest", "provenance"):
        p.add_argument("--" + name, type=Path, required=True)
    args = p.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
