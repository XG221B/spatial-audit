"""Reconstruct a common Xenium count panel using submitted GEO outcomes and cell centroids."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from spatial_audit.archive import verified_file
from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/GSE348499"
ACCESSION = "GSE348499"
METADATA_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348499&targ=gsm&form=text&view=full"
)
OUTCOMES = {"pathologic complete response (pCR)": "pCR", "residual disease (RD)": "RD"}


def parse_geo(payload):
    sections = []
    for match in re.finditer(r"(?ms)^\^SAMPLE = (GSM\d+)\s*$\n(.*?)(?=^\^|\Z)", payload):
        gsm, body = match.groups()
        titles = re.findall(r"(?m)^!Sample_title = ([^\r\n]+)", body)
        labels = re.findall(r"(?m)^!Sample_characteristics_ch1 = outcome: ([^\r\n]+)", body)
        if (
            len(titles) != 1
            or not re.fullmatch(r"P[0-9]+", titles[0])
            or len(labels) != 1
            or labels[0] not in OUTCOMES
        ):
            raise ValueError(f"Ambiguous GEO donor/outcome for {gsm}")
        patient = titles[0]
        links = {}
        for u in re.findall(r"(?m)^!Sample_supplementary_file = ([^\r\n]+)", body):
            fields = urlsplit(u)
            filename = Path(fields.path).name
            for key, suffix in [
                ("counts", "cell_feature_matrix.h5"),
                ("cells", "cells.parquet.gz"),
            ]:
                if filename == f"{gsm}_{patient}_{suffix}":
                    if fields.hostname != "ftp.ncbi.nlm.nih.gov" or key in links:
                        raise ValueError(f"Ambiguous archive asset {filename}")
                    links[key] = "https://" + fields.netloc + fields.path
        if set(links) != {"counts", "cells"}:
            raise ValueError(f"Missing explicit matrix or cell file URL for {gsm}")
        sections.append(
            {
                "gsm": gsm,
                "patient_id": patient,
                "group": OUTCOMES[labels[0]],
                "original_outcome": labels[0],
                "urls": links,
            }
        )
    if (
        not sections
        or len({s["gsm"] for s in sections}) != len(sections)
        or len({s["patient_id"] for s in sections}) != len(sections)
    ):
        raise ValueError("Missing or duplicated GEO record / patient identity")
    return sorted(sections, key=lambda s: s["gsm"])


def check_geo(snapshot, payload):
    observed = parse_geo(payload.decode("utf-8-sig"))
    expected = [
        {k: r[k] for k in ("gsm", "patient_id", "group", "original_outcome", "urls")}
        for r in snapshot["patients"]
    ]
    if observed != expected:
        raise ValueError("Current GEO design or asset mapping differs from the frozen snapshot")


def select_subset(rows):
    counts: Counter[str] = Counter()
    selected = []
    for r in sorted(rows, key=lambda row: row["gsm"]):
        if counts[r["group"]] < 3:
            counts[r["group"]] += 1
            selected.append(r)
    if counts != {"pCR": 3, "RD": 3}:
        raise ValueError("Cannot form prespecified three-per-outcome subset")
    return selected


def section(matrix_file, cells_file, row, panel):
    import pyarrow.parquet as pq

    with h5py.File(matrix_file) as handle:
        m = handle["matrix"]
        barcodes = pd.Index(m["barcodes"].asstr()[:], name="cell_id")
        names = m["features/name"].asstr()[:]
        types = m["features/feature_type"].asstr()[:]
        biological = np.flatnonzero(types == "Gene Expression")
        symbols = pd.Index(names[biological])
        # The explicit panel is only valid when each selected biological symbol
        # occurs once. Identical display names among control probes are irrelevant.
        hits = []
        for symbol in panel:
            found = np.flatnonzero(symbols == symbol)
            if len(found) != 1:
                raise ValueError(f"Ambiguous or absent gene symbol {symbol}")
            hits.append(biological[found[0]])
        shape = tuple(map(int, m["shape"][:]))
        full = sparse.csr_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape[::-1])
        library_size = np.asarray(full @ (types == "Gene Expression").astype(np.int64)).ravel()
        matrix = full[:, hits].astype(np.int32)
    with gzip.open(cells_file, "rb") as fh:
        table = pq.read_table(
            fh, columns=["cell_id", "x_centroid", "y_centroid", "transcript_counts"]
        ).to_pandas()
    cells = table.set_index("cell_id")
    if (
        not cells.index.is_unique
        or not barcodes.is_unique
        or len(barcodes) != len(cells)
        or not barcodes.isin(cells.index).all()
    ):
        raise ValueError("Xenium cell IDs and count matrix barcodes disagree")
    xy = cells.loc[barcodes, ["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError("Nonfinite cell centroid")
    if not np.array_equal(library_size, cells.loc[barcodes, "transcript_counts"].to_numpy()):
        raise ValueError(
            "All biological target counts must reproduce the vendor cell transcript total"
        )
    sid = f"{row['gsm']}_{row['patient_id']}"
    obs = pd.DataFrame(
        {"patient_id": row["patient_id"], "section_id": sid, "condition": row["group"]},
        index=barcodes,
    )
    obj = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(panel, name="symbol")))
    obj.obs["library_size"] = library_size
    obj.obsm["spatial"] = xy
    obj.uns["measurement"] = "X is raw integer transcript count on frozen panel"
    obj.uns["unit"] = "um"
    return obj


def compare(obj, path):
    from anndata.io import read_elem

    with h5py.File(path) as h:
        var = h["var"]
        ids = pd.Index(read_elem(var[var.attrs["_index"]])).get_indexer(obj.var_names)
        if (ids < 0).any():
            raise ValueError("Reference omits a panel feature")
        stored = read_elem(h["layers/counts"])
        if not sparse.issparse(stored):
            raise ValueError("Expected sparse reference count layer")
        counts = sparse.csr_matrix(stored)[:, ids]
        ob = h["obs"]
        tests = {
            "counts_equal": (counts != obj.X).nnz == 0,
            "coords_um_equal": np.array_equal(h["obsm/spatial"][:], obj.obsm["spatial"]),
            "barcodes_equal": np.array_equal(
                obj.obs_names, np.asarray(read_elem(ob[ob.attrs["_index"]]))
            ),
            "donor_group_equal": all(
                np.array_equal(obj.obs[k], np.asarray(read_elem(ob[k])))
                for k in ["patient_id", "condition"]
            ),
        }
    if not all(tests.values()):
        raise ValueError(f"Scientific arrays differ in {path.name}: {tests}")
    return tests


def prepare(
    raw, out, download=False, verify_geo=False, reference=None, evidence=None, subset_only=False
):
    out.mkdir(parents=True, exist_ok=True)
    snapshot = strict_read(DATA / "geo_samples.json")
    if verify_geo:
        with urlopen(
            Request(METADATA_URL, headers={"User-Agent": "spatial-audit/0.1"}), timeout=90
        ) as response:
            check_geo(snapshot, response.read())
    panel = (DATA / "feature_panel.txt").read_text(encoding="utf-8").splitlines()
    if not panel or len(set(panel)) != len(panel):
        raise ValueError("Empty or duplicated gene panel")
    rows = snapshot["patients"]
    subset = select_subset(rows)
    chosen = {r["gsm"] for r in subset}
    manifest_files = {}
    receipts = []
    with TemporaryDirectory(prefix="spatial-audit-xenium-", dir=out) as work:
        for row in subset if subset_only else rows:
            sid = f"{row['gsm']}_{row['patient_id']}"
            paths = {}
            originals = []
            for key, pin in row["assets"].items():
                name = pin["name"]
                url = row["urls"][key]
                if name != Path(urlsplit(url).path).name:
                    raise ValueError("Pinned file and GEO asset name disagree")
                paths[key] = verified_file(raw / sid / name, url, pin, download)
                originals.append({**pin, "url": url})
            obj = section(paths["counts"], paths["cells"], row, panel)
            receipt = {
                "gsm": row["gsm"],
                "patient_id": row["patient_id"],
                "group": row["group"],
                "cells": obj.n_obs,
                "raw_files": originals,
                "subset_member": row["gsm"] in chosen,
            }
            if reference:
                receipt["reference_comparison"] = compare(obj, reference / (sid + ".h5ad"))
            file = Path(work) / (sid + ".h5ad")
            obj.write_h5ad(file, compression="gzip")
            digest = sha256(file)
            dest = f"{sid}-{digest[:16]}.h5ad"
            file.replace(out / dest)
            manifest_files[row["gsm"]] = {"path": dest, "sha256": digest}
            receipts.append({**receipt, "prepared_file": dest, "prepared_sha256": digest})
    base = {
        "platform": "Xenium",
        "kind": "public",
        "units": "um",
        "source": snapshot["source"],
        "features": panel,
        "columns": {"group": "condition"},
        "measurement": "counts",
        "library_size_column": "library_size",
        "matrix_mode": "memory",
    }
    for is_subset in [True] if subset_only else [False, True]:
        selected = subset if is_subset else rows
        manifest = {
            **base,
            "id": ACCESSION + ("-subset" if is_subset else ""),
            "selection": (
                "First three distinct GSM records per archived response group, lexical GSM order"
                if is_subset
                else "All distinct archived Xenium donors"
            )
            + "; frozen common biological expression panel.",
            "files": [manifest_files[r["gsm"]] for r in selected],
        }
        write_json(
            out / ("manifest-subset.json" if is_subset else "manifest.json"),
            {"cohorts": [manifest]},
        )
    record = {
        "kind": "public",
        "record_type": "archive_reconstruction",
        "accession": ACCESSION,
        "n_patients": len(receipts),
        "n_observations": sum(r["cells"] for r in receipts),
        "features": panel,
        "metadata_source": METADATA_URL,
        "metadata_snapshot_sha256": snapshot["source_sha256"],
        "mapping": "Each GEO Sample_title supplies a donor ID. Sample_characteristics_ch1 outcome is matched exactly to its deposited label. Each vendor H5 barcode is joined to cells.parquet cell_id.",
        "coordinates": "Vendor x_centroid and y_centroid fields in micrometres.",
        "library_size": "Sum of all biological Gene Expression counts in each vendor HDF5 barcode column; compared for exact equality with the raw cells.parquet transcript_counts field.",
        "patients": receipts,
    }
    write_json(out / "dataset_manifest.json", record)
    if evidence:
        write_json(evidence, record)
    print(f"Prepared {len(receipts)} Xenium donors, {record['n_observations']} cells")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--download",
        action="store_true",
        help="Fetch missing matrix/Parquet files via official GEO sample URLs",
    )
    p.add_argument(
        "--verify-geo",
        action="store_true",
        help="Compare the current GEO sample-level SOFT design with the snapshot",
    )
    p.add_argument(
        "--subset-only",
        action="store_true",
        help="Prepare the frozen three-per-group comparator subset",
    )
    p.add_argument(
        "--reference",
        type=Path,
        help="Compare all selected arrays with an existing per-GSM H5AD directory",
    )
    p.add_argument(
        "--evidence", type=Path, help="Save a portable duplicate of the reconstruction evidence"
    )
    args = p.parse_args()
    prepare(
        args.raw_dir,
        args.out,
        args.download,
        args.verify_geo,
        args.reference,
        args.evidence,
        args.subset_only,
    )
