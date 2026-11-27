"""Prepare the anti-GBM versus reference renal Xenium contrast from GSE294965."""

import argparse
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from anndata.io import read_elem, sparse_dataset
from fetch_and_verify import fetch

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE294965"


def section_value(codes, names, ids):
    distinct = np.unique(codes[ids])
    if len(distinct) != 1 or distinct[0] < 0:
        raise ValueError("Biopsy assignment is incomplete or inconsistent")
    return str(names[distinct[0]])


def select_rows(node, ids, columns=None):
    blocks = []
    for positions in np.split(ids, np.flatnonzero(np.diff(ids // 8192)) + 1):
        lo, hi = int(positions[0]), int(positions[-1]) + 1
        sub = node[lo:hi]
        if columns is not None:
            sub = sub[:, columns].toarray()
        blocks.append(sub[positions - lo])
    return np.concatenate(blocks)


def prepare(cache, manifest, evidence):
    row = next(
        c for c in strict_read(ROOT / "docs/data/catalog.json")["cohorts"] if c["id"] == ACCESSION
    )
    asset = next(a for a in row["assets"] if a["name"].endswith(".h5ad"))
    source = cache / asset["name"]
    receipts = fetch({"files": [dict(asset, accession=ACCESSION)]}, cache, local=source.is_file())
    manifest.parent.mkdir(parents=True, exist_ok=True)
    dest = manifest.parent / "renal-sections"
    dest.mkdir(exist_ok=True)
    sections = []
    groups: dict[str, str] = {}
    samples = []
    with h5py.File(source) as h:
        obs = h["obs"]

        def codes(k):
            node = obs[k]
            return node["codes"][:], node["categories"].asstr()[:]

        label, group_names = codes("Disease")
        patient, patient_names = codes("Patient_Sample_ID")
        biopsy, biopsy_names = codes("Biopsy_ID")
        slide, slide_names = codes("Slide_ID")
        batch, batch_names = codes("batch")
        celltype, types = codes("celltype_l1")
        include = np.isin(label, np.flatnonzero(np.isin(group_names, ["Cntrl", "GBM"])))
        features = np.asarray(read_elem(h["var"][h["var"].attrs["_index"]]), str)
        biological = np.flatnonzero(
            [not s.lower().startswith(("negcontrol", "blank", "unassigned")) for s in features]
        )
        chosen = biological[:32]
        matrix: Any = sparse_dataset(h["layers/counts"])
        library = np.empty(len(patient), dtype=float)
        for start in range(0, len(patient), 8192):
            block = matrix[start:start+8192, :][:, biological]
            library[start:start+block.shape[0]] = np.asarray(block.sum(axis=1)).ravel()
        for code in np.unique(biopsy[include]):
            ids = np.flatnonzero(include & (biopsy == code))

            pid = section_value(patient, patient_names, ids)
            grp = section_value(label, group_names, ids)
            if groups.setdefault(pid, grp) != grp:
                raise ValueError("Patient crosses disease groups")
            sid = section_value(slide, slide_names, ids)
            b = section_value(batch, batch_names, ids)
            xy = select_rows(h["obsm/spatial"], ids)
            if not np.allclose(
                xy,
                np.column_stack([select_rows(obs[k], ids) for k in ("x_centroid", "y_centroid")]),
            ):
                raise ValueError("Spatial embedding not equal to archived Xenium centroids")
            X = select_rows(matrix, ids, chosen)
            if X.min() < 0 or not np.allclose(X, np.rint(X)):
                raise ValueError("Raw count layer required")
            C = (celltype[ids, None] == np.arange(len(types))[None, :]).astype(np.uint8)
            local = dest / f"biopsy-{int(code)}.npz"
            np.savez_compressed(local, X=X, coords=xy, C=C, features=features[chosen], library_size=library[ids])
            sections.append(
                {
                    "path": local.relative_to(manifest.parent).as_posix(),
                    "group": grp,
                    "patient_id": pid,
                    "section_id": str(biopsy_names[code]),
                    "slide_id": sid,
                    "batch": b,
                    "sha256": sha256(local),
                }
            )
            samples.append({"patient": pid, "group": grp, "slide": sid, "observations": len(ids)})
    if sorted(list(groups.values()).count(g) for g in set(groups.values())) != [6, 6]:
        raise ValueError("Reconcile anti-GBM and control donor counts with official study design")
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": ACCESSION,
                    "kind": "public",
                    "source": row["source"],
                    "platform": "Xenium",
                    "units": "um",
                    "measurement": "counts",
                    "selection": "All anti-GBM and control patients; first 32 archived biological panel targets, Xenium micrometre centroids, one-hot archived cell types.",
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
            "accession": ACCESSION,
            "files": receipts,
            "source": row["source"],
            "paper": row["paper"],
            "features": features[chosen].tolist(),
            "cells_per_biopsy": samples,
            "covariates": types.tolist(),
            "library_size": "Total raw counts across archived biological targets before the 32-gene subset.",
            "annotation_note": "Assignments are derived from the archived molecular profiles; projection can include circularity.",
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
