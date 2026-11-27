"""Use archive-split patient and specimen labels to audit lesional psoriasis versus AD."""

import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from fetch_and_verify import fetch
from prepare_visium import GEOMETRY, PANEL

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "GSE206391"


def clinical_map(text):
    records = {}
    for part in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        title = re.search(r"!Sample_title = (.+)", part)
        match = re.search(r"\[([^\]]+)\]", title[1]) if title else None
        meta = dict(re.findall(r"(?m)^!Sample_characteristics_ch1 = ([^:\n]+): (.+)$", part))
        if not match or match[1] in records:
            raise ValueError("GEO sample label missing or duplicated")
        records[match[1]] = {
            "patient": meta["patient"],
            "group": meta["disease"],
            "biopsy": meta["tissue type"],
            "gsm": part.splitlines()[0].strip(),
        }
    return records


def prepare(cache, manifest, evidence):
    catalog = strict_read(ROOT / "docs/data/catalog.json")["cohorts"]
    row = next(c for c in catalog if c["id"] == ACCESSION)
    file = next(a for a in row["assets"] if "Preprocessed_data" in a["name"])
    target = cache / file["name"]
    files = fetch({"files": [dict(file, accession=ACCESSION)]}, cache, local=target.is_file())
    meta = cache / "samples.soft"
    files += fetch(
        {
            "files": [
                {
                    "name": meta.name,
                    "url": row["source"] + "&targ=gsm&form=text&view=full",
                    "accession": ACCESSION,
                }
            ]
        },
        cache,
        local=meta.is_file(),
    )
    clinical = clinical_map(meta.read_text(encoding="utf8"))
    dest = manifest.parent / "skin-sections"
    dest.mkdir(parents=True, exist_ok=True)
    sections = []
    mapping = []
    with h5py.File(target) as h:
        obs = h["obs"]

        def labels(k):
            codes = obs[k][:]
            if np.any(codes < 0):
                raise ValueError("Missing observed category " + k)
            return obs["__categories"][k][:].astype(str)[codes]

        patient, sample, biopsy, slide, batch = map(
            labels, ["patient", "sample", "biopsy_type", "object_slide", "project"]
        )
        rc = np.column_stack([obs["array_row"][:], obs["array_col"][:]]).astype(float)
        coord = np.column_stack((rc[:, 1] * 50, rc[:, 0] * 50 * np.sqrt(3)))
        names = pd.Index(h["var/gene_name"].asstr()[:])
        panel = [g for g in PANEL if sum(names == g) == 1]
        columns = [np.flatnonzero(names == g)[0] for g in panel]
        matrix = h["layers/counts"]
        X = np.empty((len(patient), len(panel)), np.float32)
        depth = np.empty(len(patient))
        for j in range(0, len(patient), 256):
            data = matrix[j : j + 256, :]
            if not np.isfinite(data).all() or np.any(data < 0) or np.any(data != np.rint(data)):
                raise ValueError("Deposited counts must be finite nonnegative integers")
            X[j : j + len(data)] = data[:, columns]
            depth[j : j + len(data)] = data.sum(axis=1, dtype=np.float64)
        for label in sorted(set(sample)):
            archive_key = label if label in clinical else label.split("_", 1)[-1]
            if archive_key not in clinical:
                raise ValueError("Library missing clinical metadata: " + label)
            c = clinical[archive_key]
            if (
                c["group"] not in ("psoriasis", "atopic dermatitis")
                or c["biopsy"] != "lesional"
                or not c["patient"].isdigit()
            ):
                continue
            idx = sample == label
            if set(biopsy[idx]) != {"LESIONAL"} or set(patient[idx]) != {c["patient"]}:
                raise ValueError("Independent GEO and per-spot patient maps disagree")
            positions = rc[idx]
            if (
                np.any(positions % 1)
                or np.any(positions < 0)
                or np.any(positions[:, 0] >= 78)
                or np.any(positions[:, 1] >= 128)
                or np.any(positions.sum(axis=1) % 2)
                or len(np.unique(positions, axis=0)) != idx.sum()
            ):
                raise ValueError("Unexpected Visium array geometry")
            sids = np.unique(slide[idx])
            batches = np.unique(batch[idx])
            if len(sids) != 1 or len(batches) != 1:
                raise ValueError("Ambiguous slide or project within section")
            path = dest / (c["gsm"] + ".npz")
            np.savez_compressed(
                path, features=np.array(panel), coords=coord[idx], X=X[idx], library_size=depth[idx]
            )
            entry = {
                "patient_id": c["patient"],
                "group": c["group"],
                "section_id": c["gsm"],
                "slide_id": str(sids[0]),
                "batch": str(batches[0]),
            }
            sections.append(
                {
                    **entry,
                    "path": path.relative_to(manifest.parent).as_posix(),
                    "sha256": sha256(path),
                }
            )
            mapping.append({**entry, "library": label, "spots": int(idx.sum())})
    selection = "All single-donor lesional capture areas with concordant GEO and per-spot patient assignments for psoriasis or atopic dermatitis. Repeated sections retain their donor ID. NPZ library_size sums the complete archived counts layer before panel selection; square-grid audit uses the physical Visium hex-array addresses."
    spec = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "kind": "public",
        "source": row["source"],
        "platform": "Visium",
        "units": "um",
        "measurement": "counts",
        "selection": selection,
        "files": sections,
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        evidence,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "accession": ACCESSION,
            "source": row["source"],
            "paper": row["paper"],
            "files": files,
            "selection": selection,
            "physical_calibration": GEOMETRY,
            "mapping": mapping,
            "panel": panel,
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("cache", "manifest", "provenance"):
        p.add_argument("--" + key, required=True, type=Path)
    args = p.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
