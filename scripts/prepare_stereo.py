"""Rebuild the treated dMMR Stereo-seq contrast in the public STT0000036 archive.

Use the deposited bin100 expression values for the primary analysis. An optional
second manifest contains integer counts only if every archived nonzero and
every complete per-bin library total pass the inversion checks in every section.
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
ACCESSION = "STT0000036"
PITCH_UM = 50.0
NORMALIZATION_TARGET = 10_000
LOG_ATOL = 5e-6
GENES = [
    "ACTB",
    "GAPDH",
    "B2M",
    "RPLP0",
    "RPL13A",
    "RPS18",
    "RPS3",
    "RPL10",
    "EEF1A1",
    "PPIA",
    "HPRT1",
    "CALM1",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD74",
    "HLA-DRA",
    "PTPRC",
    "LST1",
    "TYROBP",
    "MS4A1",
    "CD79A",
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM",
    "PECAM1",
    "VWF",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
]
COVARIATES = ["Epi", "Endo", "CD8_Teff", "Fibro_CXCL14", "Mac_SPP1", "B_naive"]
GROUPS = {"dR": "anti-PD1 responder", "dNR": "anti-PD1 nonresponder"}
FIELDS = ["id", "Type", "Treatment", "group1", "group2"]


def select_sections(assets: dict, cache: Path) -> tuple[list, dict]:
    """Select all treated dMMR samples by archived labels, checking every bin."""
    records = []
    selected = {}
    for name in sorted(assets):
        if not name.startswith("STmetadata_") or not name.endswith(".csv"):
            continue
        df = pd.read_csv(cache / name, usecols=FIELDS)
        unique = df.drop_duplicates()
        if len(unique) != 1 or unique.isna().any().any():
            raise ValueError(f"Missing or contradictory sample labels in {name}")
        row = unique.iloc[0].to_dict()
        section = str(row["id"])
        if name != f"STmetadata_{section}.csv":
            raise ValueError(f"Sample name and id differ in {name}")
        include = row["Type"] == "dMMR" and row["Treatment"] == "pd1"
        records.append(dict(**row, bins=len(df), selected=include))
        if include:
            group = GROUPS[row["group2"]]
            gem = f"STexpression_{section}.gem.gz"
            if gem not in assets:
                raise ValueError(f"No archived matrix for {section}")
            selected[section] = (group, name, gem)
    if {row[0] for row in selected.values()} != set(GROUPS.values()):
        raise ValueError("Both response categories must be represented")
    return records, selected


def get_panel(metadata: pd.DataFrame, gem: Path, chunk_size: int = 300_000):
    """Join all matrix entries to unique bins and test any inverse-log transformation."""
    coordinates = metadata[["x", "y"]].to_numpy(dtype=float)
    if not np.isfinite(coordinates).all() or not np.equal(coordinates, np.rint(coordinates)).all():
        raise ValueError("Expected finite integer bin indices")
    xy = coordinates.astype("int64")
    if np.any(xy < 0) or len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError("Expected unique nonnegative bin coordinates")
    # Axis keys avoid a dense bounding box when tissue occupies disconnected fields.
    stride = int(xy[:, 1].max()) + 1
    keys = pd.Index(xy[:, 0] * stride + xy[:, 1])
    depth = metadata["nCount_Spatial"].to_numpy(float)
    if (
        not np.isfinite(depth).all()
        or np.any(depth < 1)
        or not np.equal(depth, np.rint(depth)).all()
    ):
        raise ValueError("Expected positive integer library sizes for all analysed bins")
    panel = np.zeros((len(xy), len(GENES)), dtype=np.float32)
    counts = np.zeros(panel.shape, dtype=np.int64)
    depth_sum = np.zeros(len(keys), dtype=np.int64)
    nnz_per_bin = np.zeros(len(keys), dtype=np.int64)
    max_log_error = 0.0
    max_round_error = 0.0
    invalid_entries = 0
    entries = 0
    target_index = pd.Index(GENES)
    for chunk in pd.read_csv(gem, sep="\t", chunksize=chunk_size, keep_default_na=False):
        xx = chunk["x"].to_numpy()
        yy = chunk["y"].to_numpy()
        if not (np.issubdtype(xx.dtype, np.integer) and np.issubdtype(yy.dtype, np.integer)):
            raise ValueError("Non-integer bin coordinate in GEM")
        ix = keys.get_indexer(xx * stride + yy)
        if (ix < 0).any() or (yy >= stride).any() or (yy < 0).any() or (xx < 0).any():
            raise ValueError("GEM coordinate missing from metadata")
        v = chunk["MIDCounts"].to_numpy(float)
        if not np.isfinite(v).all() or (v <= 0).any():
            raise ValueError("Expected positive finite sparse archive entries")
        recovered = np.expm1(v) * depth[ix] / NORMALIZATION_TARGET
        if not np.isfinite(recovered).all() or (recovered > np.iinfo(np.int32).max).any():
            raise ValueError("Unsupported expression scale in sparse archive")
        integer = np.rint(recovered).astype(np.int64)
        log_error = np.abs(np.log1p(integer * NORMALIZATION_TARGET / depth[ix]) - v)
        max_log_error = max(max_log_error, float(log_error.max()))
        max_round_error = max(max_round_error, float(np.abs(recovered - integer).max()))
        invalid_entries += int(((log_error > LOG_ATOL) | (integer <= 0)).sum())
        depth_sum += np.bincount(ix, weights=integer, minlength=len(keys)).astype(np.int64)
        nnz_per_bin += np.bincount(ix, minlength=len(keys))
        entries += len(v)
        which = target_index.get_indexer(chunk["geneID"])
        hit = which >= 0
        rows, cols = ix[hit], which[hit]
        encoded = rows * len(GENES) + cols
        if len(np.unique(encoded)) != len(rows) or np.any(panel[rows, cols] != 0):
            raise ValueError("Repeated panel gene at a single archived bin")
        panel[rows, cols] = v[hit]
        counts[rows, cols] = integer[hit]
    bad_totals = np.count_nonzero(depth_sum != depth)
    bad_detected = np.count_nonzero(nnz_per_bin != metadata["nFeature_Spatial"].to_numpy())
    proof = {
        "sparse_entries_checked": entries,
        "library_total_mismatch_bins": int(bad_totals),
        "nfeature_mismatch_bins": int(bad_detected),
        "invalid_inverse_entries": invalid_entries,
        "max_inverse_log_error": max_log_error,
        "max_count_rounding_residual": max_round_error,
        "validated_integer_reconstruction": bool(
            entries and not (invalid_entries or bad_totals or bad_detected)
        ),
    }
    return xy.astype(float) * PITCH_UM, panel, counts, proof


def save_section(path, coords, matrix, c, domain, library_size=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz", delete=False) as stream:
        temp = Path(stream.name)
    try:
        np.savez_compressed(
            temp,
            coords=coords,
            X=matrix,
            C=c,
            features=np.asarray(GENES),
            domain=domain,
            **({"library_size": library_size} if library_size is not None else {}),
        )
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def prepare(
    cache: Path,
    manifest: Path,
    evidence: Path,
    counts_manifest: Path | None,
    download: bool = False,
) -> None:
    outputs = [manifest, evidence] + ([counts_manifest] if counts_manifest is not None else [])
    if len(set(outputs)) != len(outputs):
        raise ValueError("Use distinct output paths for manifests and evidence")
    catalog = strict_read(ROOT / "docs/data/catalog.json")
    item = next(row for row in catalog["cohorts"] if row["id"] == ACCESSION)
    assets = {asset["name"]: asset for asset in item["assets"]}
    if download:
        from fetch_and_verify import fetch

        meta = [{**a, "accession": ACCESSION} for a in assets.values() if a["role"] == "metadata"]
        fetch({"files": meta}, cache)
    selection, selected = select_sections(assets, cache)
    if download:
        gems = [{**assets[row[2]], "accession": ACCESSION} for row in selected.values()]
        fetch({"files": gems}, cache)
    patients = [section.split("_")[0] for section in selected]
    if len(set(patients)) != len(patients):
        raise ValueError("Selected contrast needs one section per distinct archived patient")
    files = {f"STmetadata_{r['id']}.csv" for r in selection}
    files.update(names[2] for names in selected.values())
    receipts = [
        {
            **assets[n],
            "accession": ACCESSION,
            "bytes": (cache / n).stat().st_size,
            "sha256": sha256(cache / n),
        }
        for n in sorted(files)
    ]
    sections = []
    count_sections = []
    report = []
    for section, (group, name, gem) in selected.items():
        df = pd.read_csv(cache / name, low_memory=False)
        coords, panel, counts, proof = get_panel(df, cache / gem)
        c = df[COVARIATES].to_numpy(float)
        if not np.isfinite(c).all() or (c < 0).any() or (c > 1 + 1e-6).any():
            raise ValueError(f"Invalid deposited composition estimate in {section}")
        domain = np.asarray(df["level1"].fillna("unassigned").astype(str), dtype=str)
        out = cache / "derived" / (section + ".npz")
        save_section(out, coords, panel, c, domain)
        shared = {"patient_id": section.split("_")[0], "section_id": section, "group": group}
        sections.append(
            {
                **shared,
                "path": os.path.relpath(out, manifest.parent).replace(os.sep, "/"),
                "sha256": sha256(out),
            }
        )
        obs = {
            **shared,
            **proof,
            "bins": len(df),
            "source_expression": gem,
            "source_metadata": name,
            "section_bytes": out.stat().st_size,
            "section_sha256": sections[-1]["sha256"],
            "panel_nonzero_bins": np.count_nonzero(panel, axis=0),
            "known_domain_bins": int(df["level1"].notna().sum()),
        }
        if counts_manifest is not None and proof["validated_integer_reconstruction"]:
            dest = cache / "derived-counts" / (section + ".npz")
            save_section(dest, coords, counts, c, domain, df["nCount_Spatial"].to_numpy(np.int64))
            digest = sha256(dest)
            count_sections.append(
                {
                    **shared,
                    "path": os.path.relpath(dest, counts_manifest.parent).replace(os.sep, "/"),
                    "sha256": digest,
                }
            )
            obs["count_section_sha256"] = digest
        report.append(obs)
        print(section, len(df), proof, flush=True)
    policy = "All treated dMMR Stereo-seq sections from the study. Predefined housekeeping, immune, stromal, and epithelial marker panel. No treatment outcomes or study-derived rankings used for panel selection."
    base = {
        "id": ACCESSION,
        "kind": "public",
        "accession": ACCESSION,
        "units": "um",
        "preprocessing_evidence": evidence.name,
        "source": item["source"],
        "platform": item["platform"],
        "features": GENES,
        "covariates": COVARIATES,
        "selection": policy,
    }
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    **base,
                    "files": sections,
                    "layer": "archived log-normalized GEM MIDCounts",
                    "measurement": "log_normalized",
                }
            ]
        },
    )
    write_json(
        evidence,
        {
            "accession": ACCESSION,
            "kind": "public",
            "source": item["source"],
            "paper": item["paper"],
            "files": receipts,
            "archive_section_screening": selection,
            "sections": report,
            "preprocessing": {
                "selection": policy,
                "target_genes": GENES,
                "covariates": COVARIATES,
                "expression": "Deposited nonzero log-normalized values, with unlisted panel entries set to zero. The deposited MIDCounts field contains real-valued abundances.",
                "optional_count_check": {
                    "inverse": "round(expm1(MIDCounts) * nCount_Spatial / normalization_target)",
                    "normalization_target": NORMALIZATION_TARGET,
                    "maximum_log_roundtrip_tolerance": LOG_ATOL,
                    "per_bin_checks": "Every nonzero must round-trip, all recovered counts must sum to nCount_Spatial, and detected feature count must equal nFeature_Spatial.",
                },
                "composition": "Six deposited RCTD subtype frequencies. These estimates use the spatial transcriptome itself and can explain signal through shared measurements.",
                "frame": "Bin100 (paper Methods), x and y are local bin indices. Instrument stage, slide and batch metadata unavailable.",
                "bin_pitch_um": PITCH_UM,
                "domain": "Deposited level1 bin annotations, missing values marked unassigned.",
            },
        },
    )
    if counts_manifest is not None:
        if len(count_sections) != len(sections):
            raise ValueError(
                "One or more sections failed count reconstruction. The archived scale manifest is available."
            )
        write_json(
            counts_manifest,
            {
                "cohorts": [
                    {
                        **base,
                        "files": count_sections,
                        "layer": "verified reconstructed integer counts",
                        "measurement": "counts",
                    }
                ]
            },
        )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--counts-manifest", type=Path)
    p.add_argument(
        "--fetch",
        action="store_true",
        help="Download catalog metadata and selected expression archives",
    )
    p.add_argument("--provenance", type=Path, required=True)
    args = p.parse_args()
    prepare(
        args.cache.resolve(),
        args.manifest.resolve(),
        args.provenance.resolve(),
        args.counts_manifest.resolve() if args.counts_manifest else None,
        args.fetch,
    )


if __name__ == "__main__":
    main()
