"""Within-response donor randomization and nested-panel sensitivity on verified public inputs."""

import argparse
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import ndtr

from spatial_audit.pseudorep import bh, patient_test
from spatial_audit.schema import load_manifest, sha256, strict_read, write_json
from spatial_audit.spectral_audit import section_statistics

PUBLIC = Path(__file__).resolve().parents[1] / "results/public"


def panel_sensitivity():
    paths = [PUBLIC / f"{name}.json" for name in ("GSE220442", "GSE220442-full-panel")]
    small, large = (strict_read(p) for p in paths)
    if (
        small["provenance"]["inputs"] != large["provenance"]["inputs"]
        or any(
            small[k] != large[k]
            for k in ("accession", "n_patients", "n_sections", "n_observations")
        )
        or small["dimensions"]["pseudorep"]["features"]
        != large["dimensions"]["pseudorep"]["features"][: small["n_features"]]
    ):
        raise ValueError("Nested panels must use identical inputs and observations")
    panels = []
    for report, path in zip((small, large), paths, strict=True):
        dims = report["dimensions"]
        pseudo, nugget = dims["pseudorep"], dims["nugget"]
        panels.append(
            {
                "id": report["id"],
                "report_sha256": sha256(path),
                "features": report["n_features"],
                "status": {k: v["status"] for k, v in dims.items()},
                "pooled_rate": pseudo["spot_rate"],
                "patient_rate": pseudo["patient_rate"],
                "bands": [
                    {k: b[k] for k in ("wavelength_um", "share", "estimable_features")}
                    for b in dims["composition"]["bands"]
                ],
                "nugget_flagged": nugget["flagged"],
                "nugget_evaluated": nugget["evaluated"],
            }
        )
    a, b = (r["dimensions"] for r in (small, large))
    changed = {}
    for method, key in (("pooled", "spot_q"), ("patient", "q")):
        tests = [
            d["pseudorep"] if method == "pooled" else d["pseudorep"]["patient"] for d in (a, b)
        ]
        called = [np.array(t[key][: small["n_features"]]) < 0.05 for t in tests]
        changed[method] = [
            f
            for f, flag in zip(a["pseudorep"]["features"], called[0] != called[1], strict=True)
            if flag
        ]
    return {
        "kind": "public",
        "record_type": "panel_sensitivity",
        "accession": small["accession"],
        "source": small["source"],
        "n_patients": small["n_patients"],
        "n_observations": small["n_observations"],
        "alpha": 0.05,
        "panels": panels,
        "shared_features_with_changed_BH_call": changed,
        "matching_dimensions": [k for k in a if a[k]["status"] == b[k]["status"]],
    }


def rank_summaries(sections):
    """Precompute donor sums of pooled midranks and the asymptotic MW tie correction."""
    sizes = np.array([s.X.shape[0] for s in sections])
    features = sections[0].X.shape[1]
    ranks = np.empty((len(sizes), features))
    variances = np.empty(features)
    starts = np.r_[0, np.cumsum(sizes)[:-1]]
    total = int(sizes.sum())
    for j in range(features):
        cols = [s.X[:, j : j + 1] for s in sections]
        data = np.concatenate(
            [np.asarray(c.toarray() if sparse.issparse(c) else c).ravel() for c in cols]
        )
        if np.any(data < 0) or np.any(data != np.floor(data)):
            raise ValueError("Integer counts required for histogram midranks")
        counts = data.astype(np.int64)
        hist = np.bincount(counts).astype(float)
        mid = hist.cumsum() - (hist - 1) / 2
        ranks[:, j] = np.add.reduceat(mid[counts], starts)
        variances[j] = ((total + 1) - np.sum(hist**3 - hist) / (total * (total - 1))) / 12
    return sizes, ranks, variances


def pooled_p(mask, sizes, ranks, variances):
    n, m = int(sizes[mask].sum()), int(sizes[~mask].sum())
    u = ranks[mask].sum(axis=0) - n * (n + 1) / 2
    sd = np.sqrt(n * m * np.maximum(0, variances))
    z = np.divide(abs(u - n * m / 2) - 0.5, sd, out=np.zeros_like(u), where=sd > 0)
    return np.minimum(1.0, 2 * ndtr(-z))


def experiment(manifest, splits, draws, seed):
    (cohort,) = load_manifest(manifest)
    receipt_path = PUBLIC / "GSE348499_reconstruction.json"
    receipt = strict_read(receipt_path)
    local_receipt = strict_read(manifest.parent / "dataset_manifest.json")

    def source_facts(r):
        return {
            d["patient_id"]: (d["gsm"], d["group"], d["cells"], d["raw_files"])
            for d in r["patients"]
        }

    if (
        cohort.id != receipt["accession"]
        or cohort.features != receipt["features"]
        or source_facts(receipt) != source_facts(local_receipt)
        or sorted(i["sha256"] for i in cohort.provenance["inputs"])
        != sorted(d["prepared_sha256"] for d in local_receipt["patients"])
        or sorted((s.patient_id, s.group, s.X.shape[0]) for s in cohort.sections)
        != sorted((p["patient_id"], p["group"], p["cells"]) for p in receipt["patients"])
    ):
        raise ValueError("Cohort, panel, donor membership or reconstruction receipt mismatch")
    rng = np.random.default_rng(seed)
    arms = []
    for group in sorted({s.group for s in cohort.sections}):
        sections = sorted(
            (s for s in cohort.sections if s.group == group), key=lambda s: s.patient_id
        )
        sizes, ranks, variances = rank_summaries(sections)
        endpoints = {
            "patient_mean": np.array([np.asarray(s.X.mean(axis=0)).ravel() for s in sections])
        }
        for mode in ("raw", "vst"):
            endpoints[mode + "_corrected"] = np.array(
                [
                    section_statistics(s, (250.0, 600.0, 1400.0), nugget_mode=mode)[0][-3]
                    for s in sections
                ]
            )
        # Register complete-support features once, before drawing labels.
        tested = {key: np.isfinite(values).all(axis=0) for key, values in endpoints.items()}
        if not all(m.any() for m in tested.values()):
            raise ValueError("No complete-support features for an endpoint")
        families = {
            "pooled": cohort.features,
            **{key: np.array(cohort.features)[m].tolist() for key, m in tested.items()},
        }
        rows: list[dict] = []
        for _ in range(splits):
            mask = np.zeros(len(sizes), dtype=bool)
            mask[rng.choice(len(mask), len(mask) // 2, replace=False)] = True
            perm_seed = int(rng.integers(2**31))
            ps = {"pooled": pooled_p(mask, sizes, ranks, variances)}
            for key, values in endpoints.items():
                ps[key] = np.array(
                    patient_test(
                        values[:, tested[key]], mask, max_exact=1, draws=draws, seed=perm_seed
                    )["p"]
                )
            rows.append(
                {
                    "arm_a": [
                        s.patient_id for s, select in zip(sections, mask, strict=True) if select
                    ],
                    "permutation_seed": perm_seed,
                    "bh_calls": {k: int((bh(p) < 0.05).sum()) for k, p in ps.items()},
                    "nominal_calls": {k: int((p < 0.05).sum()) for k, p in ps.items()},
                }
            )
        summary = {}
        for k, family in families.items():
            counts = np.array([r["bh_calls"][k] for r in rows])
            summary[k] = {
                "tests": splits * len(family),
                "bh_calls": int(counts.sum()),
                "bh_fraction": float(counts.mean() / len(family)),
                "splits_with_calls": int((counts > 0).sum()),
                "nominal_fraction": sum(r["nominal_calls"][k] for r in rows)
                / (splits * len(family)),
            }
        arms.append(
            {
                "group": group,
                "donors": [s.patient_id for s in sections],
                "n_cells": int(sizes.sum()),
                "group_sizes": [len(sizes) // 2, len(sizes) - len(sizes) // 2],
                "families": families,
                "summary": summary,
                "rows": rows,
            }
        )
    return {
        "kind": "public",
        "record_type": "empirical_null",
        "accession": cohort.id,
        "source": cohort.source,
        "reconstruction_sha256": sha256(receipt_path),
        "prepared_inputs": cohort.provenance["inputs"],
        "seed": seed,
        "splits_per_arm": splits,
        "permutation_draws": draws,
        "alpha": 0.05,
        "arms": arms,
        "inference": "Random half-splits within each original outcome arm. All cells retained. "
        "Continuity- and tie-corrected pooled Mann-Whitney, or donor-level moderated permutation "
        "on mean counts and signed nugget-corrected energies in raw and Anscombe modes. "
        "BH within each split and endpoint. Repeated splits share donors and features.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, required=True, help="Full prepare_gse348499.py manifest"
    )
    parser.add_argument("--splits", type=int, default=50)
    parser.add_argument("--draws", type=int, default=2047)
    parser.add_argument("--seed", type=int, default=9721)
    args = parser.parse_args()
    if min(args.splits, args.draws) < 1 or args.seed < 0:
        parser.error("invalid randomization budget")
    write_json(PUBLIC / "panel_sensitivity.json", panel_sensitivity())
    write_json(
        PUBLIC / "empirical_null.json",
        experiment(args.manifest, args.splits, args.draws, args.seed),
    )
