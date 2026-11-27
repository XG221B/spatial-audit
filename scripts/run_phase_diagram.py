"""Generate independent synthetic cohorts on a physical lattice and scan operating characteristics."""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations, product
from math import comb, sqrt
from pathlib import Path

# Each worker runs its own moderately sized dense linear algebra; a nested pool
# of BLAS workers is slower and brings no additional statistical precision.
for envvar in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[envvar] = "1"

import numpy as np

from spatial_audit.schema import strict_read, write_json
from spatial_audit.spectral_audit import build_grid, coefficients, filter_bands

ROOT = Path(__file__).resolve().parent.parent
NGRID = (3, 5, 8, 10, 15, 20, 30)
COMPOSITION = (0.0, 0.1, 0.2, 0.25, 0.4, 0.5, 0.65, 0.75)
DISPERSION = (0.0, 0.025, 0.05, 0.1, 0.2, 0.35)
MEAN_SHIFT = (1.0, 1.1, 1.25, 1.5, 1.75, 2.0)
BATCH_ALIGNMENT = (0.0, 0.25, 0.5, 0.75, 1.0)


@lru_cache
def operators():
    """Cache the same physical finite-volume operator for every simulated section."""
    y, x = np.mgrid[:6, :6]
    xy = np.column_stack((x.ravel(), y.ravel())) * 50.0
    grid = build_grid(xy, 250.0)
    h = filter_bands(grid, np.eye(len(xy)), coefficients([250.0, 600.0, 1400.0], grid.width))[0]
    c = np.cos((x.ravel() + 0.5) * np.pi / 3)
    c = c / np.sqrt(np.mean(c * c))
    signal = np.cos((y.ravel() + 0.5) * np.pi / 3)
    signal /= np.sqrt(np.mean(signal * signal))
    hc = h @ c
    qc = hc / np.linalg.norm(hc)
    residual_h = h - np.outer(qc, qc @ h)
    trace = float(np.sum(residual_h**2) / len(xy))
    return grid, h, residual_h, trace, c, signal


def energy(counts):
    grid, h, residual_h, trace, c, _ = operators()
    field = counts.reshape(-1, counts.shape[-1]).astype(float)
    filtered = field @ h.T
    raw = (filtered**2).mean(axis=1)
    residual_energy = ((field @ residual_h.T) ** 2).mean(axis=1)
    centered = field - field.mean(axis=1, keepdims=True)
    residual = centered - np.outer((centered @ c) / (c @ c), c)
    # On this complete row-major lattice the middle site is the midpoint
    # of the two endpoint indices. Use identical path support at both lags.
    mid = (grid.far_left + grid.far_right) // 2
    g1 = (
        ((residual[:, grid.far_left] - residual[:, mid]) ** 2).mean(axis=1)
        + ((residual[:, mid] - residual[:, grid.far_right]) ** 2).mean(axis=1)
    ) / 4
    g2 = ((residual[:, grid.far_left] - residual[:, grid.far_right]) ** 2).mean(axis=1) / 2
    corrected = residual_energy - (4 * g1 - g2) / 3 * trace
    return raw.reshape(counts.shape[:-1]), corrected.reshape(counts.shape[:-1])


def assignments(n, rng):
    if comb(2 * n, n) <= 512:
        perms = list(combinations(range(2 * n), n))
        masks = np.full((len(perms), 2 * n), -1.0 / n)
        masks[np.arange(len(perms))[:, None], np.array(perms)] = 1.0 / n
        return masks, True
    idx = np.argsort(rng.random((511, 2 * n)), axis=1)[:, :n]
    masks = np.full((len(idx), 2 * n), -1.0 / n)
    masks[np.arange(len(idx))[:, None], idx] = 1.0 / n
    return masks, False


def rejection_rate(values, masks, exact, n):
    reference = np.abs(values @ masks.T)
    observed = abs(values[:, :n].mean(1) - values[:, n:].mean(1))
    exceed = (reference >= observed[:, None] - 1e-11 * np.maximum(1.0, abs(observed[:, None]))).sum(
        1
    )
    p = (exceed + int(not exact)) / (len(masks) + int(not exact))
    return float((p < 0.05).mean())


def wilson(rate, replicates):
    z = 1.959963984540054
    mid = rate + z * z / (2 * replicates)
    radius = z * sqrt((rate * (1 - rate) + z * z / (4 * replicates)) / replicates)
    return [
        max(0.0, (mid - radius) / (1 + z * z / replicates)),
        min(1.0, (mid + radius) / (1 + z * z / replicates)),
    ]


def simulate(job):
    ordinal, params, reps, seed = job
    n, share, disp, shift, alignment = params
    rng = np.random.default_rng(np.random.SeedSequence([seed, ordinal]))
    g = np.repeat([0.0, 1.0], n)
    _, _, _, _, c, signal = operators()
    # Batch label spans completely mixed to completely separated. At partial
    # alignment, realized overlap is measured for each independently drawn cohort.
    prob = 0.5 + alignment * (g - 0.5)
    masks, exact = assignments(n, rng)
    record = {
        "n_per_group": n,
        "composition_setting": share,
        "nb_dispersion": disp,
        "case_mean_multiplier": shift,
        "batch_alignment": alignment,
        "permutation": "enumerated" if exact else "randomized_plus_one",
    }
    for effect_name, amplitude in (("null", 0.0), ("alternative", 7.0)):
        batch = (rng.random((reps, 2 * n)) < prob).astype(float)
        separated = (np.all(batch[:, :n] == 0, axis=1) & np.all(batch[:, n:] == 1, axis=1)) | (
            np.all(batch[:, :n] == 1, axis=1) & np.all(batch[:, n:] == 0, axis=1)
        )
        donor = rng.normal(0, 1.0, (reps, 2 * n))
        # Biological target is a difference in structured amplitude; the null
        # fixes that amplitude with the same population distribution across arms.
        # Batch shifts that same spatial pattern. The within-band composition
        # projection does not identify a patient-level effect under batch aliasing.
        spatial_amplitude = 4.0 + donor + 2.5 * batch + amplitude * g
        composition_amplitude = 4 * np.sqrt(share / (1 - share)) * (0.7 + g)
        mu = (
            25 * (1 + (shift - 1) * g)[None, :, None]
            + spatial_amplitude[:, :, None] * signal
            + composition_amplitude[:, None] * c
        )
        mu = np.maximum(1.0, mu)
        intensity = rng.gamma(1 / disp, mu * disp) if disp else mu
        counts = rng.poisson(intensity)
        raw, corrected = energy(counts)
        for key, arr in (("raw", raw), ("corrected", corrected)):
            rate = rejection_rate(arr, masks, exact, n)
            record[f"{effect_name}_{key}"] = rate
            if key == "corrected":
                record[f"{effect_name}_interval"] = wilson(rate, reps)
        record[f"{effect_name}_separated"] = float(separated.mean())
    return record


def minimum_n(rows):
    keys = ("composition_setting", "nb_dispersion", "case_mean_multiplier", "batch_alignment")
    groups: dict[tuple, list] = {}
    for row in rows:
        groups.setdefault(tuple(row[k] for k in keys), []).append(row)
    lookups = []
    for params, sub in groups.items():
        eligible = [
            r["n_per_group"]
            for r in sub
            if r["batch_alignment"] < 1.0
            and r["alternative_interval"][0] >= 0.8
            and r["null_interval"][1] <= 0.1
        ]
        lookups.append(
            {
                "composition_setting": params[0],
                "nb_dispersion": params[1],
                "case_mean_multiplier": params[2],
                "batch_alignment": params[3],
                "minimum_tested_n_per_group": min(eligible) if eligible else None,
            }
        )
    return lookups


def render(report, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    rows = [
        r
        for r in report["rows"]
        if r["case_mean_multiplier"] == 1.5
        and r["nb_dispersion"] == 0.1
        and r["batch_alignment"] == 0.0
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    keys = ("null_raw", "null_corrected", "alternative_corrected")
    titles = ("Raw null rejection", "Adjusted null rejection", "Adjusted alternative detection")
    levels = np.linspace(0, 1, 11)
    for ax, key, title in zip(axes, keys, titles, strict=True):
        table = np.array(
            [
                [
                    next(
                        r[key]
                        for r in rows
                        if r["composition_setting"] == s and r["n_per_group"] == n
                    )
                    for n in NGRID
                ]
                for s in COMPOSITION
            ]
        )
        cf = ax.contourf(NGRID, COMPOSITION, table, levels=levels, cmap="viridis", vmin=0, vmax=1)
        xx, yy = np.meshgrid(NGRID, COMPOSITION)
        ax.scatter(xx, yy, c="white", s=8, linewidths=0.2, edgecolor="#272e45")
        ax.set(
            title=title,
            xlabel="Patients per group",
            ylabel="Composition setting",
            xticks=NGRID,
            yticks=COMPOSITION,
        )
    cbar = fig.colorbar(cf, ax=axes, shrink=0.83)
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(1))
    fig.suptitle("Synthetic operating characteristics", fontsize=14)
    output.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context({"svg.fonttype": "none", "svg.hashsalt": "spatial-audit"}):
        fig.savefig(output, dpi=300, metadata={"Date": None})
    plt.close(fig)
    with output.with_name("sample_size_lookup.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        table_writer = csv.DictWriter(handle, fieldnames=list(report["minimum_n"][0]))
        table_writer.writeheader()
        table_writer.writerows(report["minimum_n"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates", type=int, default=96)
    p.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    p.add_argument("--seed", type=int, default=45019)
    p.add_argument("--out", type=Path, default=ROOT / "results/synthetic/phase_diagram.json")
    p.add_argument("--plot", type=Path, default=ROOT / "docs/data/synthetic/phase_diagram.svg")
    args = p.parse_args()
    if args.replicates < 1 or args.workers < 1 or args.seed < 0:
        p.error("replicates and workers must be positive; seed must be nonnegative")
    params = list(product(NGRID, COMPOSITION, DISPERSION, MEAN_SHIFT, BATCH_ALIGNMENT))
    jobs = [(i, v, args.replicates, args.seed) for i, v in enumerate(params)]
    if args.workers == 1:
        rows = list(map(simulate, jobs))
    else:
        with ProcessPoolExecutor(args.workers) as ex:
            rows = list(ex.map(simulate, jobs, chunksize=3))
    report = {
        "kind": "synthetic",
        "record_type": "phase_diagram",
        "seed": args.seed,
        "scenarios": len(rows),
        "replicates_per_state": args.replicates,
        "simulated_cohorts": 2 * len(rows) * args.replicates,
        "design": "Two independent groups, one feature, 6 by 6 occupied volumes per donor, "
        "50 micrometres pitch, 250 micrometres fine band. Gamma-Poisson counts. "
        "Both effect states simulate independent cohorts. Composition setting controls "
        "amplitude, and is not the realized explained energy. All rates are per-cohort.",
        "inference": "Signed, untruncated nugget-adjusted energy; matching one- and two-step path support. "
        "Two-sided absolute difference in equal-patient mean energy. Enumerate up to "
        "512 allocations, otherwise use 511 random allocations and add one. Within "
        "each scenario the Monte Carlo allocations are shared, so Wilson intervals "
        "condition on that reference. Unrestricted tests at nonzero batch alignment "
        "are sensitivity experiments. The intervals exclude uncertainty in model choice.",
        "lookup_rule": "Smallest tested n whose lower 95% Wilson sensitivity bound reaches 0.8 and "
        "upper null-rejection bound is at most 0.1. Complete aliasing excluded. "
        "Null means no tested n qualifies. Interpolation in the figure is visual.",
        "figure_slice": {"nb_dispersion": 0.1, "case_mean_multiplier": 1.5, "batch_alignment": 0.0},
        "rows": rows,
        "minimum_n": minimum_n(rows),
    }
    write_json(args.out, report)
    render(strict_read(args.out), args.plot)
    print(f"{report['simulated_cohorts']} synthetic cohorts, {len(rows)} scenarios")


if __name__ == "__main__":
    main()
