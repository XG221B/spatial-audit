"""Paired raw/Anscombe calibration on the zero-alignment physical lattice."""

import argparse
from itertools import combinations, product
from pathlib import Path

import numpy as np
from run_phase_diagram import operators, wilson

from spatial_audit.schema import Section, write_json
from spatial_audit.spectral_audit import section_statistics

ROOT = Path(__file__).resolve().parents[1]


def run(reps=192, seed=70401):
    rng = np.random.default_rng(seed)
    _, _, _, _, c, signal = operators()
    y, x = np.mgrid[:6, :6]
    coords = np.column_stack((x.ravel(), y.ravel())) * 50.0
    rows = []
    for n, shift, disp in product((4, 5, 7), (1.0, 1.5, 2.0), (0.0, 0.1, 0.35)):
        group = np.repeat([0, 1], n)
        allocations = np.array(list(combinations(range(2 * n), n)))
        masks = np.full((len(allocations), 2 * n), -1.0 / n)
        masks[np.arange(len(masks))[:, None], allocations] = 1.0 / n
        for state, effect in (("null", 0.0), ("alternative", 7.0)):
            batch = rng.binomial(1, 0.5, (reps, 2 * n))
            amplitude = 4 + rng.normal(size=batch.shape) + 2.5 * batch + effect * group
            mu = np.maximum(
                1,
                (25 * (1 + (shift - 1) * group))[:, None]
                + amplitude[:, :, None] * signal
                + 4 * np.sqrt(0.25 / 0.75) * (0.7 + group[:, None]) * c,
            )
            counts = rng.poisson(rng.gamma(1 / disp, mu * disp) if disp else mu)
            # Each feature column is one donor realization on the identical lattice.
            # All operations through signed energy are columnwise in the audit engine.
            section = Section("simulation", state, coords, counts.reshape(-1, 36).T, c[:, None])
            for mode in ("raw", "vst"):
                stats, _ = section_statistics(
                    section, (250.0, 600.0, 1400.0), block=128, nugget_mode=mode
                )
                energy = stats[-3].reshape(reps, 2 * n)
                if not np.isfinite(energy).all():
                    raise ValueError("Incomplete lattice support")
                observed = abs(energy[:, :n].mean(axis=1) - energy[:, n:].mean(axis=1))
                reference = abs(energy @ masks.T)
                p = (reference >= (observed - 1e-11 * np.maximum(1, observed))[:, None]).mean(
                    axis=1
                )
                calls = int(np.sum(p < 0.05))
                rows.append(
                    {
                        "n_per_arm": n,
                        "case_mean_multiplier": shift,
                        "nb_dispersion": disp,
                        "state": state,
                        "nugget_mode": mode,
                        "allocations": len(masks),
                        "rejections": calls,
                        "rate": calls / reps,
                        "wilson95": wilson(calls / reps, reps),
                    }
                )
    return {
        "kind": "synthetic",
        "record_type": "vst_calibration",
        "seed": seed,
        "replicates_per_state": reps,
        "batch_alignment": 0,
        "composition_setting": 0.25,
        "alpha": 0.05,
        "design": "Independent 6x6 sections, one site per 50 um lattice position, 250 um fine band. "
        "Same count draws in both modes. Equal latent raw-scale amplitude distribution under null; "
        "alternative adds 7 raw count units of structured amplitude in cases. Independent Bernoulli(0.5) "
        "batch shifts amplitude in both groups. Poisson or gamma-Poisson observations.",
        "inference": "All balanced donor allocations, absolute difference in mean signed corrected "
        "fine-band energy. Strict p < 0.05. Wilson intervals describe replicate sampling uncertainty.",
        "rows": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=192)
    parser.add_argument("--seed", type=int, default=70401)
    args = parser.parse_args()
    if args.replicates < 1 or args.seed < 0:
        parser.error("invalid simulation budget")
    write_json(ROOT / "results/synthetic/vst_calibration.json", run(args.replicates, args.seed))
