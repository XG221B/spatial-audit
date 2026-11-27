"""Calibration checks for the vectorized lattice experiment and its published lookup."""

import runpy
from pathlib import Path

import numpy as np
import pytest

from spatial_audit.schema import Section, strict_read
from spatial_audit.spectral_audit import section_statistics

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def simulator():
    return runpy.run_path(str(ROOT / "scripts/run_phase_diagram.py"))


def test_reduced_vectorized_statistics_match_audit(simulator):
    rng = np.random.default_rng(719)
    X = rng.poisson(20.0, size=(2, 3, 36))
    raw, corr = simulator["energy"](X)
    _, _, _, _, c, _ = simulator["operators"]()
    coords = np.indices((6, 6))[::-1].reshape(2, -1).T * 50.0
    for arm in range(2):
        for patient in range(3):
            s = Section("P", "A", coords, X[arm, patient, :, None], C=c[:, None])
            stats, _ = section_statistics(s, (250.0, 600.0, 1400.0))
            assert raw[arm, patient] == pytest.approx(stats[0, 0], abs=1e-8)
            assert corr[arm, patient] == pytest.approx(stats[7, 0], abs=1e-8)


@pytest.mark.parametrize("n", [3, 5, 8, 30])
def test_allocations_are_balanced_and_centered(simulator, n):
    masks, exact = simulator["assignments"](n, np.random.default_rng(818))
    assert np.max(abs(masks.sum(axis=1))) < 1e-12
    assert ((masks > 0).sum(axis=1) == n).all()
    if exact:
        assert len(np.unique(masks, axis=0)) == len(masks)
    assert simulator["rejection_rate"](np.ones((2, 2 * n)), masks, exact, n) == 0


def test_published_scenarios_and_lookup():
    data = strict_read(ROOT / "results/synthetic/phase_diagram.json")
    assert data["kind"] == "synthetic"
    assert data["simulated_cohorts"] == 2 * data["replicates_per_state"] * len(data["rows"])
    assert data["scenarios"] == len(data["rows"]) >= 10000
    keys = [
        "n_per_group",
        "composition_setting",
        "nb_dispersion",
        "case_mean_multiplier",
        "batch_alignment",
    ]
    unique = {tuple(row[k] for k in keys) for row in data["rows"]}
    assert len(unique) == data["scenarios"]
    for row in data["rows"]:
        for state in ("null", "alternative"):
            lower, upper = row[state + "_interval"]
            assert 0 <= lower <= row[state + "_corrected"] + 1e-12 <= upper + 1e-10 <= 1.0 + 1e-9
            if row["n_per_group"] == 3:
                assert row[state + "_corrected"] == row[state + "_raw"] == 0.0
    assert all(
        row["minimum_tested_n_per_group"] is None
        for row in data["minimum_n"]
        if row["batch_alignment"] == 1.0
    )


def test_memory_record_uses_os_process_peak():
    record = strict_read(ROOT / "results/synthetic/memory.json")
    assert record["kind"] == "synthetic" and len(record["dimensions"]) == 5
    assert record["n_observations"] >= 1_000_000 and record["nnz"] > 0
    assert 0 < record["peak_rss_bytes"] < record["budget_bytes"] <= 1.5 * 2**30
    assert record["os_counter"] in {"Windows PeakWorkingSetSize", "getrusage RUSAGE_SELF ru_maxrss"}
