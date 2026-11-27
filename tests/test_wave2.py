"""Cross-check sufficient statistics and denominators of the validation experiments."""

import sys
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse

from spatial_audit.pseudorep import pooled_rank_p
from spatial_audit.schema import strict_read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_empirical_experiments import panel_sensitivity, pooled_p, rank_summaries


def test_midrank_aggregation_matches_scipy_on_unequal_tied_sections():
    rng = np.random.default_rng(123)
    x = [np.column_stack((rng.poisson(0.3, (n, 4)), np.ones(n))) for n in (14, 7, 31, 19)]
    sizes, ranks, variance = rank_summaries([SimpleNamespace(X=sparse.csr_matrix(a)) for a in x])
    for selected in combinations(range(4), 2):
        mask = np.isin(range(4), selected)
        a, b = (
            np.concatenate([row for row, keep in zip(x, label) if keep]) for label in (mask, ~mask)
        )
        expected = [pooled_rank_p(col_a, col_b) for col_a, col_b in zip(a.T, b.T)]
        np.testing.assert_allclose(pooled_p(mask, sizes, ranks, variance), expected, atol=1e-12)


def test_checked_panel_and_empirical_denominators():
    assert panel_sensitivity() == strict_read(ROOT / "results/public/panel_sensitivity.json")
    null = strict_read(ROOT / "results/public/empirical_null.json")
    for arm in null["arms"]:
        assert len(arm["rows"]) == null["splits_per_arm"]
        for row in arm["rows"]:
            assert len(set(row["arm_a"])) == arm["group_sizes"][0]
            assert set(row["arm_a"]) <= set(arm["donors"])
        for key, summary in arm["summary"].items():
            assert summary["bh_calls"] == sum(r["bh_calls"][key] for r in arm["rows"])
            assert np.isclose(summary["bh_fraction"], summary["bh_calls"] / summary["tests"])
