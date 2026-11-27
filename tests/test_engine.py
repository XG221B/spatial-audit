from dataclasses import replace
from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy import sparse

from spatial_audit.design_audit import design_audit, exact_floor, power_audit
from spatial_audit.fixtures import fixtures
from spatial_audit.pseudorep import bh, patient_test, pooled_rank_p, pseudorep_audit
from spatial_audit.report import audit
from spatial_audit.schema import Cohort, Section


@pytest.mark.parametrize("n0,n1", [(i, j) for i in (1, 2, 3, 6, 12) for j in (2, 3, 6, 19)])
def test_double_tail_floor(n0, n1):
    assert exact_floor(n0, n1)["p_floor"] == pytest.approx(min(1, 2 / comb(n0 + n1, n0)))


@pytest.mark.parametrize(
    "n0,n1,m,a",
    [
        (0, 1, 1, 0.05),
        (1, -2, 1, 0.05),
        (True, 1, 1, 0.05),
        (1, 2, 0, 0.05),
        (1, 2, 1, 1),
        (1, 2, 1, 0),
    ],
)
def test_invalid_design(n0, n1, m, a):
    with pytest.raises(ValueError):
        exact_floor(n0, n1, m, a)


def test_bh_rank_not_bonferroni():
    d = exact_floor(6, 6, 20000)
    assert 1 < d["minimum_bh_rank"] < 20000 and not d["bh_impossible"]
    assert exact_floor(3, 3, 20000)["bh_impossible"]
    assert exact_floor(1, 39, 1)["minimum_bh_rank"] is None


def toy():
    coords = np.array([[0.0, 0], [1, 0], [0, 1], [1, 1]])
    return Cohort(
        "toy",
        "synthetic",
        ["f"],
        [
            Section(
                str(i),
                str(i % 2),
                coords,
                np.arange(4.0)[:, None] + i,
                slide_id="same",
                batch="one",
                stage=(i // 2, 0),
            )
            for i in range(6)
        ],
    )


def test_independent_patients():
    c = toy()
    c.sections *= 2
    assert power_audit(c)["assignments"] == "20"
    assert design_audit(c)["status"] == "green"


def test_categorical_alias():
    c = toy()
    for s in c.sections:
        s.batch = s.group
    assert design_audit(c)["reasons"] == ["batch"]


def test_missing_metadata():
    c = toy()
    c.sections[0].slide_id = None
    assert design_audit(c)["status"] == "yellow"


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int64])
def test_stereo_constant_feature_with_large_unequal_groups(dtype):
    a, b = np.zeros(19781, dtype=dtype), np.zeros(58508, dtype=dtype)
    assert pooled_rank_p(a, b) == 1.0
    assert list(bh([0.001, pooled_rank_p(a, b)])) == pytest.approx([0.002, 1.0])
    b[::3] = 1
    assert pooled_rank_p(a, b) < 1e-10


def test_bh_fails_before_invalid_family_can_be_published():
    with pytest.raises(ValueError, match="finite p"):
        bh([0.01, np.nan])


def test_eb_invariant_to_feature_specific_units():
    rng = np.random.default_rng(512)
    y = rng.normal(size=(8, 5))
    groups = np.repeat([0, 1], 4)
    original = patient_test(y, groups)
    transformed = patient_test(y * [1, 1e3, 1e-4, 7, 3] + [10, 2, 0, 1, 8], groups)
    np.testing.assert_allclose(transformed["statistic"], original["statistic"])
    assert transformed["p"] == original["p"]


def test_effective_resolution_reflects_monte_carlo_budget():
    c = toy()
    c.features = ["f" + str(i) for i in range(2500)]
    for j in range(12):
        c.sections.append(replace(c.sections[j % 6], patient_id=f"other{j}"))
    result = power_audit(c, max_exact=1, draws=4095)
    assert result["effective_p_floor"] == 2 / 4096
    assert result["effective_minimum_bh_rank"] == 25
    assert result["minimum_bh_rank"] < result["effective_minimum_bh_rank"]


def test_joint_alias_requires_observed_acquisition_metadata():
    c = toy()
    for section in c.sections:
        section.slide_id = section.batch = section.stage = None
    assert design_audit(c)["joint_acquisition_design"]["fixed_effect_alias"] is None


def test_patient_positive_control_survives_bh():
    rng = np.random.default_rng(117)
    group = np.repeat([0, 1], 6)
    values = rng.normal(0, 0.1, size=(12, 8))
    values[:, :3] += 5 * group[:, None]
    assert np.all(np.asarray(patient_test(values, group)["q"])[:3] < 0.05)


def test_joint_acquisition_alias_can_escape_marginal_screens():
    c = toy()
    for section, (slide, batch) in zip(c.sections, [(0, 0), (1, 0), (0, 1)] * 2):
        section.slide_id, section.batch = str(slide), str(batch)
        section.group, section.stage = str(slide + batch), None
    result = design_audit(c)
    assert not any(d["complete_alias"] for d in result["categorical"].values())
    assert result["reasons"] == ["joint_acquisition_main_effects"]
    assert result["joint_acquisition_design"]["relative_residual"] < 1e-9


def test_observed_stage_separation_need_not_alias_a_linear_trend():
    c = toy()
    for i, section in enumerate(c.sections):
        section.stage = (float(i), 0.0)
        section.group = str(int(i >= 3))
    result = design_audit(c)
    assert result["coordinates"]["stage"]["linear_separation"]
    assert not result["coordinates"]["stage"]["linear_fixed_effect_alias"]
    assert not result["joint_acquisition_design"]["fixed_effect_alias"]


@pytest.mark.parametrize("size", [(2, 2), (2, 3), (3, 3), (2, 7)])
def test_oracle_and_symmetry(size):
    n0, n1 = size
    rng = np.random.default_rng(19)
    y = rng.normal(size=n0 + n1)
    labels = np.r_[np.zeros(n0), np.ones(n1)]
    perm = patient_test(y[:, None], labels)
    diffs = np.array(
        [
            y[list(idx)].mean() - np.delete(y, list(idx)).mean()
            for idx in combinations(range(n0 + n1), n1)
        ]
    )
    obs = y[-n1:].mean() - y[:n0].mean()
    p = min(1, 2 * min(np.mean(diffs >= obs - 1e-12), np.mean(diffs <= obs + 1e-12)))
    assert perm["p"] == pytest.approx([p])
    assert patient_test(y[:, None], 1 - labels)["p"] == pytest.approx(perm["p"])


def test_monte_carlo_and_zero_features():
    v = np.zeros((12, 3))
    p = patient_test(v, [0] * 3 + [1] * 9, max_exact=1, draws=9)
    assert p["method"] == "monte_carlo_plus_one" and p["q"] == [1, 1, 1]


def test_sparse_pooled_result():
    c = toy()
    original = pseudorep_audit(c)
    for s in c.sections:
        s.X = sparse.csr_matrix(s.X)
    assert original == pseudorep_audit(c)


@pytest.mark.parametrize(
    "c,expected",
    list(zip(fixtures(), [set(), {"pseudorep"}, {"composition", "nugget"}, {"design"}])),
)
def test_fixture_alarm_vector(c, expected):
    d = audit(c)["dimensions"]
    assert {k for k, v in d.items() if v["status"] != "green"} == expected
