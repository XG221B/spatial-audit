"""Independent spectral identities, pathological inputs and regression controls."""

import json
from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse

from spatial_audit.fixtures import trap_cohort
from spatial_audit.schema import Cohort, Section
from spatial_audit.spectral_audit import (
    basis,
    build_grid,
    coefficients,
    component_support,
    filter_bands,
    nugget_variance,
    section_statistics,
    spectral_audit,
    trace_gain,
)


def lattice(n):
    y, x = np.mgrid[:n, :n]
    return np.column_stack((x.ravel(), y.ravel())).astype(float)


@pytest.mark.parametrize("n", [2, 3, 4, 6, 9])
@pytest.mark.parametrize("um", [0.4, 1.0, 12.5])
def test_neumann_flux_and_spectrum(n, um):
    grid = build_grid(lattice(n) * um, 3 * um)
    L = grid.laplacian.toarray()
    assert np.allclose(L, L.T)
    assert np.allclose(L.sum(axis=1), 0)
    eig = np.linalg.eigvalsh(L)
    assert eig.min() >= -1e-10 / um**2
    assert eig.max() <= 8 / um**2
    assert np.allclose(eig[1], (2 - 2 * np.cos(np.pi / n)) / um**2)
    assert grid.left.size == 2 * n * (n - 1)


@pytest.mark.parametrize("mode", [0, 1, 2, 3, 6, 8])
@pytest.mark.parametrize("scale", [0.8, 1.4, 3.0])
def test_transfer_against_known_modes(mode, scale):
    xy = lattice(9)
    waves = np.array([3.0, 8.0, 17.0]) * scale
    grid = build_grid(xy, 3.0)
    f = np.cos(np.pi * mode * (xy[:, 0] + 0.5) / 9)[:, None]
    freq = 2 - 2 * np.cos(np.pi * mode / 9)
    time = np.log(4) / 3 * (waves / (2 * np.pi)) ** 2
    outputs = filter_bands(grid, f, coefficients(waves, grid.width))
    expected = np.exp(-time * freq) - np.exp(-4 * time * freq)
    for out, gain in zip(outputs, expected, strict=True):
        np.testing.assert_allclose(out, gain * f, atol=5e-9)


@pytest.mark.parametrize("shift", [0, 1, -31.0, 110.0])
@pytest.mark.parametrize("factor", [0.1, 1, 150.0])
def test_physical_similarity_and_translation(shift, factor):
    xy = lattice(8)
    f = np.random.default_rng(11).normal(size=(64, 2))
    waves = np.array([4.0, 10.0, 20.0])
    base = build_grid(xy, waves[0])
    changed = build_grid((xy + shift) * factor, waves[0] * factor)
    obs = filter_bands(changed, f, coefficients(waves * factor, changed.width))
    ref = filter_bands(base, f, coefficients(waves, base.width))
    np.testing.assert_allclose(obs, ref, atol=2e-10)


@pytest.mark.parametrize("n", [2, 3, 4, 7])
def test_disconnected_isolated_constant(n):
    coords = np.r_[lattice(n), lattice(n) + 50, [[100, 100]]]
    grid = build_grid(coords, 4)
    x = np.repeat([2.0, -5.0, 12.0], [n * n, n * n, 1])[:, None]
    f = filter_bands(grid, grid.average(x), coefficients([4, 10, 20], grid.width))
    assert np.max(np.abs(f)) < 1e-10


@pytest.mark.parametrize("p", [1, 2, 5])
def test_exact_projected_white_gain(p):
    coords = lattice(7)
    grid = build_grid(coords, 4)
    weights = coefficients([4, 8, 12], grid.width)
    rng = np.random.default_rng(37)
    q = np.linalg.qr(rng.normal(size=(len(coords), p)))[0]
    eig, u = np.linalg.eigh(grid.laplacian.toarray())
    t = np.log(4) / 3 * (4 / (2 * np.pi)) ** 2
    h = (u * (np.exp(-t * eig) - np.exp(-4 * t * eig))) @ u.T
    truth = np.linalg.norm(h - q @ (q.T @ h)) ** 2 / len(coords)
    assert trace_gain(grid, weights, q) == pytest.approx(truth, abs=1e-9)


@pytest.mark.parametrize("factor", [1, 1e-8, 1e8])
def test_singular_covariates_and_column_units(factor):
    u = np.random.default_rng(4).normal(size=(40, 2))
    x = np.column_stack([u[:, 0] * factor, u[:, 1], u.sum(1), np.zeros(40)])
    q = basis(x, x)
    assert q.shape == (40, 2)
    np.testing.assert_allclose(q @ (q.T @ x), x, atol=1e-6)
    np.testing.assert_allclose(q.T @ q, np.eye(2), atol=1e-12)


@pytest.mark.parametrize(
    "waves",
    [
        [1, 1, 2],
        [-1, 2, 3],
        [0, 2, 3],
        [1, np.nan, 3],
        [2, 1, 3],
        [],
        [[1, 3]],
        [1, 2, np.inf],
        [1, 2, 1e9],
    ],
)
def test_invalid_bands(waves):
    with pytest.raises(ValueError):
        coefficients(waves, 1)


@pytest.mark.parametrize("raw", [True, False])
def test_zero_singular_and_no_covariates(raw):
    coords = lattice(5)
    X = np.zeros((25, 3))
    X = X if raw else sparse.csr_matrix(X)
    sections = [Section(str(i), str(i // 2), coords, X, C=np.ones((25, 2))) for i in range(4)]
    ct = Cohort("zero", "synthetic", ["a", "b", "c"], sections)
    composition, nugget = spectral_audit(ct, (3, 6, 10))
    assert all(b["share"] is None for b in composition["bands"])
    assert nugget["flagged"] == 0
    for section in sections:
        section.C = None
    assert spectral_audit(ct, (3, 6, 10))[0]["status"] == "yellow"


def test_permutation_sparse_and_basis_invariance():
    coords = lattice(8)
    rng = np.random.default_rng(124)
    C = rng.normal(size=(64, 2))
    X = rng.normal(size=(64, 7)) + C[:, [0]]
    order = rng.permutation(len(X))
    s = Section("p", "a", coords, X, C=C)
    t = replace(
        s,
        coords=coords[order],
        X=sparse.csr_matrix(X[order]),
        C=(C @ [[2.0, 1.0], [1.0, 1.0]])[order],
    )
    before, _ = section_statistics(s, (4.0, 8.0, 16.0), block=2)
    after, _ = section_statistics(t, (4.0, 8.0, 16.0), block=1)
    np.testing.assert_allclose(after, before, atol=1e-9)


def test_two_lag_variogram_local_linear_field():
    grid = build_grid(lattice(6), 3)
    x = grid.average(lattice(6))
    np.testing.assert_allclose(nugget_variance(grid, x), 0, atol=1e-12)


def test_white_noise_nugget_estimation():
    grid = build_grid(lattice(32), 4)
    rng = np.random.default_rng(808)
    x = rng.normal(size=(1024, 64)) * 3
    assert np.mean(nugget_variance(grid, x)) == pytest.approx(9.0, rel=0.04)


def test_sample_repeat_weighting_and_observation_order():
    ct = trap_cohort("composition")
    prior = spectral_audit(ct)
    ct.sections = ct.sections + [replace(ct.sections[0], section_id="duplicate")]
    current = spectral_audit(ct)
    for before, after in zip(prior[0]["bands"], current[0]["bands"], strict=True):
        assert before["share"] == after["share"]
        assert before["supported_patients"] == after["supported_patients"]
        assert before["supported_sections"] + 1 == after["supported_sections"]
    assert prior[1]["features"] == current[1]["features"]


@pytest.mark.parametrize(
    ("trap", "share", "hf"),
    [
        ("null", "green", "green"),
        ("pseudorep", "green", "green"),
        ("composition", "red", "red"),
        ("stage", "green", "green"),
    ],
)
def test_prespecified_spectral_controls(trap, share, hf):
    comp, nugget = spectral_audit(trap_cohort(trap))
    assert comp["status"] == share
    assert nugget["status"] == hf
    assert all(0 <= b["share"] <= 1 for b in comp["bands"])


def test_unresolved_fine_band_and_unknown_nugget():
    s = Section(
        "a", "one", np.array([[0, 0], [300, 0]]), np.array([[0.0], [1.0]]), C=np.zeros((2, 1))
    )
    ct = Cohort("tiny", "synthetic", ["a"], [s, replace(s, patient_id="b", group="two")])
    comp, noise = spectral_audit(ct, (2, 3, 4))
    assert comp["status"] == "yellow"
    assert noise["status"] == "yellow"
    assert all(row["corrected"] == [None] for row in noise["patient_energies"])


def diagnostic_cohort(coords, cov=None, seed=14):
    rng = np.random.default_rng(seed)
    rows = [
        Section(
            str(i),
            "A" if i < 2 else "B",
            coords,
            rng.normal(size=(len(coords), 2)),
            C=cov,
            section_id=str(i),
        )
        for i in range(4)
    ]
    return Cohort("geometry-test", "synthetic", ["F1", "F2"], rows)


def test_small_field_returns_null_at_unsupported_wavelengths():
    coords = lattice(18) * 30
    result, nugget = spectral_audit(diagnostic_cohort(coords, coords[:, 0:1]))
    assert [b["supported_patients"] for b in result["bands"]] == [4, 0, 0]
    assert result["bands"][0]["share"] is not None
    assert all(b["share"] is None for b in result["bands"][1:])
    assert result["status"] == "yellow"
    json.dumps([result, nugget], allow_nan=False)


def test_disjoint_islands_have_no_joint_physical_span():
    xy = np.concatenate([lattice(4), lattice(4) + 20])
    grid = build_grid(xy, 4)
    masks, support, count, span = component_support(grid, xy, (4, 8, 12))
    assert count == 2 and span == pytest.approx(3)
    assert support == [0, 0, 0] and not np.any(masks)


def test_diagonal_ribbon_cannot_supply_isotropic_support():
    x = np.arange(30)
    xy = np.r_[np.column_stack((x, x)), np.column_stack((x + 1, x))]
    grid = build_grid(xy, 4)
    _, support, count, span = component_support(grid, xy, (4, 8, 12))
    assert count == 1
    assert span < 4
    assert support == [0, 0, 0]


def test_axial_diffusion_limit_and_exact_span_boundary():
    xy = lattice(6) * 50
    grid = build_grid(xy, 120)
    _, support, _, _ = component_support(grid, xy, (120, 50 * np.pi, 250))
    assert support == [0, 1, 1]


def test_six_by_six_phase_diagram_footprint_at_fine_boundary():
    xy = lattice(6) * 50
    sec = Section("P", "A", xy, xy.copy())
    vals, info = section_statistics(sec, (250.0, 600.0, 1400.0))
    assert info["bands"][0]["supported_volumes"] == 36
    assert np.isfinite(vals[7]).all()
    assert all(b["supported_volumes"] == 0 for b in info["bands"][1:])


def test_missing_support_reports_actual_patients_and_sections():
    xy = lattice(16)
    cohort = diagnostic_cohort(xy, xy[:, 0:1])
    s = cohort.sections[1]
    cohort.sections[1] = replace(s, coords=s.coords[:16], C=s.C[:16], X=s.X[:16])
    result, noise = spectral_audit(cohort, (4.0, 8.0, 12.0))
    for b in result["bands"]:
        assert b["supported_patients_by_group"] == {"A": 1, "B": 2}
        assert b["supported_sections"] == 3
        assert b["share"] is not None and not b["full_coverage"]
    assert result["status"] == "yellow"
    assert noise["patient_energies"][1]["corrected"] == [None, None]
    assert noise["features"][0]["debiased_delta"] is None
    assert noise["status"] == "yellow"
    json.dumps([result, noise], allow_nan=False)


def test_saturated_projection_has_no_interpretable_share():
    xy = lattice(6)
    cohort = diagnostic_cohort(xy, np.eye(len(xy)))
    composition, nugget = spectral_audit(cohort, (3.2, 4.0, 5.0))
    for diag in composition["sections"]:
        for b in diag["bands"]:
            assert b["saturated"]
            assert b["projected_white_fraction"] == pytest.approx(1)
            assert 1 <= b["effective_spectral_dimension"] <= len(xy) - 1
    assert all(b["share"] is None for b in composition["bands"])
    assert composition["status"] == nugget["status"] == "yellow"
    json.dumps([composition, nugget], allow_nan=False)


def test_effective_dimension_against_exact_eigenvalues():
    xy = lattice(8)
    _, info = section_statistics(Section("p", "A", xy, xy.copy()), (4.0, 5.0, 6.0))
    eig = np.linalg.eigvalsh(build_grid(xy, 4).laplacian.toarray())
    for record, wave in zip(info["bands"], [4.0, 5.0, 6.0], strict=True):
        t = np.log(4) / 3 * (wave / (2 * np.pi)) ** 2
        gain = np.exp(-t * eig) - np.exp(-4 * t * eig)
        expected = np.sum(gain**2) ** 2 / np.sum(gain**4)
        assert record["effective_spectral_dimension"] == pytest.approx(expected, abs=1e-8)


def test_matched_lag_intercept_against_second_difference_on_uneven_islands():
    xy = np.r_[lattice(7), lattice(3) + 18, [[25.0, 31.0]]]
    grid = build_grid(xy, 4)
    x = np.random.default_rng(17).normal(size=(len(grid.count), 3))
    masks = component_support(grid, xy, [4.0, 8.0, 12.0])[0]
    mask = masks[0]
    edges = mask[grid.far_left] & mask[grid.middle] & mask[grid.far_right]
    a, m, b = (idx[edges] for idx in (grid.far_left, grid.middle, grid.far_right))
    lag1 = np.mean((x[a] - x[m]) ** 2 + (x[m] - x[b]) ** 2, axis=0) / 4
    lag2 = np.mean((x[b] - x[a]) ** 2, axis=0) / 2
    np.testing.assert_allclose(nugget_variance(grid, x, mask), (4 * lag1 - lag2) / 3)


def test_signed_negative_energy_reaches_patient_inference():
    xy = lattice(8)
    cohort = diagnostic_cohort(xy)
    for s in cohort.sections:
        amp = 1 if s.group == "A" else 2
        checkerboard = (-1.0) ** (xy[:, 0] + xy[:, 1])
        s.X = amp * np.column_stack([checkerboard, -checkerboard])
    result, nugget = spectral_audit(cohort, (4.0, 5.0, 6.0))
    assert all(v < 0 for row in nugget["patient_energies"] for v in row["corrected"])
    energies = nugget["features"][0]["corrected_energy"]
    assert energies[1] == pytest.approx(4 * energies[0])
    assert nugget["features"][0]["debiased_delta"] == pytest.approx(energies[1] - energies[0])
    json.dumps([result, nugget], allow_nan=False)


def test_constant_intensity_offset_preserves_band_energies():
    xy = lattice(8)
    base = Section("p", "A", xy, np.column_stack([xy[:, 0] ** 2, xy[:, 1]]))
    shifted = replace(base, X=base.X + 1e14)
    values, _ = section_statistics(base, (4.0, 5.0, 6.0))
    changed, _ = section_statistics(shifted, (4.0, 5.0, 6.0))
    np.testing.assert_allclose(changed[[0, 1, 2, 6, 7, 9]], values[[0, 1, 2, 6, 7, 9]], atol=5e-8)
