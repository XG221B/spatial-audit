"""Numerical and end-to-end contracts for the optional study designs and spatial estimators."""

from dataclasses import replace
from itertools import product
from pathlib import Path
from runpy import run_path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import fisher_exact

from spatial_audit import (
    audit,
    audit_anndata,
    compare_methods,
    load_manifest,
    recommend_wavelengths,
)
from spatial_audit.cli import export_fixtures, main
from spatial_audit.design_audit import design_audit, exact_floor
from spatial_audit.pseudorep import patient_test, pseudorep_audit
from spatial_audit.schema import Cohort, Section, strict_read, write_json
from spatial_audit.spectral_audit import (
    build_grid,
    build_voronoi,
    coefficients,
    filter_bands,
    section_statistics,
)

ROOT = Path(__file__).resolve().parents[1]


def lattice(n=8):
    return np.indices((n, n)).reshape(2, -1).T.astype(float)


def cohort(paired=False):
    xy = lattice()
    rng = np.random.default_rng(13)
    sections = [
        Section(
            str(i) if paired else f"{g}{i}",
            g,
            xy,
            rng.poisson(6 + 4 * j, size=(len(xy), 2)).astype(float),
            C=xy[:, :1],
            stage=(i + 6 * j, 0),
            batch=str(i),
            slide_id="shared",
            section_id=f"s-{g}{i}",
        )
        for i in range(4)
        for j, g in enumerate("AB")
    ]
    return Cohort(
        "example", "synthetic", ["x", "y"], sections, design="paired" if paired else "independent"
    )


@pytest.mark.parametrize(
    ("budget", "numerators", "denominator"),
    [(70, [6, 30, 20, 6, 4], 70), (1, [2, 18, 10, 2, 6], 38)],
)
def test_batched_matches_preupgrade_scalar_test(budget, numerators, denominator):
    y = np.random.default_rng(7).normal(size=(8, 5))
    test = patient_test(y, [0] * 4 + [1] * 4, max_exact=budget, draws=37)
    np.testing.assert_array_equal(test["p"], np.asarray(numerators) / denominator)


@pytest.mark.parametrize("budget", [1, 16])
def test_paired_randomization_against_sign_flips(budget):
    y = np.array([[1, 4], [2, 0], [3, 9], [7, 6]], float)
    patients = np.repeat(np.arange(4), 2)
    labels = np.tile(["before", "after"], 4)
    draws = 43
    test = patient_test(
        y.reshape(-1, 1),
        labels,
        max_exact=budget,
        draws=draws,
        design="paired",
        patient_ids=patients,
    )
    delta = y[:, 1] - y[:, 0]
    if budget == 16:
        stats = np.array(list(product([-1, 1], repeat=4))) @ delta
        expected = np.mean(np.abs(stats) >= abs(delta.sum()))
    else:
        rng = np.random.default_rng(17)
        stats = [
            sum((2 * rng.choice([0, 1], 1, replace=False)[0] - 1) * x for x in delta)
            for _ in range(draws)
        ]
        # Sorted labels put "before" in the positive arm.
        expected = min(
            1,
            2
            * (
                1
                + min(sum(s >= -delta.sum() for s in stats), sum(s <= -delta.sum() for s in stats))
            )
            / (draws + 1),
        )
    assert test["p"] == pytest.approx([expected])
    assert test["assignments"] == "16"
    assert exact_floor(4, 4, design="paired")["p_floor"] == 2 / 16


def test_paired_validation_repetition_and_all_endpoints():
    c = cohort(True)
    c.sections.append(c.sections[0])
    r = audit(c, wavelengths=[4], alpha=0.2)
    assert r["n_patients"] == 4
    d = r["dimensions"]
    assert d["power"]["alpha"] == d["pseudorep"]["alpha"] == 0.2
    for test in (d["pseudorep"]["patient"], d["nugget"]["raw_test"], d["nugget"]["corrected_test"]):
        assert test["assignments"] == "16" and test["design"] == "paired"
    c.sections.pop(1)
    with pytest.raises(ValueError, match="both conditions"):
        c.validate()


def test_stratified_randomization_block_counts_and_fixed_labels():
    labels = ["a", "b", "b", "a", "b", "a"]
    strata = [0, 0, 0, 1, 1, 1]
    y = np.array([0.0, 2, 3, 5, 6, 1])
    from itertools import combinations

    assignments = [
        list(a + b) for a, b in product(combinations(range(3), 2), combinations(range(3, 6), 1))
    ]
    sums = y[assignments].sum(1)
    observed = y[np.array(labels) == "b"].sum()
    expected = min(1, 2 * min(np.mean(sums >= observed), np.mean(sums <= observed)))
    result = patient_test(y[:, None], labels, design="stratified", strata=strata)
    assert result["p"] == pytest.approx([expected])
    assert result["assignments"] == "9"
    assert (
        exact_floor(3, 3, design="stratified", stratum_counts=[(1, 2), (2, 1)])["assignments"]
        == "9"
    )
    c = cohort()
    for s in c.sections:
        s.slide_id = s.group
    r = audit(c, wavelengths=[4], design="stratified", stratum="slide_id")
    assert r["dimensions"]["power"]["p_floor"] == 1
    assert r["dimensions"]["pseudorep"]["patient"]["p"] == [1, 1]
    c.sections.append(replace(c.sections[0], slide_id="elsewhere"))
    with pytest.raises(ValueError, match="one observed stratum"):
        audit(c, design="stratified", stratum="slide_id")


def test_design_centroids_lopo_multilevel_and_missingness():
    c = cohort()
    c.sections[0].slide_id = None
    d = design_audit(c)
    stage = d["coordinates"]["stage"]
    assert stage["patient_centroid_separation_um"] == pytest.approx(6)
    assert stage["lopo_cv_auc"] == pytest.approx(1)
    assert d["categorical"]["batch"]["permutation_p"] == 1
    slide = d["categorical"]["slide_id"]
    assert slide["missing_fraction"] == 1 / 8
    assert slide["missingness_group_p"] == fisher_exact([[1, 3], [0, 4]]).pvalue
    c.sections.extend([c.sections[2]] * 4)
    assert stage["lopo_cv_auc"] == design_audit(c)["coordinates"]["stage"]["lopo_cv_auc"]


def test_depth_normalization_removes_only_depth_shift():
    c = cohort()
    for s in c.sections:
        factor = 10 if s.group == "B" else 1
        s.X = sparse.csr_matrix(np.tile([1.0, 2.0], (len(s.coords), 1)) * factor)
        s.library_size = np.full(len(s.coords), 50 * factor)
    d = pseudorep_audit(c, alpha=0.1)
    assert d["group_depth_ratio"] == 10
    assert d["patient_rate"] == 1 and d["depth_normalized_patient_rate"] == 0
    for s in c.sections:
        s.library_size = None
    assert pseudorep_audit(c)["depth_source"] == "registered panel sum"
    c.sections[0].X.data[0] = -1
    assert pseudorep_audit(c)["depth_normalized_patient_rate"] is None


@pytest.mark.parametrize("waves", [(4,), (4, 5), (4, 5, 6, 7), (4, 4.5, 5, 6, 7)])
def test_arbitrary_bands_share_identical_fine_endpoint(waves):
    s = cohort().sections[0]
    values, info = section_statistics(s, waves)
    reference, _ = section_statistics(s, [4])
    assert len(info["bands"]) == len(waves)
    np.testing.assert_allclose(values[-4:], reference[-4:], atol=1e-10)


def test_vst_matches_explicit_anscombe_reference():
    s = cohort().sections[0]
    raw, _ = section_statistics(s, [4, 6])
    transformed, _ = section_statistics(replace(s, X=2 * np.sqrt(s.X + 3 / 8)), [4, 6])
    vst, _ = section_statistics(s, [4, 6], nugget_mode="vst")
    np.testing.assert_allclose(vst[:4], raw[:4])
    np.testing.assert_allclose(vst[-4:], transformed[-4:])
    with pytest.raises(ValueError, match="nonnegative"):
        section_statistics(replace(s, X=-s.X), [4], nugget_mode="vst")


@pytest.mark.parametrize("n", [8, 13])
def test_trace_standard_error_matches_probe_variance(n):
    xy = lattice(n)
    section = Section("x", "a", xy, xy)
    _, diag = section_statistics(section, [4])
    grid = build_grid(xy, 4)
    p = (
        np.eye(n * n)
        if n * n <= 144
        else np.random.default_rng(819).choice([-1.0, 1.0], (n * n, 48))
    )
    (h,) = filter_bands(grid, p, coefficients([4], 0.5), 0.5)
    estimates = np.mean(h * h, axis=0)
    se = 0 if n * n <= 144 else np.std(estimates, ddof=1) / np.sqrt(48)
    assert diag["bands"][0]["residual_white_gain_se"] == pytest.approx(se)
    assert diag["trace_probes"] == (0 if n * n <= 144 else 48)


def test_leave_one_ion_out_matches_individual_covariate_projections():
    s = cohort().sections[0]
    s.C = np.column_stack([s.X.sum(1), s.C])
    original = s.C.copy()
    v, info = section_statistics(s, [4, 6], leave_one_out_tic=True)
    for j in range(2):
        cov = s.C.copy()
        cov[:, 0] -= s.X[:, j]
        reference, ref_info = section_statistics(replace(s, X=s.X[:, j : j + 1], C=cov), [4, 6])
        np.testing.assert_allclose(v[:, j], reference[:, 0])
        assert info["feature_bands"][j] == ref_info["bands"]
    np.testing.assert_equal(s.C, original)


@pytest.mark.parametrize("irregular", [True, False])
def test_voronoi_rigid_motion_invariance(irregular):
    xy = lattice()
    if irregular:
        xy += np.random.default_rng(9).uniform(-0.2, 0.2, size=xy.shape)
    ang = 0.49
    rotate = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
    moved = xy @ rotate + [37, -11]
    graph = build_voronoi(xy)
    np.testing.assert_allclose(
        graph.laplacian.toarray(), build_voronoi(moved).laplacian.toarray(), atol=1e-9
    )
    np.testing.assert_allclose(graph.laplacian @ np.ones(len(xy)), 0, atol=1e-12)
    s = Section("p", "g", xy, xy * xy)
    v, a = section_statistics(s, [4], graph_mode="voronoi")
    v2, b = section_statistics(replace(s, coords=moved), [4], graph_mode="voronoi")
    np.testing.assert_allclose(v, v2, atol=1e-9)
    assert a["bands"][0]["supported_volumes"] == b["bands"][0]["supported_volumes"] == len(xy)
    assert a["largest_2d_span_um"] == pytest.approx(b["largest_2d_span_um"])


def test_auto_bands_enforce_physical_support():
    assert recommend_wavelengths(2, 2 * np.pi) == pytest.approx([2 * np.pi])
    for args in [(1, 2), (np.nan, 5), (1, -2)]:
        with pytest.raises(ValueError):
            recommend_wavelengths(*args)
    result = audit(cohort(), wavelengths="auto")
    assert [
        b["wavelength_um"] for b in result["dimensions"]["composition"]["bands"]
    ] == pytest.approx(recommend_wavelengths(1, 7))


def test_manifest_cli_and_markdown(tmp_path, capsys):
    c = cohort(True)
    export_fixtures([c], tmp_path)
    path = tmp_path / "manifest.json"
    spec = strict_read(path)
    spec["cohorts"][0]["design"] = "paired"
    write_json(path, spec)
    assert load_manifest(path)[0].design == "paired"
    main(["validate-manifest", "--manifest", str(path)])
    assert "example" in capsys.readouterr().out
    main(
        [
            "run",
            "--manifest",
            str(path),
            "--out",
            str(tmp_path),
            "--wavelengths",
            "4",
            "6",
            "--alpha",
            "0.2",
            "--nugget-mode",
            "vst",
        ]
    )
    report = strict_read(tmp_path / "synthetic/example.json")
    assert report["dimensions"]["nugget"]["nugget_mode"] == "vst"
    assert report["dimensions"]["power"]["alpha"] == 0.2
    md = (tmp_path / "synthetic/example.md").read_text(encoding="utf-8")
    for label in (
        "Group counts",
        "Physical wavelength shares",
        "Acquisition design flags",
        "Top 10 pseudoreplication",
        "Top 10 nugget-sensitive",
    ):
        assert label in md


def test_quick_run_and_public_api(tmp_path):
    c = cohort()
    frame = pd.concat(
        [
            pd.DataFrame(
                {
                    "unit": s.patient_id,
                    "arm": s.group,
                    "lot": s.batch,
                    "slide_id": None if i == 0 else "s",
                },
                index=[f"{i}-{j}" for j in range(len(s.coords))],
            )
            for i, s in enumerate(c.sections)
        ]
    )
    adata = ad.AnnData(
        np.vstack([s.X for s in c.sections]),
        obs=frame,
        obsm={"spatial": np.vstack([s.coords for s in c.sections])},
    )
    path = tmp_path / "data.h5ad"
    adata.write_h5ad(path)
    main(
        [
            "quick-run",
            str(path),
            "--patient-col",
            "unit",
            "--group-col",
            "arm",
            "--batch-col",
            "lot",
            "--design",
            "stratified",
            "--out",
            str(tmp_path),
            "--wavelengths",
            "auto",
            "--cohort-id",
            "test",
        ]
    )
    report = strict_read(tmp_path / "local/test.json")
    other = audit_anndata(
        adata, "unit", "arm", batch_col="lot", design="stratified", wavelengths="auto"
    )
    assert report["dimensions"] == other["dimensions"]
    assert report["dimensions"]["design"]["categorical"]["slide_id"]["missing_fraction"] == 1 / 8


def test_bilingual_count_pairs_and_negation(tmp_path):
    path = tmp_path / "methods.txt"
    path.write_text(
        "Ninety-nine donors. Twenty one sections. Patients (n = twelve vs. 14). "
        "Harmony was not applied. Without mixed effects, we used FDR. "
        "共计十八名患者。未使用批次校正。采用秩和检验。",
        encoding="utf-8",
    )
    rows = compare_methods(path, strict_read(ROOT / "results/public/GSE220442.json"))["matrix"]
    assert [(r["unit"], r["claimed"]) for r in rows if r["field"] == "unit_count"] == [
        ("patients", 99),
        ("sections", 21),
        ("patients", 18),
    ]
    assert next(r for r in rows if r["field"] == "group_counts")["claimed"] == [12, 14]
    assert [(r["strategy"], r["negated"]) for r in rows if r["field"] == "methods_mention"] == [
        ("batch adjustment", True),
        ("patient-level analysis", True),
        ("multiple testing", False),
        ("batch adjustment", True),
        ("rank test", False),
    ]


def test_execute_quickstart_tutorial(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    notebook = strict_read(ROOT / "examples/quickstart_tutorial.ipynb")
    assert notebook["nbformat"] == 4
    example = tmp_path / "quickstart.py"
    example.write_text(
        "\n".join(
            "".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"
        ),
        encoding="utf-8",
    )
    run_path(str(example))
    assert list(tmp_path.iterdir()) == [example]


def test_stratified_mc_matches_seeded_block_draws():
    labels = np.array(["A", "A", "B", "B", "A", "B"])
    strata = np.repeat(["slide1", "slide2"], 3)
    y = np.array([4.0, 2.0, 8.0, 9.0, 1.0, 7.0])
    draws = 41
    rng = np.random.default_rng(17)
    distribution = []
    for _ in range(draws):
        selected = np.r_[
            rng.choice([0, 1, 2], 1, replace=False), rng.choice([3, 4, 5], 2, replace=False)
        ]
        distribution.append(y[selected].sum())
    observed = y[labels == "B"].sum()
    counts = min(sum(v >= observed for v in distribution), sum(v <= observed for v in distribution))
    test = patient_test(
        y[:, None], labels, design="stratified", strata=strata, max_exact=1, draws=draws
    )
    assert test["p"] == pytest.approx([min(1, 2 * (counts + 1) / (draws + 1))])


def test_multilevel_association_and_partial_acquisition_alias():
    c = cohort()
    for i, s in enumerate(c.sections):
        s.batch = f"{s.group}-{i // 4}"
    fields = design_audit(c)["categorical"]
    p = fields["batch"]["permutation_p"]
    # Under fixed margins each pure 4-by-2 table places two of four levels in group B.
    from math import comb

    exact = comb(4, 2) / comb(8, 4)
    assert p == pytest.approx(exact, abs=0.035)
    c.sections[0].batch = None
    partial = design_audit(c)["categorical"]["batch"]
    assert partial["status"] == "missing_or_partial"
    assert partial["observed_patient_conditions"] == 7
    assert 0 < partial["permutation_p"] <= 1


def test_invalid_options_and_zero_exposure():
    c = cohort()
    for alpha in (0, 1, float("nan")):
        with pytest.raises(ValueError, match="alpha"):
            audit(c, alpha=alpha)
    with pytest.raises(ValueError, match="noncollinear"):
        build_voronoi(np.array([[0, 0], [1, 1], [2, 2]]))
    with pytest.raises(ValueError, match="total ion current"):
        section_statistics(replace(c.sections[0], C=None), [4], leave_one_out_tic=True)
    for s in c.sections:
        s.library_size = np.zeros(len(s.coords))
    result = pseudorep_audit(c)
    assert result["group_depth_ratio"] is None and result["depth_normalized_patient_rate"] is None
