"""Contract tests for sample handling, missing models and the reported denominators."""

import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmread

from spatial_audit.schema import Cohort, Section, strict_read

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_sota_benchmark.py"
spec = spec_from_file_location("sota", SCRIPT)
assert spec and spec.loader
sota = module_from_spec(spec)
spec.loader.exec_module(sota)


def cohort(platform="Visium"):
    sections = []
    for i in range(6):
        xy = np.column_stack([np.arange(20), np.full(20, i)])
        sections.append(
            Section(
                f"p{i}",
                "A" if i < 3 else "B",
                xy,
                sparse.csr_matrix(np.ones((20, 2))),
                section_id=f"s{i}",
                domain=np.repeat(["one", "two"], 10),
            )
        )
    return Cohort("synthetic-r", "synthetic", ["a", "b"], sections, platform, "")


@pytest.mark.parametrize("size,cap", [(24, None), (24, 12), (7, 12), (12, 12)])
def test_spot_sampling(size, cap):
    idx = sota.choose_rows(size, cap, 1, 0)
    assert len(idx) == min(size, cap or size)
    assert np.all(np.diff(idx) > 0)
    assert np.array_equal(idx, sota.choose_rows(size, cap, 1, 0))


@pytest.mark.parametrize("method", list(sota.ESTIMANDS))
def test_count_platform_eligibility(method):
    reason = sota.inapplicable(cohort(), method, True)
    assert (reason is None) == (method != "Cardinal")


@pytest.mark.parametrize("method", list(sota.ESTIMANDS))
def test_msi_platform_eligibility(method):
    reason = sota.inapplicable(cohort("MALDI-MSI"), method, False)
    assert (reason is None) == (method == "Cardinal")


@pytest.mark.parametrize("method", ["edgeR", "DESeq2", "DESpace2", "nnSVG", "SPARK-X"])
def test_count_models_refuse_unspecified_scale(method):
    assert sota.inapplicable(cohort(), method, False)


@pytest.mark.parametrize("bad", [-1, 0.25])
def test_malformed_counts(bad):
    c = cohort()
    c.sections[-1].X.data[0] = bad
    assert sota.inapplicable(c, "DESeq2", True)


def test_missing_or_single_domain():
    c = cohort()
    for section in c.sections:
        section.domain = np.repeat("one", len(section.coords))
    assert sota.inapplicable(c, "DESpace2", True)
    c.sections[0].domain = None
    assert sota.inapplicable(c, "DESpace2", True)


def test_patient_metadata_and_feature_order(tmp_path):
    c = cohort()
    c.sections[1].patient_id = c.sections[0].patient_id
    sizes = sota.export_input(c, tmp_path, 12, 4)
    assert sizes == [12] * 6
    assert strict_read(tmp_path / "input.json")["features"] == c.features
    assert len(pd.read_csv(tmp_path / "patients.csv")) == 5
    assert mmread(tmp_path / "s00001.mtx").shape == (2, 12)
    assert len(pd.read_csv(tmp_path / "s00002.csv")) == 12


def test_despace_requires_shared_domain_support():
    c = cohort()
    for section in c.sections:
        section.domain = np.repeat(["shared", section.group], 10)
    assert "shared across groups" in sota.inapplicable(c, "DESpace2", True)


def test_unassigned_domains_do_not_establish_a_shared_anatomical_contrast():
    c = cohort()
    for s in c.sections:
        s.domain = np.repeat(["Layer A", " Noise "], 10)
    assert "shared across groups" in sota.inapplicable(c, "DESpace2", True)
    assert (sota.domain_values(c.sections[0]) == "").sum() == 10


def test_registered_full_depth_and_invalid_depth():
    s = cohort().sections[0]
    s.library_size = np.arange(3, 23)
    np.testing.assert_array_equal(sota.depths(s, np.array([1, 4]), np.array([2, 2])), [4, 7])
    s.library_size = np.zeros(20)
    with pytest.raises(ValueError, match="at least"):
        sota.depths(s, np.array([1, 4]), np.array([2, 2]))


def test_nnsvg_accepts_explicit_deposited_logcounts():
    c = cohort()
    c.provenance["measurement"] = "log_normalized"
    assert sota.inapplicable(c, "nnSVG", False) is None
    assert sota.inapplicable(c, "DESeq2", False) is not None
    assert "directly" in sota.normalization(c, "nnSVG")


@pytest.mark.parametrize("method", list(sota.ESTIMANDS))
def test_declared_universe_and_estimand(method):
    panel = {
        "p": [0.01, None, 0.9],
        "q": [0.03, None, 1.0],
        "n_tests": 2,
        "n_family": 3,
        "n_called": 1,
    }
    row = {"method": method, "status": "ok", "panels": [panel], "spot_minus_method": None}
    sota.summarize(row, ["x", "y", "z"], 2 / 3)
    assert row["call_rate"] == pytest.approx(1 / 3)
    assert row["tested_fraction"] == pytest.approx(2 / 3)
    if method in sota.MARGINAL:
        assert row["spot_minus_method"] == pytest.approx(1 / 3)
    else:
        assert row["spot_minus_method"] is None


def test_unavailable_r_is_not_a_zero_call_rate():
    meta, rows = sota.benchmark(cohort(), True, None, ["edgeR", "Cardinal"], 1, None, 17)
    assert meta["kind"] == "synthetic"
    assert [r["status"] for r in rows] == ["unavailable", "not_applicable"]
    assert all(r["call_rate"] is None and r["wall_seconds"] is None for r in rows)


def test_timeout_reaps_worker_and_preserves_diagnostic(tmp_path, monkeypatch):
    launch = subprocess.Popen
    processes: list = []

    def slow_worker(args, **kwargs):
        process = launch([sys.executable, "-c", "import time;time.sleep(10)"], **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr(sota.subprocess, "Popen", slow_worker)
    result = sota.invoke(sys.executable, "edgeR", tmp_path, timeout=0.12)
    assert result["status"] == "timeout"
    assert result["wall_seconds"] >= 0.12
    assert result["peak_rss_bytes"] > 0
    assert processes[0].poll() is not None
    assert (tmp_path / "run_sota_benchmark.R").read_bytes() == sota.WORKER.read_bytes()


def test_malformed_worker_output():
    row = {"method": "edgeR", "status": "ok", "panels": [{"p": [], "q": [], "n_family": 2}]}
    with pytest.raises(ValueError, match="hypothesis family"):
        sota.summarize(row, ["a", "b"], 0.5)


def test_empty_tests_remain_missing():
    row = {
        "method": "DESeq2",
        "status": "ok",
        "panels": [{"p": [None], "q": [None], "n_family": 1, "n_tests": 0, "n_called": 0}],
        "spot_minus_method": None,
    }
    sota.summarize(row, ["a"], 0.5)
    assert row["status"] == "no_estimable_features"
    assert row["call_rate"] is None
