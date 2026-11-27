import importlib.util
from pathlib import Path

import pytest

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_receipt_append_and_failure_atomicity(tmp_path):
    script = load_script("fetch_and_verify")

    def register(name):
        return {
            "files": [
                {"accession": "public-test", "url": "https://example.org/" + name, "name": name}
            ]
        }

    for name in ["meta", "matrix", "meta"]:
        (tmp_path / name).write_text(name, encoding="utf-8")
        script.fetch(register(name), tmp_path, local=True)
    manifest = tmp_path / "dataset_manifest.json"
    old = strict_read(manifest)
    assert len(old["files"]) == 2
    invalid = register("meta")
    invalid["files"][0]["sha256"] = "unmatched"
    with pytest.raises(ValueError):
        script.fetch(invalid, tmp_path, local=True)
    assert strict_read(manifest) == old
    for f in old["files"]:
        assert f["sha256"] == sha256(tmp_path / f["name"])


def test_metadata_record_cannot_become_a_dashboard_audit(tmp_path):
    suite = load_script("run_benchmark_suite")
    res = tmp_path / "results"
    write_json(
        res / "synthetic/memory.json", {"kind": "synthetic", "dimensions": ["power", "design"]}
    )
    suite.publish(res, tmp_path / "docs")
    assert strict_read(tmp_path / "docs/audit_summary.json")["reports"] == []


def test_archive_md5_checked_when_available(tmp_path):
    import hashlib

    script = load_script("fetch_and_verify")
    (tmp_path / "matrix").write_bytes(b"fixture")
    src = {
        "name": "matrix",
        "accession": "test",
        "url": "https://example.org/matrix",
        "archive_checksum": "md5:" + hashlib.md5(b"fixture").hexdigest(),
    }
    assert script.fetch({"files": [src]}, tmp_path, local=True)[0]["reference_digest_matched"]
    with pytest.raises(ValueError, match="Archive md5 mismatch"):
        script.fetch(
            {"files": [dict(src, archive_checksum="md5:" + "0" * 32)]}, tmp_path, local=True
        )


def test_mislabeled_artifact_rejected(tmp_path):
    suite = load_script("run_benchmark_suite")
    write_json(tmp_path / "results/public/z.json", {"kind": "synthetic"})
    with pytest.raises(ValueError, match="Provenance"):
        suite.publish(tmp_path / "results", tmp_path / "docs")


def test_no_dataset_id_hallucination():
    catalog = strict_read(ROOT / "docs/data/catalog.json")["cohorts"]
    assert len({r["id"] for r in catalog if r["scope"] == "human"}) >= 20
    for c in catalog:
        assert c["source"].startswith("https://") and c["contrast"] and c["method"]


def test_doi_accession_and_portable_report_ids():
    match = load_script("run_benchmark_suite").match_accession
    row = {"id": "10.5281/zenodo.6578047", "source": "https://zenodo.org/records/6578047"}
    report = {"id": "zenodo-6578047-FZ", "source": row["source"]}
    assert match(report, [row]) == row["id"]
    assert match(dict(report, source="other", accession=row["id"]), [row]) == row["id"]
    assert match(dict(report, source="other"), [row]) is None
    with pytest.raises(ValueError, match="Multiple"):
        match(report, [row, dict(row, id="ambiguous")])


def test_gate_reconstructs_native_bh_with_an_unestimable_feature():
    verify = load_script("verify_release_gate").verify_native_fit
    family = {"n_tests": 2, "n_family": 3, "n_called": 1}
    panel = dict(family, p=[0.01, None, 0.9], q=[0.03, None, 1.0])
    row = dict(family, panels=[panel], call_rate=1 / 3, wall_seconds=1.0, peak_rss_bytes=128)
    verify(row, 3, 0.05)
    panel["q"][0] = 0.02
    with pytest.raises(ValueError, match="BH reconstruction"):
        verify(row, 3, 0.05)
