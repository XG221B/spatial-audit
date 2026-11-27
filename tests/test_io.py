import json
from pathlib import Path

import numpy as np
import pytest

from spatial_audit.claim_check import compare_methods
from spatial_audit.cli import export_fixtures
from spatial_audit.fixtures import fixtures
from spatial_audit.schema import archive_accession, clean, load_manifest, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]


def test_accession_is_read_from_the_authoritative_geo_parameter():
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442"
    assert archive_accession(None, url) == "GSE220442"
    assert archive_accession(None, url.replace("ncbi.nlm.nih.gov", "example.org")) is None
    with pytest.raises(ValueError, match="disagrees"):
        archive_accession("GSE348499", url)


@pytest.mark.parametrize("literal", ["Infinity", "-Infinity", "NaN", "1e999"])
def test_reject_nonfinite_input(tmp_path, literal):
    path = tmp_path / "bad.json"
    path.write_text('{"value":' + literal + "}", encoding="utf-8")
    with pytest.raises(ValueError):
        strict_read(path)


def test_nonfinite_output(tmp_path):
    dest = tmp_path / "a.json"
    write_json(dest, {"arr": np.array([np.nan, np.inf, -np.inf, 0.0]), "id": np.int64(4)})
    assert strict_read(dest) == {"arr": [None, None, None, 0], "id": 4}
    assert clean(np.bool_(True)) is True


def test_archive_evidence_name_and_noise_scale_round_trip(tmp_path):
    from spatial_audit.report import audit

    export_fixtures([next(iter(fixtures()))], tmp_path)
    path = tmp_path / "manifest.json"
    m = strict_read(path)
    m["cohorts"][0].update(preprocessing_evidence="study_inputs.json", measurement="log_normalized")
    write_json(path, m)
    report = audit(load_manifest(path)[0])
    assert report["preprocessing_evidence"] == "study_inputs.json"
    assert report["dimensions"]["nugget"]["measurement_scope"]["measurement"] == "log_normalized"
    m["cohorts"][0]["preprocessing_evidence"] = "../elsewhere.json"
    write_json(path, m)
    with pytest.raises(ValueError, match="portable"):
        load_manifest(path)


def test_roundtrip_export(tmp_path):
    c = next(iter(fixtures()))
    export_fixtures([c], tmp_path)
    spec = tmp_path / "manifest.json"
    again = load_manifest(spec)[0]
    assert again.features == c.features and len(again.sections) == len(c.sections)
    assert np.array_equal(again.sections[0].X, c.sections[0].X)
    data = strict_read(spec)
    data["cohorts"][0]["files"][0]["sha256"] = "invalid"
    write_json(spec, data)
    with pytest.raises(ValueError, match="Checksum"):
        load_manifest(spec)


def test_repeated_data_file_cannot_create_a_donor(tmp_path):
    export_fixtures([next(iter(fixtures()))], tmp_path)
    path = tmp_path / "manifest.json"
    data = strict_read(path)
    files = data["cohorts"][0]["files"]
    files.append(dict(files[0], patient_id="copy"))
    write_json(path, data)
    with pytest.raises(ValueError, match="Repeated data file"):
        load_manifest(path)


@pytest.mark.parametrize("change", ["units", "shared_patient", "missing_source", "unknown_format"])
def test_invalid_manifest(tmp_path, change):
    export_fixtures([next(iter(fixtures()))], tmp_path)
    path = tmp_path / "manifest.json"
    data = strict_read(path)
    spec = data["cohorts"][0]
    if change == "units":
        spec["units"] = "pixels"
    if change == "shared_patient":
        other = next(e for e in spec["files"] if e["group"] != spec["files"][0]["group"])
        other["patient_id"] = spec["files"][0]["patient_id"]
    if change == "missing_source":
        spec.update(kind="public", source="")
    if change == "unknown_format":
        spec["files"][0].update(path="matrix.xyz")
        del spec["files"][0]["sha256"]
    write_json(path, data)
    with pytest.raises(ValueError):
        load_manifest(path)


def test_claim_scope_not_verdict(tmp_path):
    p = tmp_path / "methods.txt"
    p.write_text("20 samples. 6 patients. 72 cells. Harmony was applied.", encoding="utf-8")
    report = strict_read(ROOT / "results/public/GSE220442.json")
    r = compare_methods(p, report)
    states = {v.get("unit"): v["comparison"] for v in r["matrix"]}
    assert (
        states["samples"] == "unit_unmapped"
        and states["patients"] == "numeric_agreement_check_scope"
    )
    assert states[None] == "manual_scope_review" and r["methods_sha256"]


def test_unreadable_pdf(tmp_path):
    from pypdf import PdfWriter

    pdf = tmp_path / "empty.pdf"
    w = PdfWriter()
    w.add_blank_page(width=200, height=200)
    w.write(pdf)
    with pytest.raises(ValueError, match="transcript"):
        compare_methods(pdf, {})


def test_no_nonatomic_write_on_failure(tmp_path):
    path = tmp_path / "value.json"
    write_json(path, {"old": True})
    with pytest.raises(TypeError):
        write_json(path, object())
    assert json.loads(path.read_text()) == {"old": True}
    assert list(tmp_path.iterdir()) == [path]
