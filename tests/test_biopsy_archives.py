"""Regression checks for biopsies nested on slides and checksummed ZIP ranges."""

import io
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_colitis as colitis
import prepare_renal as renal


def test_multiple_donors_on_one_renal_slide():
    ids = np.array([1, 3])
    assert renal.section_value(np.array([0, 1, 0, 1]), np.array(["A", "B"]), ids) == "B"
    with pytest.raises(ValueError, match="Biopsy assignment"):
        renal.section_value(np.array([0, 1, 0, 0]), np.array(["A", "B"]), ids)


def test_renal_subsetting_preserves_discontiguous_rows_and_feature_order():
    rows = np.array([0, 2, 8191, 8192, 8198])
    columns = np.array([3, 0])
    values = np.arange(8200 * 5).reshape(8200, 5)
    np.testing.assert_array_equal(
        renal.select_rows(sparse.csr_matrix(values), rows, columns), values[rows][:, columns]
    )


def test_colitis_resolves_donor_and_therapy_from_archived_fields():
    sample = """^SAMPLE = GSM100
!Sample_title = GI Biobank Corona Recovery 17 (GI 6266)
!Sample_characteristics_ch1 = section cut: Biopsy
!Sample_characteristics_ch1 = activity: Inflamed
!Sample_characteristics_ch1 = disease: Ulcerative Colitis
!Sample_characteristics_ch1 = immunotherapy: N/A
!Sample_supplementary_file_1 = ftp://example.org/selected.zip
"""
    assert colitis.biopsy_design(sample)[0]["patient_id"] == "GI6266"
    assert colitis.biopsy_design(sample.replace("Biopsy", "Resection")) == []
    with pytest.raises(ValueError, match="repeated"):
        colitis.biopsy_design(sample + sample.replace("GSM100", "GSM101"))


def test_remote_archive_rejects_changed_version(monkeypatch):
    class Response(io.BytesIO):
        status = 206

        def __init__(self, payload, lo, hi, version):
            super().__init__(payload)
            self.headers = {"Content-Range": f"bytes {lo}-{hi}/100000", "Last-Modified": version}

    replies = iter([Response(b"abcd", 99996, 99999, "a"), Response(bytes(65536), 0, 65535, "b")])
    monkeypatch.setattr(colitis, "urlopen", lambda *args, **kwargs: next(replies))
    with colitis.RemoteZip("https://example.org/archive.zip") as remote:
        remote.seek(-4, 2)
        assert remote.read() == b"abcd"
        remote.seek(0)
        with pytest.raises(ValueError, match="changed"):
            remote.read(4)
