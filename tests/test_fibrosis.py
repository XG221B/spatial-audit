"""Crosswalk and bounded XDR projection contracts with synthetic byte fixtures."""

import gzip
import hashlib
import struct
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_crosswalk_selects_one_core_and_retains_acquisition_differences(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from prepare_fibrosis import select_samples

    observed = pd.DataFrame(
        {
            "sample": ["caseB", "caseA", "control"],
            "patient": ["case", "case", "ctrl"],
            "disease_status": ["Disease", "Disease", "Control"],
            "tma": ["TMA2", "TMA1", "TMA1"],
            "run": ["Run2", "Run1", "Run1"],
        }
    )
    summary = pd.DataFrame(
        {
            "patient": ["case", "ctrl"],
            "disease_status": ["Disease", "Control"],
            "tmas": ["1 & 2", "1"],
            "runs": ["1", "1"],
        }
    ).set_index("patient")
    result = select_samples(observed, summary).set_index("sample")
    assert result.loc["caseA", "retained"] and not result.loc["caseB", "retained"]
    assert not result.loc["caseB", "run_in_donor_summary"]
    observed.loc[0, "disease_status"] = "Control"
    with pytest.raises(ValueError, match="clinical IDs or status"):
        select_samples(observed, summary)


def test_streaming_release_projection_retains_original_count_order_and_depth(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    import project_fibrosis_rds as projection

    n, m = 5, 6
    counts = np.arange(n * m).reshape(n, m)
    payload = bytearray(b"fixture")
    offsets = []
    for name in ["counts", "data", "scale.data"]:
        payload.extend(struct.pack(">II", 0x40009, len(name)) + name.encode())
        offsets.append((len(payload), name))
        payload.extend(struct.pack(">II", 0x20E, n * m))
        payload.extend(counts.astype(">f8").tobytes())
        payload.extend(b"attribute pairlist")
    source = tmp_path / "source.gz"
    source.write_bytes(gzip.compress(bytes(payload)))
    monkeypatch.setattr(projection, "DIGEST", hashlib.sha256(source.read_bytes()).hexdigest())
    for k, v in {"N_TARGETS": m, "N_CELLS": n, "PANEL_SIZE": 2, "SLOTS": offsets}.items():
        monkeypatch.setattr(projection, k, v)
    projection.project(source, tmp_path / "result")
    np.testing.assert_array_equal(
        np.fromfile(tmp_path / "result/counts.i32", dtype="<i4").reshape(n, 2), counts[:, :2]
    )
    np.testing.assert_array_equal(
        np.fromfile(tmp_path / "result/depth.f64", dtype="<f8"), counts.sum(axis=1)
    )
    assert (tmp_path / "result/metadata.rds").read_bytes().count(b"attribute pairlist") == 3
