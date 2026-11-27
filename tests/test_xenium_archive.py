"""Clinical-outcome, deposited file and barcode linkage tests for the Xenium adapter."""

import gzip
import runpy
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from spatial_audit.schema import strict_read

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def adapter():
    return runpy.run_path(str(ROOT / "scripts/prepare_gse348499.py"))


def record():
    return """^SAMPLE = GSM10073485
!Sample_title = P10
!Sample_characteristics_ch1 = outcome: pathologic complete response (pCR)
!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10073nnn/GSM10073485/suppl/GSM10073485_P10_cell_feature_matrix.h5
!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10073nnn/GSM10073485/suppl/GSM10073485_P10_cells.parquet.gz
"""


def test_geo_parser_uses_outcome_not_title(adapter):
    (r,) = adapter["parse_geo"](
        record().replace("pathologic complete response (pCR)", "residual disease (RD)")
    )
    assert r["patient_id"] == "P10" and r["group"] == "RD"
    assert all(u.startswith("https://ftp.ncbi.nlm.nih.gov/") for u in r["urls"].values())


@pytest.mark.parametrize("mode", ["unknown", "missing_file", "wrong_domain", "two_outcomes"])
def test_geo_fails_on_undetermined_metadata(adapter, mode):
    s = record()
    if mode == "unknown":
        s = s.replace("pathologic complete response (pCR)", "unknown")
    elif mode == "missing_file":
        s = s.replace("_cells.parquet.gz", "_images.tif")
    elif mode == "wrong_domain":
        s = s.replace("ftp.ncbi.nlm.nih.gov", "archive.invalid")
    else:
        s = s.replace(
            "!Sample_title",
            "!Sample_characteristics_ch1 = outcome: residual disease (RD)\n!Sample_title",
        )
    with pytest.raises(ValueError):
        adapter["parse_geo"](s)


def test_gene_counts_join_reordered_centroids(tmp_path, adapter):
    pq = pytest.importorskip("pyarrow.parquet")
    pa = pytest.importorskip("pyarrow")
    with h5py.File(tmp_path / "cells.h5", "w") as h:
        m = h.create_group("matrix")
        X = sparse.csr_matrix([[1, 0, 7], [9, 3, 4]])
        for k in ["data", "indices", "indptr"]:
            m[k] = getattr(X, k)
        m["shape"] = X.shape[::-1]
        m["barcodes"] = [b"first-1", b"next-1"]
        m["features/name"] = [b"A", b"Control", b"B"]
        m["features/feature_type"] = [b"Gene Expression", b"Control Probe", b"Gene Expression"]
    with gzip.open(tmp_path / "cells.parquet.gz", "wb") as fh:
        pq.write_table(
            pa.Table.from_pandas(
                pd.DataFrame(
                    {
                        "cell_id": ["next-1", "first-1"],
                        "x_centroid": [15, 22],
                        "y_centroid": [23, 6],
                        "transcript_counts": [13, 8],
                    }
                )
            ),
            fh,
        )
    row = {"gsm": "GSM10073485", "patient_id": "P10", "group": "pCR"}
    out = adapter["section"](tmp_path / "cells.h5", tmp_path / "cells.parquet.gz", row, ["B", "A"])
    np.testing.assert_array_equal(out.X.toarray(), [[7, 1], [4, 9]])
    np.testing.assert_array_equal(out.obsm["spatial"], [[22, 6], [15, 23]])
    np.testing.assert_array_equal(out.obs["library_size"], [8, 13])
    with pytest.raises(ValueError, match="gene symbol"):
        adapter["section"](tmp_path / "cells.h5", tmp_path / "cells.parquet.gz", row, ["Control"])


def test_reference_crosschecks_for_all_donors(adapter):
    sample = strict_read(ROOT / "data/GSE348499/geo_samples.json")["patients"]
    evidence = strict_read(ROOT / "results/public/GSE348499_reconstruction.json")
    keep = {s["gsm"] for s in adapter["select_subset"](sample)}
    assert {s["gsm"] for s in evidence["patients"] if s["subset_member"]} == keep
    assert len(evidence["patients"]) == len(sample) == 44
    assert all(all(r["reference_comparison"].values()) for r in evidence["patients"])
    assert sum(r["cells"] for r in evidence["patients"]) == evidence["n_observations"]
