"""Biopsy assignment, count-layer and coordinate checks for the Xenium colitis archive."""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_uc_xenium as uc


def test_post_treatment_and_unassigned_cells_excluded():
    clinical = {"HS45": "UC", "HS46": "UC", "HS31": "HC"}
    assert uc.baseline_patient("HS45_POST_VDZ_NR", clinical) is None
    assert uc.baseline_patient("unassigned", clinical) is None
    assert uc.baseline_patient("HS46__PRE_VDZ_NR", clinical) == ("HS46", "UC")
    assert uc.baseline_patient("HS31_HC", clinical) == ("HS31", "HC")
    with pytest.raises(ValueError, match="disagrees"):
        uc.baseline_patient("HS46_HC", clinical)


@pytest.mark.parametrize("bad", ["none", "counts", "centroid", "prefix"])
def test_annotation_join_through_raw_barcode_across_runs(tmp_path, bad):
    matrices = [[4, 2], [3, 1], [2, 4], [1, 3]]
    counts = sparse.csr_matrix(matrices)
    rows = pd.DataFrame(
        {"x_centroid": [1.0, 2.0, 15011.0, 15012.0], "y_centroid": [0.0, 1.0, 3.0, 4.0]},
        index=["Sept-c1", "Sept-c2", "Nov-c1", "Nov-c2"],
    )
    raw = ad.AnnData(X=counts, obs=rows, var=pd.DataFrame(index=["G1", "G2"]))
    rp, ap = tmp_path / "raw.h5ad", tmp_path / "annotated.h5ad"
    raw.write_h5ad(rp)
    order = [2, 1, 3, 0]
    frame = rows.iloc[order].copy()
    frame.index = pd.Index(["Nov-c1", "c2", "Nov-c2", "c1"], name="cell_id")
    frame["batch"] = pd.Categorical(["Nov", "Sept", "Nov", "Sept"])
    frame["24_01_17_HS"] = pd.Categorical(["HS31_HC", "HS46__PRE_VDZ_NR", "unassigned", "HS31_HC"])
    frame["24_01_17_Condition"] = pd.Categorical(["HC", "PRE_VDZ_NR", "unassigned", "HC"])
    frame["Patient_ID_cores_combined"] = pd.Categorical(
        ["HS31_HC_r1", "HS46_VDZ_r2", "NA", "HS31_HC_r1"]
    )
    frame["25_06_11_Common_Coarse_Xenium_Combinedv3"] = pd.Categorical(["T", "B", "B", "T"])
    if bad == "centroid":
        frame.iloc[0, 0] += 1
    if bad == "prefix":
        frame.index = pd.Index(["Sept-c1", "c2", "Nov-c2", "c1"], name="cell_id")
    ann = ad.AnnData(X=counts[order].copy(), obs=frame, var=raw.var)
    ann.layers["raw_counts"] = ann.X.copy()
    if bad == "counts":
        ann.layers["raw_counts"].data[0] += 1
    ann.write_h5ad(ap)
    if bad != "none":
        with pytest.raises(ValueError):
            uc.project(rp, ap, {"HS31": "HC", "HS46": "UC"}, tmp_path, tmp_path)
        return
    genes, sections, details = uc.project(rp, ap, {"HS31": "HC", "HS46": "UC"}, tmp_path, tmp_path)
    assert genes == ["G1", "G2"] and len(sections) == 3
    assert details["included_cells"] == 3 and details["included_donors"] == 2
    entry = next(s for s in sections if s["batch"] == "Nov")
    with np.load(tmp_path / entry["path"], allow_pickle=False) as data:
        np.testing.assert_array_equal(data["coords"], [[11.0, 3.0]])
        np.testing.assert_array_equal(data["X"], [[2, 4]])
        np.testing.assert_array_equal(data["library_size"], [6])
