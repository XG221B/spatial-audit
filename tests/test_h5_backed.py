"""Backed X and non-X input fidelity across h5py dense, CSR and CSC encodings."""

from __future__ import annotations

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from spatial_audit.h5_backed import H5Matrix, classic_visium_coordinates
from spatial_audit.schema import load_manifest, strict_read, value_blocks, write_json


def file_spec(root, form, layer, mismatch=False):
    X = np.arange(48 * 6, dtype=np.float64).reshape(48, 6)
    X[np.indices(X.shape).sum(0) % 3 == 0] = 0
    raw = {"dense": lambda x: x, "csr": sparse.csr_matrix, "csc": sparse.csc_matrix}[form](X)
    ids = np.arange(48) % 4
    obs = pd.DataFrame(
        {
            "donor": ids.astype(str),
            "arm": (ids // 2).astype(str),
            "lab": np.array([0, 1, 0, 1])[ids].astype(str),
        },
        index=np.arange(len(X)).astype(str),
    )
    if mismatch:
        obs.loc["4", "lab"] = "X"
    obj = ad.AnnData(X=raw, obs=obs, var=pd.DataFrame(index=[f"f{i}" for i in range(6)]))
    if layer:
        obj.layers["counts"] = raw.copy()
        obj.X = X * 0
    obj.obsm["spatial"] = np.column_stack((ids, np.arange(48)))
    path = root / "cohort.h5ad"
    obj.write_h5ad(path)
    # Invalid unrelated layer would fail an eager materialization of all layers.
    with h5py.File(path, "r+") as h:
        junk = h["layers"].create_group("unrelated")
        junk.attrs["encoding-type"] = "unknown-unrelated-layer"
        junk.attrs["encoding-version"] = "1.0"
    spec = {
        "id": "synthetic-h5",
        "units": "um",
        "kind": "synthetic",
        "features": ["f5", "f3", "f0"],
        "columns": {"patient_id": "donor", "group": "arm", "batch": "lab"},
        "files": [{"path": path.name}],
    }
    if layer:
        spec["layer"] = "counts"
    manifest = root / "manifest.json"
    write_json(manifest, {"cohorts": [spec]})
    return X, manifest


def test_explicit_primary_matrix_key(tmp_path):
    _, path = file_spec(tmp_path, "csr", False)
    default = load_manifest(path)[0].sections[0].X.mean(axis=0)
    spec = strict_read(path)
    spec["cohorts"][0]["layer"] = "X"
    write_json(path, spec)
    np.testing.assert_array_equal(default, load_manifest(path)[0].sections[0].X.mean(axis=0))


@pytest.mark.parametrize(
    "variant", ["good", "duplicate", "parity", "fraction", "outside", "nonfinite"]
)
def test_visium_nominal_address_projection(variant):
    rc = np.array([[0, 0], [0, 2], [1, 1], [1, 3]], float)
    if variant != "good":
        rc[0] = {
            "duplicate": [0, 2],
            "parity": [0, 1],
            "fraction": [1.5, 3],
            "outside": [78, 2],
            "nonfinite": [0, np.nan],
        }[variant]
        with pytest.raises(ValueError, match="address"):
            classic_visium_coordinates(None, rc)
        return
    xy = classic_visium_coordinates(np.zeros_like(rc), rc)
    np.testing.assert_allclose(np.linalg.norm(xy[1:3] - xy[0], axis=1), [100, 100])
    angle = 0.63
    rotated = xy @ np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]) + 5
    np.testing.assert_array_equal(classic_visium_coordinates(rotated), rotated)
    with pytest.raises(ValueError, match="pitch"):
        classic_visium_coordinates(1.8 * rotated)
    with pytest.raises(ValueError, match="hexagonal"):
        classic_visium_coordinates(np.array([[0, 0], [0, 100], [100, 0], [100, 100]]))


@pytest.mark.parametrize("form", ["dense", "csr", "csc"])
@pytest.mark.parametrize("layer", [True, False])
def test_backed_selected_blocks(tmp_path, form, layer):
    x, manifest = file_spec(tmp_path, form, layer)
    (cohort,) = load_manifest(manifest)
    for sec in cohort.sections:
        assert isinstance(sec.X, H5Matrix)
        rows = np.flatnonzero(np.arange(48) % 4 == int(sec.patient_id))
        expect = x[rows][:, [5, 3, 0]]
        np.testing.assert_allclose(sec.X.mean(), expect.mean(0, keepdims=True))
        assert all(np.isfinite(a).all() for a in value_blocks(sec.X))
        idx = [8, 2, 5, 2]
        got = sec.X[idx, [2, 0, 2]]
        np.testing.assert_array_equal(
            got.toarray() if sparse.issparse(got) else got, expect[idx][:, [2, 0, 2]]
        )


def test_changed_backing_file_detected(tmp_path):
    _, manifest = file_spec(tmp_path, "csr", True)
    (ct,) = load_manifest(manifest)
    with h5py.File(tmp_path / "cohort.h5ad", "r+") as h:
        h.create_dataset("new-data", data=np.ones(2000))
    with pytest.raises(ValueError, match="changed"):
        ct.sections[0].X[:, :1]


def test_inconsistent_batch_annotation(tmp_path):
    _, manifest = file_spec(tmp_path, "dense", False, True)
    with pytest.raises(ValueError, match="constant within"):
        load_manifest(manifest)
