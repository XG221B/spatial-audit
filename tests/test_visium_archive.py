"""GEO reconstruction checks for barcode joins, units, normalized X and provenance."""

import hashlib
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
    return runpy.run_path(str(ROOT / "scripts/prepare_gse220442.py"))


def tiny(folder):
    X = np.array([[0, 2, 8], [5, 0, 0], [0, 7, 0]], dtype=np.int32)
    store = sparse.csr_matrix(X)
    with h5py.File(folder / "filtered_feature_bc_matrix.h5", "w") as handle:
        m = handle.create_group("matrix")
        for k in ("data", "indptr", "indices"):
            m[k] = getattr(store, k)
        m["shape"] = X.shape[::-1]
        m["features/id"] = [b"id-a", b"id-b", b"id-c"]
        m["barcodes"] = [b"ONE-1", b"TWO-1", b"THREE-1"]
    (folder / "scalefactors_json.json").write_text(
        '{"spot_diameter_fullres": 110}', encoding="utf-8"
    )
    (folder / "tissue_positions_list.csv").write_text(
        "TWO-1,1,1,1,180,260\nTHREE-1,0,0,2,0,0\nONE-1,1,0,0,100,200\n", encoding="utf-8"
    )
    m = pd.DataFrame(
        {
            "archived": ["TWO-DONOR", "ONE-DONOR"],
            "patientID": ["DONOR", "DONOR"],
            "category": ["AD", "AD"],
            "Layer": ["Layer 6", "White Matter"],
        }
    )
    panel = pd.DataFrame({"ensembl_id": ["id-c", "id-a"]}, index=["C", "A"])
    return m, panel


def test_source_alignment_and_full_library_exposure(tmp_path, adapter):
    metadata, panel = tiny(tmp_path)
    section, rec = adapter["make_section"](tmp_path, "DONOR", metadata, panel)
    assert section.obs_names.tolist() == ["ONE-1", "TWO-1"]
    np.testing.assert_array_equal(section.layers["counts"].toarray(), [[8, 0], [0, 5]])
    np.testing.assert_allclose(section.X.toarray(), np.log1p([[8000, 0], [0, 10000]]), atol=1e-6)
    np.testing.assert_allclose(section.obsm["spatial"], [[0, 0], [50, 50 * np.sqrt(3)]])
    assert np.linalg.norm(np.diff(section.obsm["spatial"],axis=0)) == pytest.approx(100)
    assert rec["legacy_image_derived_pitch_um"] == 50
    assert section.obs["ct_White_Matter"].tolist() == [1, 0]


@pytest.mark.parametrize("failure", ["label", "barcode", "feature", "diameter"])
def test_archive_rejects_ambiguous_design(tmp_path, adapter, failure):
    metadata, panel = tiny(tmp_path)
    if failure == "label":
        metadata.loc[0, "category"] = "control"
    elif failure == "barcode":
        metadata.loc[0, "archived"] = "OTHER-DONOR"
    elif failure == "feature":
        panel.loc["C", "ensembl_id"] = "absent"
    else:
        (tmp_path / "scalefactors_json.json").write_text(
            '{"spot_diameter_fullres":0}', encoding="utf-8"
        )
    with pytest.raises(ValueError):
        adapter["make_section"](tmp_path, "DONOR", metadata, panel)


def test_frozen_panel_and_array_equivalence():
    evidence = strict_read(ROOT / "results/public/GSE220442_reconstruction.json")
    panel = ROOT / "data/GSE220442/feature_panel.tsv"
    assert hashlib.sha256(panel.read_bytes()).hexdigest() == evidence["feature_panel_sha256"]
    names = pd.read_csv(panel, sep="\t")
    assert len(names) == evidence["panel_size"] and names.ensembl_id.is_unique
    assert sum(r["spots"] for r in evidence["sections"]) == evidence["n_observations"]
    assert all(all(r["reference_comparison"].values()) for r in evidence["sections"])
    assert len({r["patient_id"] for r in evidence["sections"]}) == evidence["n_patients"]


def test_full_panel_uses_reconstructed_inputs():
    evidence = strict_read(ROOT / "results/public/GSE220442_reconstruction.json")
    full = strict_read(ROOT / "results/public/GSE220442-full-panel.json")
    assert {r["sha256"] for r in full["provenance"]["inputs"]} == {
        r["prepared_sha256"] for r in evidence["sections"]
    }
    panel = pd.read_csv(ROOT / "data/GSE220442/feature_panel.tsv", sep="\t")
    assert full["dimensions"]["pseudorep"]["features"] == panel["symbol"].tolist()
    assert full["n_features"] == len(panel) and full["kind"] == "public"
