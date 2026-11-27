"""Offline provenance and imzML regression tests, including synthetic binary spectra."""

import hashlib
import json
import runpy
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
CATALOG = json.loads((ROOT / "docs/data/catalog.json").read_bytes())
MSI = runpy.run_path(str(ROOT / "scripts/prepare_msi.py"))
BIB = runpy.run_path(str(ROOT / "scripts/verify_dois.py"))
STEREO = runpy.run_path(str(ROOT / "scripts/prepare_stereo.py"))
KUPPE = runpy.run_path(str(ROOT / "scripts/prepare_kuppe.py"))
PASC = runpy.run_path(str(ROOT / "scripts/prepare_pasc.py"))
IPF = runpy.run_path(str(ROOT / "scripts/prepare_ipf.py"))
BIOSTUDIES = runpy.run_path(str(ROOT / "scripts/prepare_biostudies_lung.py"))


def test_catalog_and_primary_paper_coverage():
    rows = CATALOG["cohorts"]
    references = BIB["bibliography"](ROOT / "paper/references.bib")
    assert len({row["id"] for row in rows}) == len(rows)
    human = [row for row in rows if row["scope"] == "human"]
    assert len(human) >= 20
    assert {row["platform"] for row in human} == {
        "Visium",
        "Xenium",
        "Stereo-seq",
        "CosMx SMI",
        "MALDI-MSI",
        "MERFISH",
    }
    for row in human:
        assert row["contrast"] and row["claim"] and row["method_locator"] and row["unit_note"]
        assert references[row["citation_key"]]["doi"].casefold() == row["doi"].casefold()
    assert BIB["check_snapshot"](
        ROOT / "paper/references.bib", ROOT / "paper/generated/doi_verification.json"
    ) == len(references)


def test_bibliography_rejects_unparsed_suffix(tmp_path):
    path = tmp_path / "bibliography.bib"
    path.write_bytes((ROOT / "paper/references.bib").read_bytes() + b"@article{broken")
    with pytest.raises(ValueError):
        BIB["bibliography"](path)


def tiny_pair(tmp_path):
    # Synthetic positive-mode centroid spectrum with two peaks inside
    # one target window and one peak contributing only to the TIC.
    uid = bytes(range(16))
    mz = np.asarray([100.0, 100.001, 200.0], dtype="<f8")
    counts = np.asarray([3.0, 5.0, 9.0], dtype="<f4")
    payload = uid + mz.tobytes() + counts.tobytes()
    ibd = tmp_path / "synthetic.ibd"
    ibd.write_bytes(payload)

    def cv(code, value=""):
        return f'<cvParam accession="{code}" value="{value}"/>'

    def array(ref, offset, size):
        return (
            f'<binaryDataArray><referenceableParamGroupRef ref="{ref}"/>'
            + cv("IMS:1000102", offset)
            + cv("IMS:1000103", 3)
            + cv("IMS:1000104", size)
            + "</binaryDataArray>"
        )

    xml = (
        "<mzML><fileContent>"
        + cv("IMS:1000080", uid.hex())
        + cv("IMS:1000090", hashlib.md5(payload).hexdigest())
        + "</fileContent>"
        + '<referenceableParamGroup id="mz">'
        + cv("MS:1000514")
        + cv("MS:1000523")
        + cv("MS:1000576")
        + cv("IMS:1000101")
        + "</referenceableParamGroup>"
        + '<referenceableParamGroup id="i">'
        + cv("MS:1000515")
        + cv("MS:1000521")
        + cv("MS:1000576")
        + cv("IMS:1000101")
        + "</referenceableParamGroup>"
        + "<scanSettings>"
        + cv("IMS:1000046", 50)
        + cv("IMS:1000047", 50)
        + "</scanSettings>"
        + '<spectrumList count="1"><spectrum>'
        + cv("MS:1000130")
        + "<scanList><scan>"
        + cv("IMS:1000050", 2)
        + cv("IMS:1000051", 3)
        + "</scan></scanList><binaryDataArrayList>"
        + array("mz", 16, 24)
        + array("i", 40, 12)
        + "</binaryDataArrayList></spectrum></spectrumList></mzML>"
    )
    path = tmp_path / "synthetic.imzML"
    path.write_text(xml, encoding="utf-8")
    return path, ibd


def test_imaging_window_and_tic(tmp_path):
    xml, ibd = tiny_pair(tmp_path)
    coords, x, c = MSI["extract"](xml, ibd, np.array([100.0, 300.0]), (50.0, 50.0))
    np.testing.assert_array_equal(coords, [[50.0, 100.0]])
    np.testing.assert_array_equal(x, [[8.0, 0.0]])
    np.testing.assert_array_equal(c, [[17.0, 0.0, 0.0]])


@pytest.mark.parametrize("corruption", ["uuid", "bytes", "pitch", "length"])
def test_msi_checks_input_integrity(tmp_path, corruption):
    xml, ibd = tiny_pair(tmp_path)
    if corruption in {"uuid", "bytes"}:
        data = bytearray(ibd.read_bytes())
        data[0 if corruption == "uuid" else -1] ^= 1
        ibd.write_bytes(data)
    elif corruption == "length":
        xml.write_text(xml.read_text().replace('value="40"', 'value="1000000"'))
    spacing = (100.0, 50.0) if corruption == "pitch" else (50.0, 50.0)
    with pytest.raises(ValueError):
        MSI["extract"](xml, ibd, np.array([100.0]), spacing)


def test_mass_panel_rejects_duplicate_information(tmp_path):
    maf = tmp_path / "test.tsv"
    maf.write_text("mass_to_charge\n" + "100\n" * 40)
    with pytest.raises(ValueError, match="insufficient"):
        MSI["mass_panel"](maf)


def stereo_inputs():
    # A genuine bin gap is preserved when coordinates are converted to micrometers.
    metadata = pd.DataFrame(
        {"x": [1, 4], "y": [1, 3], "nCount_Spatial": [7, 2], "nFeature_Spatial": [2, 1]}
    )
    records = pd.DataFrame(
        {
            "geneID": ["ACTB", "other", "CD3D"],
            "x": [1, 1, 4],
            "y": [1, 1, 3],
            "MIDCounts": np.log1p(10000 * np.array([2 / 7, 5 / 7, 1.0])),
        }
    )
    return metadata, records


def test_stereo_checks_every_exported_entry(tmp_path):
    metadata, records = stereo_inputs()
    gem = tmp_path / "fixture.gem.gz"
    records.to_csv(gem, sep="\t", index=False)
    xy, archived, counts, check = STEREO["get_panel"](metadata, gem, chunk_size=1)
    np.testing.assert_array_equal(xy, [[50, 50], [200, 150]])
    np.testing.assert_array_equal(counts.sum(axis=1), [2, 2])
    assert archived[0, 0] == np.float32(records.iloc[0]["MIDCounts"])
    assert check["sparse_entries_checked"] == 3
    assert check["validated_integer_reconstruction"]


@pytest.mark.parametrize("change", ["depth", "fractional", "nfeature"])
def test_stereo_refuses_unvalidated_counts(tmp_path, change):
    metadata, records = stereo_inputs()
    if change == "depth":
        records = records.iloc[[0, 2]]
    elif change == "fractional":
        records.loc[1, "MIDCounts"] += 0.1
    else:
        metadata.loc[0, "nFeature_Spatial"] += 1
    gem = tmp_path / "fixture.gem.gz"
    records.to_csv(gem, sep="\t", index=False)
    _, _, _, check = STEREO["get_panel"](metadata, gem)
    assert not check["validated_integer_reconstruction"]


@pytest.mark.parametrize("change", ["duplicate", "unmatched"])
def test_stereo_rejects_invalid_spatial_join(tmp_path, change):
    metadata, records = stereo_inputs()
    if change == "duplicate":
        records = pd.concat([records, records.iloc[[0]]], ignore_index=True)
    else:
        records.loc[0, "y"] += 10
    gem = tmp_path / "fixture.gem.gz"
    records.to_csv(gem, sep="\t", index=False)
    with pytest.raises(ValueError):
        STEREO["get_panel"](metadata, gem, chunk_size=2)


def test_backed_gene_extraction_preserves_requested_order():
    data = np.arange(15).reshape(3, 5)
    indices = [3, 0, 4, 1]
    np.testing.assert_array_equal(KUPPE["ordered_columns"](data, indices), data[:, indices])


@pytest.mark.parametrize("axes", [[0, 1], [1, 0]])
@pytest.mark.parametrize("flip", [[False, False], [False, True], [True, False], [True, True]])
def test_vendor_lattice_verifies_rotated_displays(axes, flip):
    idx = pd.Index(["a", "b", "c", "d"])
    grid = pd.DataFrame([[0, 0], [1, 1], [0, 2], [2, 2]], index=idx)
    rc = grid.to_numpy() / 2
    display = np.where(flip, 1 - rc[:, axes], rc[:, axes])
    xy, transform = KUPPE["physical"](grid, idx, display)
    assert transform == {"array_axes": axes, "reflected_display_axes": flip}
    np.testing.assert_allclose(np.linalg.norm(xy[0] - xy[1]), 100)
    np.testing.assert_allclose(np.linalg.norm(xy[0] - xy[2]), 100)
    display[1, 0] += 0.05
    with pytest.raises(ValueError):
        KUPPE["physical"](grid, idx, display)


@pytest.mark.parametrize(
    "title,disease,valid",
    [
        ("Control1", "control", True),
        ("COVID1", "COVID", True),
        ("Control1", "COVID", False),
        ("COVID1", "control", False),
    ],
)
def test_pasc_sample_group_crosscheck(tmp_path, title, disease, valid):
    src = tmp_path / "sample.soft"
    src.write_text(
        f"!Sample_title = {title}\n!Sample_characteristics_ch1 = disease state: {disease}\n"
        "!Sample_geo_accession = GSM1\n!Sample_organism_ch1 = Homo sapiens\n"
    )
    if valid:
        assert PASC["read_sample"](src, "GSM1")[0] == title
    else:
        with pytest.raises(ValueError):
            PASC["read_sample"](src, "GSM1")


def test_ipf_independent_pixel_pitch_calibration():
    rc = np.array([[0, 0], [0, 2], [1, 1], [1, 3]])
    xy = np.column_stack((50 * rc[:, 1], 50 * np.sqrt(3) * rc[:, 0]))
    measured, proof = IPF["checked_visium_lattice"](rc, xy * 2 + 7, 130)
    np.testing.assert_allclose(measured, xy)
    np.testing.assert_allclose(proof["image_pitch_p05_p50_p95_um"], 100)
    with pytest.raises(ValueError, match="pitch"):
        IPF["checked_visium_lattice"](rc, xy, 200)
    # A registered image can stretch independently of the manufactured grid.
    measured, proof = IPF["checked_visium_lattice"](rc, xy, 200, rc.copy())
    np.testing.assert_allclose(measured, xy)
    assert not proof["image_pitch_within_tolerance"]
    assert proof["vendor_addresses_checked"] == len(rc)
    with pytest.raises(ValueError, match="manufacturer"):
        IPF["checked_visium_lattice"](rc, xy, 200, rc[::-1])


def test_ipf_donor_design_preserves_adjacent_sections():
    rows = []
    for group, n in [("Control", 4), ("IPF", 3)]:
        for i in range(1, n + 1):
            for tech in ["Visium", "CytAssist"] if i < 3 else ["Visium"]:
                donor = f"Donor_{group}_{i}"
                rows.append(
                    {
                        "sampleID": donor + tech,
                        "patient": donor,
                        "technology": tech,
                        "treatment": "Healthy" if group == "Control" else "IPF",
                        "VisiumSlides": "slide",
                        "array_col": 0,
                        "array_row": 0,
                        "Niche_NMF": "a",
                    }
                )
    frame = pd.DataFrame(rows)
    assert len(IPF["check_design"](frame)) == 11
    with pytest.raises(ValueError):
        IPF["check_design"](frame.iloc[1:])


@pytest.mark.parametrize("bad_value", [False, True])
def test_ipf_block_count_projection_rejects_noninteger_layer(tmp_path, bad_value):
    x = np.array([[1.0, 3.0, 4.0, 0.0], [0.0, 2.0, 5.0, 1.0], [6.0, 0.0, 0.0, 1.0]])
    if bad_value:
        x[1, 1] += 0.75
    m = sparse.csr_matrix(x)
    with h5py.File(tmp_path / "counts.h5", "w") as f:
        f.attrs["encoding-type"] = "csr_matrix"
        f.attrs["shape"] = m.shape
        for k in ["data", "indptr", "indices"]:
            f.create_dataset(k, data=getattr(m, k))
        if bad_value:
            with pytest.raises(ValueError, match="Noninteger"):
                IPF["counts_and_depth"](f, [2, 0], block=1)
        else:
            extracted, depth = IPF["counts_and_depth"](f, [2, 0], block=1)
            np.testing.assert_array_equal(extracted, x[:, [2, 0]])
            np.testing.assert_array_equal(depth, x.sum(axis=1))


def test_biostudies_sample_aliases_and_slide_keys(tmp_path):
    meta = pd.DataFrame(
        {
            "sample_id": ["s-A1", "s-B1", "t-C1"],
            "slide_id": ["s", "s", "t"],
            "slide_ca": ["A1", "B1", "C1"],
            "subject_alias": ["d1", "d2", "d2"],
            "condition": ["control", "IPF", "IPF"],
            "replicate": ["1", "1", "2"],
            "tissue_alias": ["b1", "b2", "b2"],
            "species": ["human"] * 3,
            "tissue": ["lung"] * 3,
        }
    )
    p = tmp_path / "cohort.tsv"
    meta.to_csv(p, sep="\t", index=False)
    assert BIOSTUDIES["samples"](p)["subject_alias"].nunique() == 2
    meta.loc[2, "condition"] = "control"
    meta.to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError):
        BIOSTUDIES["samples"](p)


def test_public_lung_studies_share_the_target_candidates():
    assert PASC["GENES"] == IPF["GENES"] == BIOSTUDIES["GENES"]
    assert len(PASC["GENES"]) == len(set(PASC["GENES"]))
