"""Test the archived array geometry, deposited identities and count feature matching."""

import gzip
import io
import tarfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_visium.py"
spec = spec_from_file_location("visium_adapter", SCRIPT)
assert spec and spec.loader
adapter = module_from_spec(spec)
spec.loader.exec_module(adapter)


def spot_map():
    rows = np.repeat(np.arange(78), 64)
    cols = np.tile(np.arange(0, 128, 2), 78) + rows % 2
    return pd.DataFrame(
        {
            "barcode": [f"bc{i}" for i in range(rows.size)],
            "in_tissue": 1,
            "row": rows,
            "col": cols,
            "image_y": 0,
            "image_x": 0,
        }
    )


@pytest.mark.parametrize("header", [True, False])
@pytest.mark.parametrize("compressed", [True, False])
def test_staggered_spot_pitch(header, compressed):
    raw = spot_map().to_csv(index=False, header=header).encode()
    if compressed:
        raw = gzip.compress(raw)
    xy = adapter.positions_um(raw, np.array(["bc64", "bc0", "bc1"]))
    assert np.linalg.norm(xy[0] - xy[1]) == pytest.approx(100)
    assert np.linalg.norm(xy[2] - xy[1]) == pytest.approx(100)
    assert np.linalg.norm(xy[2] - xy[0]) == pytest.approx(100)


@pytest.mark.parametrize(
    "error", ["duplicate", "parity", "outside", "absent", "empty", "matrix_duplicate", "overlap"]
)
def test_reject_unverified_coordinates(error):
    frame = spot_map()
    barcodes = ["bc0", "bc1"]
    if error == "duplicate":
        frame.loc[1, "barcode"] = "bc0"
    elif error == "parity":
        frame.loc[0, "col"] = 1
    elif error == "outside":
        frame.loc[0, "in_tissue"] = 0
    elif error == "absent":
        barcodes[0] = "unmatched"
    elif error == "empty":
        frame = frame.iloc[:10]
    elif error == "matrix_duplicate":
        barcodes = ["bc0", "bc0"]
    elif error == "overlap":
        frame.loc[1, "col"] = 0
    with pytest.raises(ValueError):
        adapter.positions_um(frame.to_csv(index=False).encode(), np.array(barcodes))


def test_feature_ambiguity_and_non_rna_exclusion():
    buf = io.BytesIO()
    with h5py.File(buf, "w") as h:
        g = h.create_group("matrix")
        g["features/name"] = [b"A", b"A", b"B"]
        g["features/feature_type"] = [b"Gene Expression", b"Antibody Capture", b"Gene Expression"]
        g["barcodes"] = [b"bc0", b"bc1"]
        g["shape"] = [3, 2]
        g["data"] = [4, 8, 3, 1]
        g["indices"] = [0, 1, 2, 2]
        g["indptr"] = [0, 3, 4]
    positions = spot_map().to_csv(index=False).encode()
    x, _ = adapter.read_panel(buf, positions, ["B", "A"])
    np.testing.assert_array_equal(x, [[3, 4], [1, 0]])
    with h5py.File(buf, "r+") as h:
        h["matrix/features/feature_type"][1] = b"Gene Expression"
    with pytest.raises(ValueError, match="exactly one"):
        adapter.read_panel(buf, positions, ["A"])


def test_donor_repeats_and_selection():
    def record(i, donor, disease, status):
        return "\n".join(
            [
                f"^SAMPLE = GSM{i}",
                "!Sample_characteristics_ch1 = sequencing run / batch: batch",
                f"!Sample_characteristics_ch1 = disease state: {disease}",
                f"!Sample_characteristics_ch1 = disease status: {status}",
                f"!Sample_characteristics_ch1 = individual: {donor}",
                f"!Sample_supplementary_file = https://example.org/GSM{i}.h5",
            ]
        )

    rows = adapter.geo_records(
        "\n".join(
            [
                record(2, "P1", "Psoriasis (PsO)", "Lesional"),
                record(1, "P1", "Psoriasis (PsO)", "Lesional"),
                record(3, "P1", "PsA", "Lesional"),
                record(4, "H1", "control", "Healthy"),
                record(5, "P2", "Psoriasis (PsO)", "Non-lesional"),
            ]
        )
    )
    assert [row["section_id"] for row in rows if row["selected"]] == ["GSM1", "GSM4"]


@pytest.fixture
def covid(monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    import prepare_covid_visium

    return prepare_covid_visium


def write_reference(path, table):
    raw = table.to_csv(index=False, header=False).encode()
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("spatial/tissue_positions_list.csv")
        member.size = len(raw)
        archive.addfile(member, io.BytesIO(raw))


def test_reference_tissue_masks_are_independent_of_lattice(covid, tmp_path):
    table = spot_map()
    p, q = tmp_path / "reference1.gz", tmp_path / "reference2.gz"
    write_reference(p, table)
    table.loc[0, "in_tissue"] = 0
    write_reference(q, table.iloc[::-1])
    coords = adapter.positions_um(covid.shared_reference([p, q]), np.array(["bc0", "bc1"]))
    np.testing.assert_allclose(coords, [[0, 0], [100, 0]])
    table.loc[[0, 1], "barcode"] = ["bc1", "bc0"]
    write_reference(q, table)
    with pytest.raises(ValueError, match="disagree"):
        covid.shared_reference([p, q])


def test_vendor_reference_rejects_incompatible_chemistry(tmp_path):
    p = tmp_path / "wrong.h5"
    with h5py.File(p, "w") as h:
        h.attrs["chemistry_description"] = "Spatial 3' v4"
    with pytest.raises(ValueError, match="chemistry"):
        adapter.read_panel(p, b"", ["A"], chemistry="Spatial 3' v1")
