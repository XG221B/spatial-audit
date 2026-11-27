"""Calibration and archival identity checks for the GSE327879 adapter."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_cosmx.py"
spec = spec_from_file_location("cosmx_adapter", SCRIPT)
assert spec and spec.loader
adapter = module_from_spec(spec)
spec.loader.exec_module(adapter)


def calibration_tables():
    field = np.repeat([4, 7, 9, 13], 2)
    xy = np.tile([[81.0, 21], [95, 43]], (4, 1))
    origin = np.array([[20, 190], [50, 310], [110, 250], [160, 100]])
    repeated = np.repeat(origin, 2, axis=0)
    meta = pd.DataFrame(
        {
            "fov": field,
            "CenterX_local_px": xy[:, 0],
            "CenterY_local_px": xy[:, 1],
            "CenterX_global_px": repeated[:, 0] + xy[:, 0],
            "CenterY_global_px": repeated[:, 1] - xy[:, 1],
        }
    )
    positions = pd.DataFrame(
        {"FOV": field[::2], "X_mm": origin[:, 0] * 0.0002 + 4, "Y_mm": origin[:, 1] * 0.00035 - 6}
    )
    return meta, positions


def test_archival_calibration_with_offsets_and_anisotropy():
    m, p = calibration_tables()
    pitch, error = adapter.physical_pitch(m, p.iloc[::-1])
    assert pitch == pytest.approx([0.2, 0.35])
    assert error < 1e-12


@pytest.mark.parametrize(
    "defect", ["within-field", "between-field", "wrong-key", "nonfinite", "constant"]
)
def test_reject_affine_and_key_failures(defect):
    m, p = calibration_tables()
    if defect == "within-field":
        m.loc[1, "CenterY_local_px"] += 1
    elif defect == "between-field":
        p.loc[2, "X_mm"] += 0.01
    elif defect == "wrong-key":
        p.loc[0, "FOV"] = 99
    elif defect == "nonfinite":
        m.loc[2, "CenterY_local_px"] = np.nan
    else:
        m["CenterX_global_px"] = m["CenterX_local_px"] + 100
    with pytest.raises(ValueError):
        adapter.physical_pitch(m, p)


@pytest.mark.parametrize("first,last", [(2, 2), (3, 1)])
def test_match_rejects_missing_or_duplicate_cell(first, last):
    m = pd.DataFrame({"fov": [1, 1], "cell_ID": [1, 2], "patient_id": ["a", "a"]})
    x = pd.DataFrame({"fov": [1, 1], "cell_ID": [first, last], "gene": [6, 2]})
    with pytest.raises((ValueError, pd.errors.MergeError)):
        adapter.match_cells(m, x)


def test_shared_cell_numbers_on_distinct_fields_join_independently():
    m = pd.DataFrame({"fov": [11, 12], "cell_ID": [1, 1], "patient_id": ["a", "b"]})
    x = pd.DataFrame({"fov": [12, 11], "cell_ID": [1, 1], "gene": [6, 2]})
    table = adapter.match_cells(m, x).set_index("patient_id")
    assert table.loc["b", "gene"] == 6
    assert table.loc["a", "gene"] == 2
