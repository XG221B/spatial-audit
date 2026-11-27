"""Read the selected H5AD layer and observation columns without opening unrelated arrays."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse


def _signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def _contiguous_chunks(major, budget):
    """Group sorted primary indices by bounded physical range in the source."""
    buckets = major // budget
    return np.split(major, np.flatnonzero(np.diff(buckets)) + 1)


def _compressed_range(node, major, minor, form, full_shape):
    """Decode one bounded CSR-row or CSC-column interval via its HDF5 indptr."""
    start, stop = int(major[0]), int(major[-1]) + 1
    ptr = np.asarray(node["indptr"][start : stop + 1], dtype=np.int64)
    head, tail = int(ptr[0]), int(ptr[-1])
    values = node["data"][head:tail]
    inner = node["indices"][head:tail]
    if form == "csr_matrix":
        raw = sparse.csr_matrix((values, inner, ptr - head), shape=(stop - start, full_shape[1]))
        return raw[major - start][:, minor]
    raw = sparse.csc_matrix((values, inner, ptr - head), shape=(full_shape[0], stop - start))
    return raw[:, major - start][minor, :]


@dataclass(frozen=True)
class H5Matrix:
    """A patient section view. Only an explicitly requested block is materialized."""

    path: Path
    key: str
    rows: np.ndarray
    cols: np.ndarray
    signature: tuple[int, int]
    full_shape: tuple[int, int]
    form: str

    @property
    def shape(self):
        return len(self.rows), len(self.cols)

    def __getitem__(self, key):
        import h5py

        row_key, col_key = key
        wanted_rows = np.atleast_1d(self.rows[row_key])
        wanted_cols = np.atleast_1d(self.cols[col_key])
        rows, ir = np.unique(wanted_rows, return_inverse=True)
        cols, ic = np.unique(wanted_cols, return_inverse=True)
        if _signature(self.path) != self.signature:
            raise ValueError(f"H5AD changed during audit: {self.path.name}")
        if not len(rows) or not len(cols):
            return sparse.csr_matrix((len(wanted_rows), len(wanted_cols)))
        with h5py.File(self.path, "r") as handle:
            node = handle[self.key]
            if self.form == "array":
                width = max(1, 8_388_608 // (len(cols) * node.dtype.itemsize))
                blocks = []
                col_sel: slice | np.ndarray
                if len(cols) == int(cols[-1] - cols[0] + 1):
                    col_sel = slice(int(cols[0]), int(cols[-1]) + 1)
                else:
                    col_sel = cols
                for span in _contiguous_chunks(rows, min(4096, width)):
                    start, end = int(span[0]), int(span[-1]) + 1
                    blocks.append(node[start:end, col_sel][span - start])
                return np.concatenate(blocks)[ir][:, ic]
            # Bound dense-equivalent work per decoded major interval. CSR only
            # reads selected row spans; CSC only reads selected column spans.
            transpose = self.form == "csc_matrix"
            major, minor = (cols, rows) if transpose else (rows, cols)
            width = self.full_shape[0] if transpose else self.full_shape[1]
            budget = max(1, min(4096, 1_048_576 // max(1, width)))
            pieces = [
                _compressed_range(node, span, minor, self.form, self.full_shape)
                for span in _contiguous_chunks(major, budget)
            ]
            selected = (
                sparse.hstack(pieces, format="csr")
                if transpose
                else sparse.vstack(pieces, format="csr")
            )
            return selected[ir][:, ic]

    def iter_values(self):
        if self.form == "csc_matrix":
            chunk = max(1, min(64, 1_048_576 // self.shape[0]))
            for col in range(0, self.shape[1], chunk):
                yield self[:, col : col + chunk].data
        else:
            chunk = max(1, min(4096, 1_048_576 // self.shape[1]))
            for row in range(0, self.shape[0], chunk):
                values = self[row : row + chunk, :]
                yield values.data if sparse.issparse(values) else values

    def mean(self, axis=0):
        if axis != 0:
            raise ValueError("backed feature summaries require axis=0")
        sums = np.zeros(self.shape[1])
        if self.form == "csc_matrix":
            chunk = max(1, min(64, 1_048_576 // self.shape[0]))
            for col in range(0, self.shape[1], chunk):
                part = self[:, col : col + chunk]
                sums[col : col + chunk] = np.asarray(part.sum(axis=0, dtype=float)).ravel()
        else:
            chunk = max(1, min(4096, 1_048_576 // self.shape[1]))
            for row in range(0, self.shape[0], chunk):
                part = self[row : row + chunk, :]
                sums += np.asarray(part.sum(axis=0, dtype=float)).ravel()
        return sums[None, :] / self.shape[0]


def value_blocks(x):
    if isinstance(x, H5Matrix):
        yield from x.iter_values()
    else:
        yield x.data if sparse.issparse(x) else x


def classic_visium_coordinates(xy, addresses=None):
    """Calibrate one GSE220442 section on the classic 100-um 6.5-mm Visium lattice."""
    from scipy.spatial import cKDTree

    if addresses is not None:
        rc = np.asarray(addresses, float)
        if (
            rc.ndim != 2
            or rc.shape[1] != 2
            or not np.isfinite(rc).all()
            or np.any(rc != np.floor(rc))
            or np.any(rc < 0)
            or np.any(rc[:, 0] >= 78)
            or np.any(rc[:, 1] >= 128)
            or np.any(rc.sum(axis=1) % 2)
            or len(np.unique(rc, axis=0)) != len(rc)
        ):
            raise ValueError("Invalid or duplicate classic Visium row/column addresses")
        return np.column_stack((50 * rc[:, 1], (50 * np.sqrt(3.0)) * rc[:, 0]))
    coords = np.asarray(xy, float)
    if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) < 3 or not np.isfinite(coords).all():
        raise ValueError("Need validated spatial coordinates or classic Visium array addresses")
    # The nearest-neighbour pitch alone cannot verify the shape of a 2-D
    # hexagonal lattice. Recover a basis from one physical nearest edge and
    # require every centre to map back to an integer address within 0.5 um.
    distance, neighbor = cKDTree(coords).query(coords, k=2)
    nearest = distance[:, 1]
    if np.any(nearest < 99.5) or abs(float(np.median(nearest)) - 100) > 0.5:
        raise ValueError("spatial coordinates do not have the expected 100-um Visium pitch")
    index = int(np.argmin(abs(nearest - 100)))
    edge = (coords[neighbor[index, 1]] - coords[index]) / nearest[index]
    along = edge * 100
    across = 50 * edge + 50 * np.sqrt(3.0) * np.array([-edge[1], edge[0]])
    lattice = np.vstack((along, across))
    offset = coords - coords[index]
    address = np.rint(offset @ np.linalg.inv(lattice))
    if np.linalg.matrix_rank(address) < 2 or np.any(
        np.linalg.norm(offset - address @ lattice, axis=1) > 0.5
    ):
        raise ValueError("spatial coordinates do not reconstruct a hexagonal Visium array")
    return coords


def read_sections(path, entry, spec, features):
    import h5py
    import pandas as pd
    from anndata.io import read_elem

    from .schema import Section

    layer = spec.get("layer") or "X"
    key = "X" if layer == "X" else "layers/" + layer
    signature = _signature(path)
    with h5py.File(path, "r") as handle:
        node = handle[key]
        if isinstance(node, h5py.Dataset):
            form = "array"
            shape = tuple(node.shape)
        else:
            form = node.attrs["encoding-type"]
            if form not in {"csr_matrix", "csc_matrix"}:
                raise ValueError(f"Unsupported H5AD matrix encoding {form}")
            shape = tuple(map(int, node.attrs["shape"]))
        var = handle["var"]
        names = pd.Index(read_elem(var[var.attrs["_index"]])).astype(str)
        if not names.is_unique:
            raise ValueError("Duplicate feature IDs")
        features = features or names.tolist()[: spec.get("feature_limit")]
        cols = names.get_indexer(features)
        if (cols < 0).any():
            raise ValueError("Missing features across sections")
        if len(set(cols)) != len(cols):
            raise ValueError("Duplicate requested features")
        n = shape[0]
        obs = handle["obs"]

        def obs_col(name):
            return np.asarray(read_elem(obs[name]))

        def field(name, default=None):
            if name in entry:
                return np.full(n, entry[name], dtype=object)
            column = spec.get("columns", {}).get(name, name)
            if column in obs:
                values = obs_col(column)
                if name in {"slide_id", "batch"}:
                    return np.array([None if pd.isna(v) else str(v) for v in values], dtype=object)
                if pd.isna(values).any():
                    raise ValueError(f"Missing values in {column}")
                return values.astype(str)
            return np.full(n, default, dtype=object)

        frame = pd.DataFrame(
            {
                k: field(k, path.stem if k == "section_id" else None)
                for k in ("patient_id", "group", "section_id")
            }
        )
        if frame.isna().any().any():
            raise ValueError("Missing patient_id, group or section_id mapping")
        if spec.get("infer_section_ids") and spec["columns"]["section_id"] not in obs:
            indices = frame.groupby(["patient_id", "group"], sort=False).ngroup()
            frame["section_id"] = "section-" + indices.astype(str)
        slides, batches = field("slide_id"), field("batch")
        stages = (
            np.column_stack([obs_col(k) for k in spec["stage_columns"]])
            if spec.get("stage_columns")
            else None
        )
        domain = obs_col(spec["domain"]) if spec.get("domain") else None
        library_size = (
            obs_col(spec["library_size_column"]) if spec.get("library_size_column") else None
        )
        C = None
        if spec.get("covariates"):
            C = np.column_stack([obs_col(k) for k in spec["covariates"]]).astype(float)
        if spec.get("covariate_key"):
            stored: Any = read_elem(handle["obsm"][spec["covariate_key"]])
            C = np.asarray(stored.toarray() if sparse.issparse(stored) else stored, float)
        classic = spec.get("coordinate_calibration") == "classic_visium_100um"
        addresses = None
        if classic and ("array_row" in obs or "array_col" in obs):
            if "array_row" not in obs or "array_col" not in obs:
                raise ValueError("Visium calibration requires both array_row and array_col")
            addresses = np.column_stack((obs_col("array_row"), obs_col("array_col")))
        coords = (
            np.asarray(read_elem(handle["obsm"][spec.get("coord_key", "spatial")]), float)
            if addresses is None
            else None
        )
        output = []

        def unique(values, rows):
            choices = pd.unique(values[rows])
            if len(choices) != 1:
                raise ValueError("Acquisition metadata must be constant within section")
            return choices[0]

        mode = spec.get("matrix_mode", "backed")
        if mode not in {"backed", "memory"}:
            raise ValueError("matrix_mode must be backed or memory")
        for (pid, grp, sid), rows in frame.groupby(list(frame), sort=False).indices.items():
            stage = entry.get("stage")
            if stages is not None:
                pairs = np.unique(stages[rows], axis=0)
                if len(pairs) != 1:
                    raise ValueError("Split acquisition fields with distinct stage coordinates")
                stage = pairs[0].astype(float)
            view = H5Matrix(path, key, rows, cols, signature, shape, form)
            selected: Any = view
            if mode == "memory":
                block = view[:, :]
                selected = block.tocsc() if sparse.issparse(block) else block
            if addresses is not None:
                section_coords = classic_visium_coordinates(None, addresses[rows])
            else:
                assert coords is not None
                section_coords = coords[rows]
                if classic:
                    section_coords = classic_visium_coordinates(section_coords)
            output.append(
                Section(
                    str(pid),
                    str(grp),
                    section_coords,
                    selected,
                    None if C is None else C[rows],
                    unique(slides, rows),
                    unique(batches, rows),
                    stage,
                    str(sid),
                    None if domain is None else domain[rows],
                    None if library_size is None else library_size[rows],
                )
            )
    return features, output
