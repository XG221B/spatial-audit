"""Finite-volume, physical-band energy diagnostics on occupied square control volumes.

Rows of X are point measurements. A volume contains the arithmetic mean of the rows in
it. Edge conductance is shared-face length / centre distance = 1. Dividing flux by
volume area h² gives the reflecting (zero outward flux) finite-volume Laplacian.
"""

from dataclasses import dataclass, replace
from itertools import combinations

import numpy as np
from scipy import sparse
from scipy.spatial import ConvexHull, Delaunay, cKDTree
from scipy.special import ive

from .schema import Cohort


@dataclass
class Grid:
    ids: np.ndarray
    count: np.ndarray
    left: np.ndarray
    right: np.ndarray
    far_left: np.ndarray
    middle: np.ndarray
    far_right: np.ndarray
    laplacian: sparse.csr_matrix
    width: float
    centers: np.ndarray
    voronoi: bool = False

    def average(self, values):
        """Bin one small feature block without forming a dense observation matrix."""
        if sparse.issparse(values):
            values = values.toarray()
        array = np.asarray(values, dtype=float)
        if array.ndim == 1:
            array = array[:, None]
        return np.column_stack(
            [
                np.bincount(self.ids, weights=column, minlength=len(self.count)) / self.count
                for column in array.T
            ]
        )


def build_grid(coords, finest_um):
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) == 0:
        raise ValueError("coordinates must be a non-empty n × 2 array")
    if not np.isfinite(coords).all() or not np.isfinite(finest_um) or finest_um <= 0:
        raise ValueError("finite coordinates and positive wavelengths are required")
    points = coords - coords.min(axis=0)
    # Unique sorted sites make the estimated pitch invariant to row permutations
    # and repeated observations at a location. The deterministic subsample is
    # used only for the pitch; graph construction and energy use all observations.
    sites = np.unique(points, axis=0)
    pitch = 0.0
    if len(sites) > 1:
        step = max(1, len(sites) // 4096)
        distance = cKDTree(sites).query(sites[::step], k=2)[0][:, 1]
        pitch = float(np.median(distance))
    width = max(finest_um / 8, pitch)
    scaled = points / width
    if np.max(scaled) > 2**45 or (scaled.max(axis=0) + 4).prod() >= 2**62:
        raise ValueError("coordinates exceed representable grid extent at this wavelength")
    # The half-cell shift places the minimum site at its volume centre.
    ix, iy = np.floor(scaled + 0.5).astype(np.int64).T
    stride = int(ix.max()) + 3
    keys, ids, count = np.unique(iy * stride + ix, return_inverse=True, return_counts=True)
    one_a, one_b, two_a, two_b, middle = [], [], [], [], []
    for delta in (1, stride):
        near = np.searchsorted(keys, keys + delta)
        first = np.flatnonzero(keys[np.minimum(near, len(keys) - 1)] == keys + delta)
        one_a.append(first)
        one_b.append(near[first])
        far = np.searchsorted(keys, keys[first] + 2 * delta)
        exists = keys[np.minimum(far, len(keys) - 1)] == keys[first] + 2 * delta
        two_a.append(first[exists])
        middle.append(near[first[exists]])
        two_b.append(far[exists])
    a, b, c, d = (np.concatenate(x) for x in (one_a, one_b, two_a, two_b))
    neighbours = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(len(keys), len(keys))
    )
    return Grid(
        ids,
        count,
        a,
        b,
        c,
        np.concatenate(middle),
        d,
        (sparse.diags(np.asarray(neighbours.sum(1)).ravel()) - neighbours).tocsr() / width**2,
        width,
        np.column_stack((keys % stride, keys // stride)).astype(float) * width,
    )


def recommend_wavelengths(pitch_um, max_span_um):
    """Log-space up to three physical bands between the axial limit and tissue span."""
    if not np.isfinite([pitch_um, max_span_um]).all() or min(pitch_um, max_span_um) <= 0:
        raise ValueError("Positive finite pitch and footprint span are required")
    shortest = np.pi * pitch_um
    if shortest > max_span_um * (1 + 1e-10):
        raise ValueError("Footprint cannot resolve a physical band at the measured pitch")
    if np.isclose(shortest, max_span_um):
        return (float(max_span_um),)
    return tuple(map(float, np.geomspace(shortest, max_span_um, 3)))


def _minimum_width(points):
    if len(points) < 3 or np.linalg.matrix_rank(points - points.mean(0)) < 2:
        return 0.0
    hull = points[ConvexHull(points).vertices]
    edges = np.roll(hull, -1, axis=0) - hull
    normals = np.column_stack((-edges[:, 1], edges[:, 0])) / np.linalg.norm(edges, axis=1)[:, None]
    return float(min(np.ptp(hull @ normal) for normal in normals))


def _clip_cell(polygon, normal, midpoint):
    """Clip a convex polygon by a bisector, retaining the nearer half-plane."""
    signed = (polygon - midpoint) @ normal
    inside = signed <= 0
    if inside.all():
        return polygon
    vertices = []
    for i, point in enumerate(polygon):
        if inside[i] != inside[i - 1]:
            weight = signed[i - 1] / (signed[i - 1] - signed[i])
            vertices.append(polygon[i - 1] + weight * (point - polygon[i - 1]))
        if inside[i]:
            vertices.append(point)
    return np.array(vertices)


def build_voronoi(coords):
    """No-flux FV graph on Voronoi cells clipped to the convex hull of sampled sites."""
    coords = np.asarray(coords, float)
    # Preserve first-site order so seeded probes follow sites under rigid motion.
    _, first, inverse, counts = np.unique(
        coords, axis=0, return_index=True, return_inverse=True, return_counts=True
    )
    order = np.argsort(first)
    ids = np.argsort(order)[inverse]
    points = coords[first[order]] - coords[first[order][0]]
    if len(points) < 3 or np.linalg.matrix_rank(points) < 2:
        raise ValueError("Voronoi graph requires at least three noncollinear sites")
    hull = points[ConvexHull(points).vertices]
    offsets, indices = Delaunay(points).vertex_neighbor_vertices
    pitch = float(np.median(cKDTree(points).query(points, k=2)[0][:, 1]))
    left, right, conductance, areas = [], [], [], []
    triples = []
    for i, point in enumerate(points):
        neighbors = indices[offsets[i] : offsets[i + 1]]
        cell = hull
        for j in neighbors:
            cell = _clip_cell(cell, points[j] - point, (point + points[j]) / 2)
        local = cell - point
        next_vertex = np.roll(local, -1, axis=0)
        areas.append(
            abs(np.sum(local[:, 0] * next_vertex[:, 1] - local[:, 1] * next_vertex[:, 0])) / 2
        )
        adjacent = []
        for j in neighbors:
            direction = points[j] - point
            dist = np.linalg.norm(direction)
            on_face = np.abs(local @ (direction / dist) - dist / 2) <= 1e-9 * dist
            lengths = np.linalg.norm(next_vertex - local, axis=1)
            face = lengths[on_face & np.roll(on_face, -1)].sum()
            if face <= 1e-10 * dist:
                continue
            adjacent.append(j)
            if j > i:
                left.append(i)
                right.append(j)
                conductance.append(face / dist)
        for a, b in combinations(adjacent, 2):
            if np.linalg.norm(points[a] + points[b] - 2 * point) <= 1e-8 * pitch:
                triples.append((a, i, b))
    a, b = np.asarray(left, int), np.asarray(right, int)
    matrix = sparse.csr_matrix(
        (np.tile(conductance, 2), (np.r_[a, b], np.r_[b, a])), shape=(len(points), len(points))
    )
    lap = sparse.diags(1 / np.asarray(areas)) @ (
        sparse.diags(np.asarray(matrix.sum(1)).ravel()) - matrix
    )
    far, mid, end = np.asarray(triples, int).reshape(-1, 3).T
    return Grid(ids, counts[order], a, b, far, mid, end, lap.tocsr(), pitch, points, True)


def coefficients(wavelengths, width, tol=2e-10):
    """Expand two heat kernels exactly in the Chebyshev basis, truncating their tails."""
    waves = np.asarray(wavelengths, dtype=float)
    if (
        waves.ndim != 1
        or not waves.size
        or not np.isfinite(waves).all()
        or (waves <= 0).any()
        or (np.diff(waves) <= 0).any()
    ):
        raise ValueError("specify one or more strictly increasing positive finite wavelengths")
    # Enclose all spectra of the four-neighbour FV stencil in [0, 8/h²].
    scaled_time = 4 * np.log(4) / (3 * (2 * np.pi * width / waves) ** 2)
    # Beyond this point the scaled Bessel tail is negligible; limit work for
    # unreasonable scales with an explicit diagnostic instead of a silent cutoff.
    proposed = float(16 + 14 * np.sqrt(1 + 4 * scaled_time.max()))
    if not np.isfinite(proposed) or proposed > 2048:
        raise ValueError("wavelength range requires >2048 Chebyshev terms; choose closer bands")
    k = np.arange(int(np.ceil(proposed)) + 1)
    coef = ive(k[:, None], scaled_time) - ive(k[:, None], 4 * scaled_time)
    coef[1:] *= 2
    coef *= ((-1.0) ** k)[:, None]
    last = np.where(np.max(abs(coef), axis=1) > tol / len(k))[0]
    coef = coef[: max(2, int(last[-1]) + 1) if last.size else 2].copy()
    coef[0] -= ((-1.0) ** np.arange(len(coef))) @ coef
    return coef


def filter_bands(grid, values, weights, reference_width=None):
    """Shared sparse Chebyshev recurrence for all wavelengths and feature columns."""
    values = np.asarray(values, dtype=float)
    previous = values - values.mean(axis=0)
    width = grid.width if reference_width is None else reference_width
    current = grid.laplacian @ previous * (width**2 / 4) - previous
    out = [weights[0, b] * previous + weights[1, b] * current for b in range(weights.shape[1])]
    for row in weights[2:]:
        following = (width**2 / 2) * (grid.laplacian @ current) - 2 * current - previous
        for b in range(weights.shape[1]):
            out[b] += row[b] * following
        previous, current = current, following
    return out


def basis(x, original):
    """Normalized covariate columns make the projector invariant to basis units."""
    norms = np.linalg.norm(x, axis=0)
    keep = norms > 1e-10 * np.maximum(np.linalg.norm(original, axis=0), np.finfo(float).tiny)
    if not keep.any():
        return np.empty((len(x), 0))
    u, s, _ = np.linalg.svd(x[:, keep] / norms[keep], full_matrices=False)
    return u[:, s > max(x.shape) * np.finfo(float).eps * s[0]]


def nugget_variance(grid, x, mask=None):
    """Two-lag variogram intercept for locally smooth signal plus uncorrelated noise."""
    a, m, b = grid.far_left, grid.middle, grid.far_right
    if mask is not None:
        keep = mask[a] & mask[m] & mask[b]
        a, m, b = a[keep], m[keep], b[keep]
    if not len(a):
        return None
    # Equivalent to (4*gamma_1-gamma_2)/3 for the same matched straight
    # paths at both lags. This form avoids cancellation for linear trends.
    return np.mean(((x[a] - x[m]) - (x[m] - x[b])) ** 2, axis=0) / 6


def component_support(grid, coords, waves):
    """Conservative 2-D geometric screen on each occupied connected component."""
    count, labels = sparse.csgraph.connected_components(grid.laplacian, directed=False)
    obs_labels = labels[grid.ids]
    lo = np.full((count, 2), np.inf)
    hi = np.full((count, 2), -np.inf)
    np.minimum.at(lo, obs_labels, coords)
    np.maximum.at(hi, obs_labels, coords)
    footprints = np.min(hi - lo, axis=1)
    # Axis-aligned bounds alone can falsely accept a narrow diagonal ribbon.
    # Also measure the transverse span on the occupied grid's minor PCA axis.
    sorted_idx = np.argsort(labels, kind="stable")
    splits = np.searchsorted(labels[sorted_idx], np.arange(count + 1))
    for k in np.flatnonzero(footprints >= min(waves) * (1 - 1e-10)):
        points = grid.centers[sorted_idx[splits[k] : splits[k + 1]]]
        centered = points - points.mean(axis=0)
        _, axes = np.linalg.eigh(centered.T @ centered)
        footprints[k] = min(footprints[k], np.ptp(centered @ axes[:, 0]))
    if grid.voronoi:
        footprints = np.array([_minimum_width(grid.centers[labels == k]) for k in range(count)])
    masks = []
    counts = []
    for wave in waves:
        keep = (wave >= np.pi * grid.width * (1 - 1e-10)) & (footprints >= wave * (1 - 1e-10))
        masks.append(keep[labels])
        counts.append(int(np.sum(keep)))
    return masks, counts, count, float(max(footprints, default=0))


def transfer_diagnostics(grid, weights, qs, masks, component_counts, reference_width=None):
    """Trace of in-band white noise, transfer participation ratio and fit saturation."""
    n = len(grid.count)
    exact = n <= 144
    probes = np.eye(n) if exact else np.random.default_rng(819).choice([-1.0, 1.0], (n, 48))
    filtered = filter_bands(grid, probes, weights, reference_width)
    records = []
    for h, q, mask, n_components in zip(filtered, qs, masks, component_counts, strict=True):
        size = int(mask.sum())
        record = {
            "supported_volumes": size,
            "supported_components": n_components,
            "projection_rank": q.shape[1],
            "effective_spectral_dimension": None,
            "projected_white_fraction": None,
            "saturated": False,
            "residual_white_gain": None,
            "residual_white_gain_se": None,
        }
        if size:
            z = h[mask]
            total = np.sum(z * z)
            if total > np.finfo(float).tiny:
                gram = z.T @ z
                fourth = (
                    np.sum(gram * gram)
                    if exact
                    else np.sum(np.triu(gram, 1) ** 2) * 2 / (len(gram) * (len(gram) - 1))
                )
                divisor = 1 if exact else probes.shape[1]
                trace = total / divisor
                dim = float(
                    np.clip(
                        trace**2 / max(fourth, np.finfo(float).tiny), 1, max(1, size - n_components)
                    )
                )
                projected = q.T @ z
                fraction = float(np.clip(np.sum(projected * projected) / total, 0, 1))
                residual = z - q @ projected
                record.update(
                    residual_white_gain_se=(
                        0.0
                        if exact
                        else float(np.std(np.mean(residual**2, axis=0), ddof=1) / np.sqrt(divisor))
                    ),
                    effective_spectral_dimension=dim,
                    projected_white_fraction=fraction,
                    saturated=bool(q.shape[1] >= 0.8 * dim and fraction >= 0.8),
                    residual_white_gain=float(np.sum(residual * residual) / (size * divisor)),
                )
        records.append(record)
    return records


def trace_gain(grid, weights, q, reference_width=None):
    """Mean fine-band white energy after projecting out q."""
    n = len(grid.count)
    k = weights.shape[1]
    qs = [q] + [np.empty((n, 0)) for _ in range(k - 1)]
    return transfer_diagnostics(
        grid, weights, qs, [np.ones(n, bool)] * k, [1] * k, reference_width
    )[0]["residual_white_gain"]


def section_statistics(
    section,
    waves,
    block=16,
    shared_weights=None,
    *,
    nugget_mode="raw",
    graph_mode="grid",
    leave_one_out_tic=False,
):
    if nugget_mode not in {"raw", "vst"} or graph_mode not in {"grid", "voronoi"}:
        raise ValueError("Unknown graph or nugget mode")
    grid = (
        build_grid(section.coords, waves[0])
        if graph_mode == "grid"
        else build_voronoi(section.coords)
    )
    reference_width = waves[0] / 8
    if grid.voronoi:
        # Gershgorin bound for the area-normalized flux operator.
        reference_width = np.sqrt(4 / max(grid.laplacian.diagonal()))
        shared_weights = None
    weight = coefficients(waves, reference_width) if shared_weights is None else shared_weights
    if not leave_one_out_tic:
        return _section_values(section, waves, block, grid, weight, reference_width, nugget_mode)
    if section.C is None or not section.C.shape[1]:
        raise ValueError(
            "Leave-one-ion-out projection requires total ion current in covariate column zero"
        )
    columns, diagnostics = [], []
    for j in range(section.X.shape[1]):
        ion = section.X[:, j : j + 1]
        ion = ion.toarray() if sparse.issparse(ion) else np.asarray(ion)
        covariates = np.array(section.C, dtype=float, copy=True)
        tic = covariates[:, 0]
        if np.any(ion < 0) or np.any(ion[:, 0] > tic + 1e-6 * np.abs(tic)):
            raise ValueError("Ion abundance must be nonnegative and included in total ion current")
        covariates[:, 0] = np.maximum(tic - ion[:, 0], 0)
        values, info = _section_values(
            replace(section, X=ion, C=covariates),
            waves,
            block,
            grid,
            weight,
            reference_width,
            nugget_mode,
        )
        columns.append(values)
        diagnostics.append(info)
    info = diagnostics[0].copy()
    info["feature_bands"] = [d["bands"] for d in diagnostics]
    info["bands"] = [dict(b) for b in info["bands"]]
    for i, band in enumerate(info["bands"]):
        for key in band:
            vals = [d["bands"][i][key] for d in diagnostics]
            band[key] = (
                all(vals) if key == "saturated" else vals[0] if len(set(vals)) == 1 else None
            )
    for i, band in enumerate(info["bands"]):
        band["saturated_features"] = sum(bool(d["bands"][i]["saturated"]) for d in diagnostics)
    info["projection_ranks"] = [b["projection_rank"] for b in info["bands"]]
    info["fine_noise_gain"] = info["bands"][0]["residual_white_gain"]
    return np.column_stack(columns), info


def _section_values(section, waves, block, grid, weight, reference_width, nugget_mode):
    k = len(waves)
    n = len(grid.count)
    masks, counts, n_components, max_span = component_support(grid, section.coords, waves)
    c = (
        grid.average(section.C)
        if section.C is not None and section.C.shape[1]
        else np.empty((n, 0))
    )
    qs = [np.empty((int(m.sum()), 0)) for m in masks]
    qc = qs[0]
    if c.shape[1]:
        for i, fc in enumerate(filter_bands(grid, c, weight, reference_width)):
            mask = masks[i]
            if mask.any():
                local = c[mask] - c[mask].mean(axis=0)
                qs[i] = basis(fc[mask], local)
                if i == 0:
                    qc = basis(local, local)
    records = transfer_diagnostics(grid, weight, qs, masks, counts, reference_width)
    gain = records[0]["residual_white_gain"]
    p = section.X.shape[1]
    values = np.full((2 * k + 4, p), np.nan)
    fine_ok = gain is not None and not records[0]["saturated"]
    for start in range(0, p, block):
        idx = slice(start, start + block)
        x = grid.average(section.X[:, idx])
        if nugget_mode == "vst":
            if np.any(x < 0):
                raise ValueError("Anscombe stabilization requires nonnegative volume means")
            fine_x = 2 * np.sqrt(x + 3 / 8)
            fine_band = filter_bands(grid, fine_x, weight[:, :1], reference_width)[0]
        else:
            fine_x = x
        bands = filter_bands(grid, x, weight, reference_width)
        for i, mask in enumerate(masks):
            if not mask.any():
                continue
            z = bands[i][mask]
            # A constant offset is a Neumann zero mode, so it must not change
            # which features qualify as measurable. Equal values on the
            # entire supported mask also imply zero within-component energy.
            local_x = x[mask]
            centered_x = local_x - local_x.mean(axis=0)
            input_norm = np.linalg.norm(centered_x, axis=0)
            near_zero = (input_norm == 0) | (np.linalg.norm(z, axis=0) <= 1e-12 * input_norm)
            z[:, near_zero] = 0
            fit = qs[i].T @ z
            values[i, idx] = np.mean(z * z, axis=0)
            if c.shape[1] and not records[i]["saturated"]:
                values[k + i, idx] = np.sum(fit * fit, axis=0) / len(z)
            if i == 0 and fine_ok:
                fine_z = fine_band[mask] if nugget_mode == "vst" else z
                residual = fine_z - qs[i] @ (qs[i].T @ fine_z)
                values[-4, idx] = np.mean(residual * residual, axis=0)
                values[-2, idx] = fine_x[mask].mean(axis=0)
        if fine_ok:
            mask = masks[0]
            if qc.shape[1]:
                fine_x[mask] -= qc @ (qc.T @ (fine_x[mask] - fine_x[mask].mean(axis=0)))
            eta = nugget_variance(grid, fine_x, mask)
            if eta is not None:
                values[-1, idx] = eta * gain
                values[-3, idx] = values[-4, idx] - values[-1, idx]
    info = {
        "volumes": n,
        "chebyshev_degree": len(weight) - 1,
        "spectral_bound_per_um2": 8 / reference_width**2,
        "pitch_um": float(grid.width),
        "edges": len(grid.left),
        "connected_components": n_components,
        "largest_2d_span_um": max_span,
        "fov_span_um": np.ptp(section.coords, axis=0).tolist(),
        "nyquist_wavelength_um": 2 * grid.width,
        "axial_limit_um": np.pi * grid.width,
        "fine_noise_gain": gain,
        "trace_probes": 0 if n <= 144 else 48,
        "trace_method": "exact" if n <= 144 else "seeded Rademacher",
        "projection_ranks": [q.shape[1] for q in qs],
        "two_step_edges": int(np.sum(masks[0][grid.far_left])),
        "occupancy_range": [int(grid.count.min()), int(grid.count.max())],
        "bands": records,
    }
    return values, info


def available_mean(values, axis=0):
    array = np.asarray(values, float)
    counts = np.sum(np.isfinite(array), axis=axis)
    return np.divide(
        np.nansum(array, axis=axis), counts, out=np.full(np.shape(counts), np.nan), where=counts > 0
    )


def finite_value(value):
    return float(value) if np.isfinite(value) else None


def spectral_audit(
    cohort: Cohort,
    wavelengths=(250.0, 600.0, 1400.0),
    *,
    nugget_mode="raw",
    graph_mode="grid",
    leave_one_out_tic=False,
):
    if isinstance(wavelengths, str) and wavelengths == "auto":
        geometries = [np.unique(s.coords, axis=0) for s in cohort.sections]
        pitches = [np.median(cKDTree(x).query(x, k=2)[0][:, 1]) for x in geometries if len(x) > 1]
        if not pitches:
            raise ValueError("Cannot estimate pitch from single-site sections")
        wavelengths = recommend_wavelengths(max(pitches), min(map(_minimum_width, geometries)))
    waves = tuple(float(x) for x in wavelengths)
    if not waves:
        raise ValueError("At least one wavelength is required")
    shared_weights = coefficients(waves, max(waves[0] / 8, 1e-100))
    patients: dict = {}
    diagnostics = []
    for section in cohort.sections:
        values, info = section_statistics(
            section,
            waves,
            shared_weights=shared_weights,
            nugget_mode=nugget_mode,
            graph_mode=graph_mode,
            leave_one_out_tic=leave_one_out_tic,
        )
        diagnostics.append(
            dict(
                section_id=section.section_id,
                patient_id=section.patient_id,
                group=section.group,
                has_covariates=bool(section.C is not None and section.C.shape[1]),
                **info,
            )
        )
        patients.setdefault((section.patient_id, section.group), []).append(values)
    if not patients:
        raise ValueError("cohort has no sections")
    labels = np.array([key[1] for key in patients])
    groups = sorted(set(labels))
    if len(groups) != 2:
        raise ValueError("spectral group contrast requires two patient groups")
    data = np.array([available_mean(rows) for rows in patients.values()])
    bands = []
    for i, wave in enumerate(waves):
        good = [d for d in diagnostics if d["bands"][i]["supported_volumes"]]
        interpretable = [d for d in good if d["has_covariates"] and not d["bands"][i]["saturated"]]
        subset = [
            key for key, rows in patients.items() if np.isfinite(data_for_band(rows, i)).any()
        ]
        share = None
        valid = 0
        if {g for _, g in subset} == set(groups):
            # Estimate numerator and denominator on the identical section subset.
            pairs = np.array([data_for_band(patients[k], i) for k in subset])
            totals, explained = available_mean(pairs)
            eligible = totals > np.finfo(float).tiny
            if leave_one_out_tic:
                covered = np.isfinite(pairs[:, 1])
                for group in groups:
                    eligible &= covered[[g == group for _, g in subset]].any(axis=0)
            valid = int(eligible.sum())
            if valid:
                share = float(np.clip(np.mean(explained[eligible] / totals[eligible]), 0, 1))
        bands.append(
            {
                "wavelength_um": wave,
                "share": share,
                "estimable_features": valid,
                "supported_sections": len(good),
                "supported_patients": len({d["patient_id"] for d in good}),
                "supported_patients_by_group": {
                    str(g): len({d["patient_id"] for d in good if d["group"] == g}) for g in groups
                },
                "interpretable_sections": len(interpretable),
                "interpretable_patients": len({p for p, _ in subset}),
                "full_coverage": len(interpretable) == len(diagnostics)
                and not any(d["bands"][i].get("saturated_features", 0) for d in good),
                "saturated_sections": sum(d["bands"][i]["saturated"] for d in good),
            }
        )
    high = any(b["full_coverage"] and b["share"] is not None and b["share"] >= 0.5 for b in bands)
    uncertain = any(
        not b["full_coverage"] or b["share"] is None or b["share"] >= 0.25 for b in bands
    )
    composition = {
        "bands": bands,
        "sections": diagnostics,
        "total_patients": len({p for p, _ in patients}),
        "total_patients_by_group": {str(g): int(np.sum(labels == g)) for g in groups},
        "total_sections": len(diagnostics),
        "covariate_sections": sum(d["has_covariates"] for d in diagnostics),
        "graph_mode": graph_mode,
        "leave_one_out_tic": leave_one_out_tic,
        "chebyshev_degree": max(d["chebyshev_degree"] for d in diagnostics),
        "shared_spectral_bound_per_um2": max(d["spectral_bound_per_um2"] for d in diagnostics),
        "status": "red" if high else "yellow" if uncertain else "green",
        "message": "Feature-averaged explained energy on matched supporting sections, equally weighted within patient and across patients. A connected 2-D footprint must span a wavelength, with wavelength at least pi times the occupied square-grid pitch. Regridding can alter native geometry. Missing covariates, saturated fits and partial coverage limit interpretation.",
    }
    if graph_mode == "voronoi":
        composition["message"] = (
            "Voronoi control volumes are clipped to the convex hull with reflecting boundaries. "
            "Energy shares use equal site weights. Hull edges may bridge gaps in sampling. "
            "Support requires a two-dimensional span of at least one wavelength."
        )
    a, b = (available_mean(data[labels == g]) for g in groups)
    features = []
    for j, name in enumerate(cohort.features):
        raw = b[-4, j] - a[-4, j] if np.isfinite(data[:, -4, j]).all() else np.nan
        adj = b[-3, j] - a[-3, j] if np.isfinite(data[:, -3, j]).all() else np.nan
        mean_shift = b[-2, j] - a[-2, j]
        bias_shift = b[-1, j] - a[-1, j]
        scale = max(a[-4, j], b[-4, j])
        reference = scale if nugget_mode == "vst" else max(a[0, j], b[0, j])
        material = abs(raw) > max(0.25 * scale, 1e-10 * max(reference, 1e-100))
        flag = bool(
            np.isfinite(adj)
            and material
            and raw * mean_shift > 0
            and raw * bias_shift > 0
            and abs(adj) <= 0.5 * abs(raw)
        )
        features.append(
            {
                "feature": str(name),
                "raw_delta": finite_value(raw),
                "debiased_delta": finite_value(adj),
                "raw_energy": [float(a[-4, j]), float(b[-4, j])] if np.isfinite(raw) else None,
                "corrected_energy": [float(a[-3, j]), float(b[-3, j])]
                if np.isfinite(adj)
                else None,
                "mean_delta": finite_value(mean_shift) if np.isfinite(raw) else None,
                "estimated_noise_delta": finite_value(bias_shift) if np.isfinite(adj) else None,
                "mean_noise_sensitive": flag,
            }
        )
    flagged = sum(bool(f["mean_noise_sensitive"]) for f in features)
    evaluated = sum(f["debiased_delta"] is not None for f in features)
    patient_energies = [
        {
            "patient_id": str(key[0]),
            "group": str(key[1]),
            "raw": list(map(finite_value, row[-4])),
            "corrected": list(map(finite_value, row[-3])),
        }
        for key, row in zip(patients, data, strict=True)
    ]
    incomplete = any(not d["two_step_edges"] or d["bands"][0]["saturated"] for d in diagnostics)
    nugget = {
        "nugget_mode": nugget_mode,
        "groups": list(map(str, groups)),
        "features": features,
        "flagged": flagged,
        "patient_energies": patient_energies,
        "evaluated": evaluated,
        "status": "red"
        if flagged and flagged / max(1, evaluated) >= 0.2 and not incomplete
        else "yellow"
        if flagged or evaluated < len(features) or incomplete
        else "green",
        "message": "Fine-band residual energy minus the estimated white-noise contribution is signed. Matched two-step spatial differences model locally smooth signal plus uncorrelated homoscedastic volume noise. Rough signal, heterogeneous noise and unequal occupancy affect this sensitivity screen.",
    }
    return composition, nugget


def data_for_band(section_rows, band):
    data = np.asarray(section_rows)
    k = (data.shape[1] - 4) // 2
    values = data[:, [band, k + band]]
    valid = np.isfinite(values[:, 1])
    return available_mean(np.where(valid[:, None, :], values, np.nan))
