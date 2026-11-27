"""Finite-volume, physical-band energy diagnostics on occupied square control volumes.

Rows of X are point measurements. A volume contains the arithmetic mean of the rows in
it. Edge conductance is shared-face length / centre distance = 1. Dividing flux by
volume area h² gives the reflecting (zero outward flux) finite-volume Laplacian.
"""

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
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


def coefficients(wavelengths, width, tol=2e-10):
    """Expand two heat kernels exactly in the Chebyshev basis, truncating their tails."""
    waves = np.asarray(wavelengths, dtype=float)
    if (
        waves.shape != (3,)
        or not np.isfinite(waves).all()
        or (waves <= 0).any()
        or (np.diff(waves) <= 0).any()
    ):
        raise ValueError("specify three strictly increasing positive finite wavelengths")
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
    out = [weights[0, b] * previous + weights[1, b] * current for b in range(3)]
    for row in weights[2:]:
        following = (width**2 / 2) * (grid.laplacian @ current) - 2 * current - previous
        for b in range(3):
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
    probes = np.eye(n) if exact else np.random.default_rng(819).choice([-1.0, 1.0], (n, 12))
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
    qs = [q, np.empty((n, 0)), np.empty((n, 0))]
    return transfer_diagnostics(
        grid, weights, qs, [np.ones(n, bool)] * 3, [1, 1, 1], reference_width
    )[0]["residual_white_gain"]


def section_statistics(section, waves, block=16, shared_weights=None):
    grid = build_grid(section.coords, waves[0])
    reference_width = waves[0] / 8
    weight = coefficients(waves, reference_width) if shared_weights is None else shared_weights
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
    values = np.full((10, p), np.nan)
    fine_ok = gain is not None and not records[0]["saturated"]
    for start in range(0, p, block):
        idx = slice(start, start + block)
        x = grid.average(section.X[:, idx])
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
                values[3 + i, idx] = np.sum(fit * fit, axis=0) / len(z)
            if i == 0 and fine_ok:
                residual = z - qs[i] @ fit
                values[6, idx] = np.mean(residual * residual, axis=0)
                values[8, idx] = x[mask].mean(axis=0)
        if fine_ok:
            mask = masks[0]
            if qc.shape[1]:
                x[mask] -= qc @ (qc.T @ (x[mask] - x[mask].mean(axis=0)))
            eta = nugget_variance(grid, x, mask)
            if eta is not None:
                values[9, idx] = eta * gain
                values[7, idx] = values[6, idx] - values[9, idx]
    info = {
        "volumes": n,
        "pitch_um": float(grid.width),
        "edges": len(grid.left),
        "connected_components": n_components,
        "largest_2d_span_um": max_span,
        "fov_span_um": np.ptp(section.coords, axis=0).tolist(),
        "nyquist_wavelength_um": 2 * grid.width,
        "axial_limit_um": np.pi * grid.width,
        "fine_noise_gain": gain,
        "trace_probes": 0 if n <= 144 else 12,
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


def spectral_audit(cohort: Cohort, wavelengths=(250.0, 600.0, 1400.0)):
    waves = tuple(float(x) for x in wavelengths)
    if len(waves) != 3:
        raise ValueError("three wavelengths are required")
    shared_weights = coefficients(waves, max(waves[0] / 8, 1e-100))
    patients: dict = {}
    diagnostics = []
    for section in cohort.sections:
        values, info = section_statistics(section, waves, shared_weights=shared_weights)
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
                "interpretable_patients": len(subset),
                "full_coverage": len(interpretable) == len(diagnostics),
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
        "total_patients": len(patients),
        "total_patients_by_group": {str(g): int(np.sum(labels == g)) for g in groups},
        "total_sections": len(diagnostics),
        "covariate_sections": sum(d["has_covariates"] for d in diagnostics),
        "chebyshev_degree": len(shared_weights) - 1,
        "shared_spectral_bound_per_um2": 8 / (waves[0] / 8) ** 2,
        "status": "red" if high else "yellow" if uncertain else "green",
        "message": "Feature-averaged explained energy on matched supporting sections, equally weighted within patient and across patients. A connected 2-D footprint must span a wavelength, with wavelength at least pi times the occupied square-grid pitch. Regridding can alter native geometry. Missing covariates, saturated fits and partial coverage limit interpretation.",
    }
    a, b = (available_mean(data[labels == g]) for g in groups)
    features = []
    for j, name in enumerate(cohort.features):
        raw = b[6, j] - a[6, j] if np.isfinite(data[:, 6, j]).all() else np.nan
        adj = b[7, j] - a[7, j] if np.isfinite(data[:, 7, j]).all() else np.nan
        mean_shift = b[8, j] - a[8, j]
        bias_shift = b[9, j] - a[9, j]
        scale = max(a[6, j], b[6, j])
        material = abs(raw) > max(0.25 * scale, 1e-10 * max(a[0, j], b[0, j], 1e-100))
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
                "raw_energy": [float(a[6, j]), float(b[6, j])] if np.isfinite(raw) else None,
                "corrected_energy": [float(a[7, j]), float(b[7, j])] if np.isfinite(adj) else None,
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
            "raw": list(map(finite_value, row[6])),
            "corrected": list(map(finite_value, row[7])),
        }
        for key, row in zip(patients, data, strict=True)
    ]
    incomplete = any(not d["two_step_edges"] or d["bands"][0]["saturated"] for d in diagnostics)
    nugget = {
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
    values = np.asarray(section_rows)[:, [band, 3 + band]]
    valid = np.isfinite(values[:, 1])
    return available_mean(np.where(valid[:, None, :], values, np.nan))
