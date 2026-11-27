"""Design diagnostics for independent biological units."""

from collections import Counter, defaultdict
from math import comb, floor, log2, prod

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.linalg import lsqr
from scipy.stats import binomtest, fisher_exact, rankdata

from .permutation import cohort_test_options, exchangeability_blocks


def exact_floor(
    n0: int, n1: int, m: int = 1, alpha: float = 0.05, *, design="independent", stratum_counts=None
) -> dict:
    if any(type(i) is not int or i < 1 for i in (n0, n1, m)) or not 0 < alpha < 1:
        raise ValueError("Positive integer counts and 0 < alpha < 1 required")
    if design == "paired":
        if n0 != n1:
            raise ValueError("Paired groups must contain equal numbers of patients")
        M = 2**n0
    elif design == "stratified":
        if not stratum_counts or any(
            type(a) is not int or type(b) is not int or min(a, b) < 0 or a + b == 0
            for a, b in stratum_counts
        ):
            raise ValueError("Valid per-stratum counts are required")
        if tuple(map(sum, zip(*stratum_counts))) != (n0, n1):
            raise ValueError("Stratum counts disagree with the group totals")
        M = prod(comb(a + b, b) for a, b in stratum_counts)
    elif design == "independent":
        M = comb(n0 + n1, n0)
    else:
        raise ValueError("Unknown permutation design")
    p = min(1.0, 2 / M)
    # Strict BH q<alpha; use rational arithmetic at exact boundaries.
    from fractions import Fraction

    rank = floor(min(Fraction(1), Fraction(2, M)) * m / Fraction(str(alpha))) + 1
    status = "red" if rank > m else ("yellow" if rank > 1 else "green")
    return {
        "status": status,
        "message": "Two-sided doubled-tail permutation resolution.",
        "design": design,
        "n0": n0,
        "n1": n1,
        "assignments": str(M),
        "p_floor": p,
        "features": m,
        "alpha": alpha,
        "minimum_bh_rank": rank if rank <= m else None,
        "bh_impossible": rank > m,
        "nominal_impossible": p >= alpha,
    }


def _stage(points, labels):
    x = np.asarray(points, float)
    y = np.asarray(labels, int)
    scale = np.std(x, axis=0)
    z = (x - x.mean(axis=0)) / np.where(scale > 0, scale, 1)
    z = np.column_stack((np.ones(len(z)), z))
    fit = linprog(
        np.zeros(3),
        A_ub=-z * (2 * y - 1)[:, None],
        b_ub=-np.ones(len(y)),
        bounds=[(None, None)] * 3,
        method="highs",
    )
    score = z @ fit.x if fit.success else z @ np.linalg.lstsq(z, y, rcond=None)[0]
    ranks = rankdata(score)
    n0, n1 = np.bincount(y)
    auc = (ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    rank_alias = np.linalg.matrix_rank(np.column_stack((z, y))) == np.linalg.matrix_rank(z)
    return {
        "linear_separation": bool(fit.success),
        "descriptive_auc": max(auc, 1 - auc),
        "linear_fixed_effect_alias": bool(rank_alias),
    }


def _joint_design(sections, groups):
    """Test estimability after jointly accounting for available acquisition main effects."""
    n = len(sections)
    blocks = [sparse.csr_matrix(np.ones((n, 1)))]
    fields = []
    for name in ("batch", "slide_id"):
        levels = [getattr(s, name) for s in sections]
        if any(v is None for v in levels):
            continue
        _, idx = np.unique(levels, return_inverse=True)
        blocks.append(sparse.csr_matrix((np.ones(n), (np.arange(n), idx))))
        fields.append(name)
    if all(s.stage is not None for s in sections):
        xy = np.asarray([s.stage for s in sections], float)
        xy -= xy.mean(axis=0)
        sd = xy.std(axis=0)
        blocks.append(sparse.csr_matrix(xy / np.where(sd > 0, sd, 1)))
        fields.append("stage")
    if not fields:
        return {
            "fields": [],
            "status": "missing",
            "relative_residual": None,
            "fixed_effect_alias": None,
        }
    y = np.asarray([s.group == groups[1] for s in sections], float)
    design = sparse.hstack(blocks, format="csr")
    fit = lsqr(design, y, atol=1e-12, btol=1e-12)
    residual = np.linalg.norm(y - design @ fit[0]) / np.linalg.norm(y - y.mean())
    converged = fit[1] in (0, 1, 2, 4, 5)
    return {
        "fields": fields,
        "status": "observed",
        "relative_residual": float(residual),
        "fixed_effect_alias": bool(residual < 1e-9) if converged else None,
    }


_CATEGORICAL_DRAWS = 999
_CATEGORICAL_SEED = 819


def _categorical_mc(visits, labels, paired, draws=_CATEGORICAL_DRAWS):
    """Shuffle whole patient profiles, preserving multi-level visits and paired membership."""
    levels = {v: i for i, v in enumerate(sorted(set.union(*visits.values())))}
    rows, cols, data = [], [], []
    for i, observed in enumerate(visits.values()):
        for v in observed:
            rows.append(i)
            cols.append(levels[v])
            data.append(1 / len(observed))
    profiles = sparse.csr_matrix((data, (cols, rows)), shape=(len(levels), len(visits)))
    y = np.asarray([g == labels[1] for _, g in visits], float)
    sizes = np.asarray(profiles.sum(axis=1)).ravel()
    expected = sizes * y.mean()
    variance = expected * (1 - y.mean())

    def chisq(assignments):
        return np.sum((profiles @ assignments - expected[:, None]) ** 2 / variance[:, None], axis=0)

    observed_stat = chisq(y[:, None])[0]
    rng = np.random.default_rng(_CATEGORICAL_SEED)
    if paired:
        blocks = exchangeability_blocks(y, "paired", [p for p, _ in visits])
        pairs = np.array([rows for rows, _ in blocks])
    exceed = 0
    for start in range(0, draws, 32):
        count = min(32, draws - start)
        samples = np.tile(y[:, None], (1, count))
        if paired:
            flips = rng.integers(2, size=(len(pairs), count))
            samples[pairs[:, 0]] = flips
            samples[pairs[:, 1]] = 1 - flips
        else:
            samples = rng.permuted(samples, axis=0)
        exceed += int((chisq(samples) >= observed_stat - 1e-10).sum())
    return (exceed + 1) / (draws + 1)


def _patient_coordinates(sections, points, labels, section_diagnostics):
    visits: dict = defaultdict(list)
    for s, xy in zip(sections, points, strict=True):
        visits[s.patient_id, s.group].append(xy)
    # Center before averaging so repeated identical coordinates stay exactly constant.
    x = np.array([v[0] + np.mean(np.asarray(v) - v[0], axis=0) for v in visits.values()])
    y = np.array([g == labels[1] for _, g in visits], int)
    pids = np.array([p for p, _ in visits])
    separation = np.linalg.norm(x[y == 0].mean(0) - x[y == 1].mean(0))
    auc = None
    # Constant centroids give only the training-fold group proportion.
    if np.any(x != x[0]):
        predictions = np.zeros(len(y))
        # Each fold excludes all conditions and sections from the held-out patient.
        for pid in set(pids):
            train = pids != pid
            center, scale = x[train].mean(0), x[train].std(0)
            z = (x - center) / np.where(scale > 0, scale, 1)
            z = np.column_stack((np.ones(len(z)), z))
            predictions[~train] = z[~train] @ np.linalg.lstsq(z[train], y[train], rcond=None)[0]
        n0, n1 = np.bincount(y)
        auc = float((rankdata(predictions)[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1))
    patient_fit = section_diagnostics if len(x) == len(sections) else _stage(x, y)
    return {
        "patient_centroid_separation_um": float(separation),
        "patient_centroid_linear_separation": patient_fit["linear_separation"],
        "patient_descriptive_auc": patient_fit["descriptive_auc"],
        "lopo_cv_auc": auc,
    }


def design_audit(cohort):
    labels = sorted({s.group for s in cohort.sections})
    details: dict = {}
    reasons = []
    warning = False
    for field in ("slide_id", "batch"):
        # Each patient has total weight one, including patients spanning several slides.
        visits: dict = defaultdict(set)
        missing = any(getattr(s, field) is None for s in cohort.sections)
        for s in cohort.sections:
            if getattr(s, field) is not None:
                visits[s.patient_id, s.group].add(getattr(s, field))
        if missing or not visits:
            missing_units: dict = defaultdict(bool)
            for s in cohort.sections:
                missing_units[s.patient_id, s.group] |= getattr(s, field) is None
            p_missing = None
            if visits:
                if cohort.design == "paired":
                    differences = [
                        int(missing_units[p, labels[1]]) - int(missing_units[p, labels[0]])
                        for p in {p for p, _ in missing_units}
                    ]
                    discordant = sum(abs(v) for v in differences)
                    p_missing = (
                        binomtest(differences.count(1), discordant).pvalue if discordant else 1.0
                    )
                else:
                    counts = [
                        [
                            sum(g == group and v == flag for (_, g), v in missing_units.items())
                            for flag in (False, True)
                        ]
                        for group in labels
                    ]
                    p_missing = float(fisher_exact(counts).pvalue)
            details[field] = {
                "status": "missing_or_partial",
                "missing_fraction": sum(getattr(s, field) is None for s in cohort.sections)
                / len(cohort.sections),
                "missingness_group_p": p_missing,
                "missingness_unit": "patient-condition with any missing section",
            }
            warning = True
            if cohort.design == "paired":
                complete = {p for p, _ in visits if all((p, group) in visits for group in labels)}
                visits = {key: levels for key, levels in visits.items() if key[0] in complete}
            if len({g for _, g in visits}) < 2:
                continue
        table: dict = defaultdict(lambda: np.zeros(2))
        sizes: Counter = Counter()
        for (_, g), levels in visits.items():
            for level in levels:
                table[level][labels.index(g)] += 1 / len(levels)
                sizes[level] += 1
        t = np.array(list(table.values()))
        prob = t / t.sum()
        baseline = prob.sum(axis=0)
        independent = prob.sum(axis=1)[:, None] * baseline
        mi = sum(v * log2(v / b) for v, b in zip(prob.flat, independent.flat) if v > 0)
        entropy = -sum(v * log2(v) for v in baseline)
        separate = bool(np.all((t > 0).sum(axis=1) == 1))
        fisher = None
        if (
            t.shape == (2, 2)
            and cohort.design != "paired"
            and all(len(v) == 1 for v in visits.values())
        ):
            fisher = float(fisher_exact(t.astype(int)).pvalue)
        details[field] = {
            "status": "observed",
            "label": cohort.provenance.get("batch_label", "acquisition batch")
            if field == "batch"
            else "slide",
            "levels": len(t),
            "normalized_mi": mi / entropy,
            "complete_alias": None if missing else separate,
            "observed_levels_group_specific": separate,
            "fisher_p": fisher,
            "permutation_p": _categorical_mc(visits, labels, cohort.design == "paired")
            if len(t) > 2
            else None,
            "permutation_draws": _CATEGORICAL_DRAWS if len(t) > 2 else None,
            "permutation_seed": _CATEGORICAL_SEED if len(t) > 2 else None,
            "missing_fraction": 0.0,
            "all_levels_patient_specific": max(sizes.values()) == 1,
            "observed_patient_conditions": len(visits),
            **details.get(field, {}),
        }
        if separate and not missing:
            reasons.append(field)
        warning |= mi / entropy >= 0.5
    points = []
    flags = {}
    for field in ("stage", "bbox"):
        measured = all(s.stage is not None for s in cohort.sections)
        if field == "stage" and not measured:
            flags[field] = {"status": "missing"}
            warning = True
            continue
        xy = [
            s.stage if field == "stage" else (s.coords.min(axis=0) + s.coords.max(axis=0)) / 2
            for s in cohort.sections
        ]
        yy = [labels.index(s.group) for s in cohort.sections]
        diag = _stage(xy, yy)
        diag.update(_patient_coordinates(cohort.sections, xy, labels, diag))
        flags[field] = diag
        if field == "stage" or not measured:
            points = [
                {"x": float(pt[0]), "y": float(pt[1]), "group": s.group, "patient": s.patient_id}
                for s, pt in zip(cohort.sections, xy)
            ]
        if field == "stage" and diag["linear_fixed_effect_alias"]:
            reasons.append(field)
        warning |= diag["linear_separation"]
        warning |= diag["descriptive_auc"] >= 0.85
    joint = _joint_design(cohort.sections, labels)
    if joint["fixed_effect_alias"] and not reasons:
        reasons.append("joint_acquisition_main_effects")
    code = None
    severity = "yellow" if warning else "green"
    message = "Local bounding boxes are coordinate-frame proxies."
    if reasons:
        code = "UNIDENTIFIABLE_DESIGN"
        severity = "red"
        message = "Group is aliased with measured acquisition fixed effects. Collect a crossed design to estimate a separate group effect."
    elif flags["stage"].get("linear_separation"):
        code = "NO_STAGE_OVERLAP"
        reasons.append("stage")
        severity = "red" if len({s.patient_id for s in cohort.sections}) > 8 else "yellow"
        message = "Stage regions have disjoint support. Assess extrapolation and sample shared stage regions."
    return {
        "status": severity,
        "code": code,
        "message": message,
        "categorical": details,
        "coordinates": flags,
        "joint_acquisition_design": joint,
        "points": points,
        "points_frame": "stage"
        if all(s.stage is not None for s in cohort.sections)
        else "local_bbox",
        "reasons": reasons,
    }


def power_audit(cohort, max_exact=50000, draws=4095, alpha=0.05):
    patients = list(dict.fromkeys((s.patient_id, s.group) for s in cohort.sections))
    counts = Counter(g for _, g in patients)
    groups = sorted(counts)
    options = cohort_test_options(cohort, [p for p, _ in patients])
    blocks = exchangeability_blocks([g for _, g in patients], **options)
    result = dict(
        exact_floor(
            counts[groups[0]],
            counts[groups[1]],
            len(cohort.features),
            alpha,
            design=cohort.design,
            stratum_counts=[(len(rows) - k, k) for rows, k in blocks],
        ),
        groups=groups,
    )
    if cohort.design == "paired":
        result["n_pairs"] = len(blocks)
    elif cohort.design == "stratified":
        result["stratum"] = cohort.stratum
        result["stratum_counts"] = [
            {"level": options["strata"][rows[0]], "n0": len(rows) - k, "n1": k}
            for rows, k in blocks
        ]
    if max_exact < 1 or draws < 1:
        raise ValueError("Positive permutation budgets required")
    simulated = int(result["assignments"]) > max_exact
    if simulated:
        from fractions import Fraction

        resolution = min(Fraction(1), Fraction(2, draws + 1))
        rank = floor(resolution * len(cohort.features) / Fraction(str(result["alpha"]))) + 1
        effective_rank = rank if rank <= len(cohort.features) else None
        effective_p = float(resolution)
    else:
        effective_rank, effective_p = result["minimum_bh_rank"], result["p_floor"]
    result.update(
        effective_p_floor=effective_p,
        effective_minimum_bh_rank=effective_rank,
        effective_bh_impossible=effective_rank is None,
        executed_mode="monte_carlo_plus_one" if simulated else "exact",
        theoretical_status=result["status"],
    )
    if simulated:
        result["status"] = (
            "red" if effective_rank is None else "yellow" if effective_rank > 1 else "green"
        )
    return result
