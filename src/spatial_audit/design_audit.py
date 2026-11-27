"""Design diagnostics for independent biological units."""

from collections import Counter, defaultdict
from math import comb, floor, log2

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.linalg import lsqr
from scipy.stats import fisher_exact, rankdata


def exact_floor(n0: int, n1: int, m: int = 1, alpha: float = 0.05) -> dict:
    if any(type(i) is not int or i < 1 for i in (n0, n1, m)) or not 0 < alpha < 1:
        raise ValueError("Positive integer counts and 0 < alpha < 1 required")
    M = comb(n0 + n1, n0)
    p = min(1.0, 2 / M)
    # Strict BH q<alpha; use rational arithmetic at exact boundaries.
    from fractions import Fraction

    rank = floor(min(Fraction(1), Fraction(2, M)) * m / Fraction(str(alpha))) + 1
    status = "red" if rank > m else ("yellow" if rank > 1 else "green")
    return {
        "status": status,
        "message": "Two-sided doubled-tail permutation resolution.",
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
            details[field] = {"status": "missing_or_partial"}
            warning = True
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
        if t.shape == (2, 2) and all(len(v) == 1 for v in visits.values()):
            fisher = float(fisher_exact(t.astype(int)).pvalue)
        details[field] = {
            "status": "observed",
            "label": cohort.provenance.get("batch_label", "acquisition batch") if field == "batch" else "slide",
            "levels": len(t),
            "normalized_mi": mi / entropy,
            "complete_alias": separate,
            "fisher_p": fisher,
            "all_levels_patient_specific": max(sizes.values()) == 1,
        }
        if separate:
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
        flags[field] = diag
        if field == "stage" or not measured:
            points = [
                {"x": float(pt[0]), "y": float(pt[1]), "group": s.group, "patient": s.patient_id}
                for s, pt in zip(cohort.sections, xy)
            ]
        if diag["linear_separation"]:
            small_sample = len({s.patient_id for s in cohort.sections}) <= 8
            if field == "stage" and (diag["linear_fixed_effect_alias"] or not small_sample):
                reasons.append(field)
            else:
                warning = True
        warning |= diag["descriptive_auc"] >= 0.85
    joint = _joint_design(cohort.sections, labels)
    if joint["fixed_effect_alias"] and not reasons:
        reasons.append("joint_acquisition_main_effects")
    return {
        "status": "red" if reasons else ("yellow" if warning else "green"),
        "code": "UNIDENTIFIABLE_DESIGN" if reasons else None,
        "message": "High risk marks categorical aliasing or no observed stage overlap. The linear fixed-effect rank diagnostic is reported separately. Local bounding boxes are coordinate-frame proxies.",
        "categorical": details,
        "coordinates": flags,
        "joint_acquisition_design": joint,
        "points": points,
        "points_frame": "stage"
        if all(s.stage is not None for s in cohort.sections)
        else "local_bbox",
        "reasons": reasons,
    }


def power_audit(cohort, max_exact=50000, draws=4095):
    patients = {s.patient_id: s.group for s in cohort.sections}
    counts = Counter(patients.values())
    groups = sorted(counts)
    result = dict(
        exact_floor(counts[groups[0]], counts[groups[1]], len(cohort.features)), groups=groups
    )
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
