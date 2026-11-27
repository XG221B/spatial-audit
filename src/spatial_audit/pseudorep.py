"""Pooled rank and patient marginal-mean sensitivity diagnostics."""

from collections import defaultdict
from itertools import batched

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .h5_backed import value_blocks
from .permutation import allocations, assignment_count, cohort_test_options, exchangeability_blocks


def bh(p):
    values = np.asarray(p, dtype=np.float64)
    if not np.isfinite(values).all() or np.any((values < 0) | (values > 1)):
        raise ValueError("A finite p value in [0, 1] is required for every registered hypothesis")
    return multipletests(values, method="fdr_bh")[1]


def pooled_rank_p(group_a, group_b):
    a, b = (np.asarray(x, dtype=np.float64) for x in (group_a, group_b))
    if a.min() == a.max() == b.min() == b.max():
        return 1.0
    value = float(mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").pvalue)
    if not np.isfinite(value):
        raise ValueError("Nonfinite rank-test probability for a nonconstant feature")
    return value


def patient_test(
    values,
    labels,
    max_exact=50000,
    draws=4095,
    seed=17,
    *,
    design="independent",
    patient_ids=None,
    strata=None,
):
    """Use joint sign flips of paired differences or a moderated t refitted per allocation."""
    y = np.array(values, dtype=float, copy=True)
    lab = np.asarray(labels)
    if y.ndim != 2 or len(y) != len(lab) or not np.isfinite(y).all() or len(set(lab)) != 2:
        raise ValueError("Finite patient-by-feature matrix and two groups required")
    n = len(lab)
    inside = lab == sorted(set(lab))[1]
    n1 = int(inside.sum())
    n0 = n - n1
    if min(n0, n1) < 2 or max_exact < 1 or draws < 1:
        raise ValueError("Invalid replicate or permutation budget")
    blocks = exchangeability_blocks(lab, design, patient_ids, strata)
    if design == "paired":
        for pair, _ in blocks:
            y[pair] -= y[pair].mean(axis=0)
    y = y - y.mean(axis=0)
    pooled_sd = np.linalg.norm(y, axis=0) / np.sqrt(n - 1)
    y /= np.where(pooled_sd > 0, pooled_sd, 1.0)
    total = np.square(y).sum(axis=0)
    df = n - 2
    factor = 1 / n0 + 1 / n1
    log_df_offset = np.log(df / 2) - digamma(df / 2)
    sampling_variance = polygamma(1, df / 2)

    def statistics(indices):
        # Preserve sequential summation within an allocation; batch across allocations.
        delta = np.sum(y[indices], axis=1) * factor
        if design == "paired":
            return delta
        v = np.maximum(total - delta**2 / factor, 0) / df
        positive = v > np.maximum(total, 1e-300) * 1e-12
        counts = positive.sum(axis=1)
        logs = np.log(np.where(positive, v, 1.0))
        means = logs.sum(axis=1) / np.maximum(1, counts)
        centered = np.where(positive, logs - means[:, None], 0)
        excess = (centered**2).sum(axis=1) / np.maximum(1, counts - 1) - sampling_variance
        prior_df = np.full(len(v), 1e5)
        for row in np.flatnonzero((counts >= 2) & (excess > 2.1e-5)):
            prior_df[row] = 2 * brentq(
                lambda x, target: polygamma(1, x) - target, 1e-5, 1e6, args=(excess[row],)
            )
        prior = np.exp(means + log_df_offset + digamma(prior_df / 2) - np.log(prior_df / 2))
        fallback = counts < 2
        prior_df[fallback] = 4
        prior[fallback] = np.where(positive, v, 0).sum(axis=1)[fallback]
        scale = np.sqrt((df * v + (prior_df * prior)[:, None]) / (df + prior_df[:, None]) * factor)
        return delta / np.maximum(scale, np.maximum(np.sqrt(total), 1e-150) * 1e-12)

    observed = statistics(np.flatnonzero(inside)[None, :])[0]
    M = assignment_count(blocks)
    exact = M <= max_exact
    permutations = allocations(blocks, exact, draws, np.random.default_rng(seed))
    upper = np.full(y.shape[1], 0 if exact else 1)
    lower = upper.copy()
    tol = 1e-10 * np.maximum(1, abs(observed))
    batch_size = max(1, min(128, 1_048_576 // max(1, n1 * y.shape[1])))
    for batch in batched(permutations, batch_size):
        t = statistics(np.asarray(batch))
        upper += (t >= observed - tol).sum(axis=0)
        lower += (t <= observed + tol).sum(axis=0)
    p = np.minimum(1.0, 2 * np.minimum(upper, lower) / (M if exact else draws + 1))
    return {
        "method": "exact" if exact else "monte_carlo_plus_one",
        "assignments": str(M),
        "evaluations": M if exact else draws,
        "seed": None if exact else seed,
        "design": design,
        "statistic_type": "standardized paired mean difference"
        if design == "paired"
        else "moderated t",
        "variance_prior_scale": None
        if design == "paired"
        else "label-blind featurewise marginal SD",
        "p": p.tolist(),
        "q": bh(p).tolist(),
        "statistic": observed.tolist(),
    }


def pseudorep_audit(cohort, max_exact=50000, draws=4095, alpha=0.05):
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie between zero and one")
    sections = cohort.sections
    groups = sorted({s.group for s in sections})
    by_patient: dict = defaultdict(list)
    for s in sections:
        by_patient[s.patient_id, s.group].append(np.asarray(s.X.mean(axis=0)).ravel())
    values = [np.mean(rows, axis=0) for rows in by_patient.values()]
    labels = [g for _, g in by_patient]
    options = cohort_test_options(cohort, [pid for pid, _ in by_patient])
    test = patient_test(values, labels, max_exact, draws, **options)
    depth = depth_diagnostic(cohort, by_patient, max_exact, draws, alpha, options)
    p = []
    for j in range(len(cohort.features)):
        samples = []
        for g in groups:
            col = [s.X[:, j : j + 1] for s in sections if s.group == g]
            samples.append(
                np.concatenate([(x.toarray() if sparse.issparse(x) else x).ravel() for x in col])
            )
        p.append(pooled_rank_p(*samples))
    q = bh(p)
    spot = float(np.mean(q < alpha))
    patient = float(np.mean(np.asarray(test["q"]) < alpha))
    gap = spot - patient
    return {
        "status": "red" if gap >= 0.25 else ("yellow" if gap >= 0.1 else "green"),
        "message": "Difference in BH call rates across sampling units; true discovery status is unknown for public cohorts.",
        **depth,
        "alpha": alpha,
        "spot_rate": spot,
        "patient_rate": patient,
        "gap": gap,
        "ratio": spot / patient if patient else None,
        "zero_denominator": patient == 0,
        "groups": groups,
        "features": cohort.features,
        "spot_p": p,
        "spot_q": q.tolist(),
        "patient": test,
    }


def depth_diagnostic(cohort, by_patient, max_exact, draws, alpha, options):
    """Divide section mean abundance by mean exposure, then equally weight sections per donor."""
    measurement = (cohort.provenance.get("measurement") or "").strip().lower()
    count_scales = {
        "count",
        "counts",
        "raw",
        "raw_counts",
        "raw counts",
        "synthetic poisson counts",
    }
    reason = None
    if measurement and measurement not in count_scales:
        reason = f"unavailable for measurement scale: {measurement}"
    else:
        for section in cohort.sections:
            for block in value_blocks(section.X):
                if np.any(block < 0):
                    reason = "unavailable for signed measurements"
                    break
                if not measurement and np.any(block != np.floor(block)):
                    reason = "unspecified measurement scale with noninteger values"
                    break
            if reason:
                break
    if reason:
        return {
            "depth_normalization_applicable": False,
            "depth_normalized_patient_rate": None,
            "group_depth_ratio": None,
            "depth_source": reason,
        }
    external = [s.library_size is not None for s in cohort.sections]
    source = "supplied library size" if all(external) else "registered panel sum"
    normalized: dict = defaultdict(list)
    depths: dict = defaultdict(list)
    for s in cohort.sections:
        key = (s.patient_id, s.group)
        means = by_patient[key][len(depths[key])]
        depth = float(np.mean(s.library_size)) if all(external) else float(means.sum())
        depths[key].append(depth)
        if depth > 0:
            normalized[key].append(means / depth)
    labels = np.asarray([g for _, g in by_patient])
    group_depths = np.array([np.mean(depths[key]) for key in by_patient])
    a, b = (group_depths[labels == g].mean() for g in sorted(set(labels)))
    result = {
        "depth_normalization_applicable": True,
        "depth_source": source,
        "group_depth_ratio": float(b / a) if a else None,
        "depth_normalized_patient_rate": None,
    }
    if all(len(normalized[key]) == len(rows) for key, rows in by_patient.items()):
        values = [np.mean(normalized[key], axis=0) for key in by_patient]
        test = patient_test(values, labels, max_exact, draws, **options)
        result["depth_normalized_patient_rate"] = float(np.mean(np.asarray(test["q"]) < alpha))
    return result
