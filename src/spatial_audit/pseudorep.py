"""Pooled rank and patient marginal-mean sensitivity diagnostics."""

from collections import defaultdict
from itertools import combinations
from math import comb

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


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


def patient_test(values, labels, max_exact=50000, draws=4095, seed=17):
    """Randomization test of a moderated signed t. Refit the prior for every allocation."""
    y = np.asarray(values, float)
    lab = np.asarray(labels)
    if y.ndim != 2 or len(y) != len(lab) or not np.isfinite(y).all() or len(set(lab)) != 2:
        raise ValueError("Finite patient-by-feature matrix and two groups required")
    n = len(lab)
    inside = lab == sorted(set(lab))[1]
    n1 = int(inside.sum())
    n0 = n - n1
    if min(n0, n1) < 2 or max_exact < 1 or draws < 1:
        raise ValueError("Invalid replicate or permutation budget")
    y = y - y.mean(axis=0)
    # Pool marginal scales without group labels. A change of units for one ion
    # must not alter the fitted cross-feature variance prior or the other tests.
    pooled_sd = np.linalg.norm(y, axis=0) / np.sqrt(n - 1)
    y /= np.where(pooled_sd > 0, pooled_sd, 1.0)
    total = np.square(y).sum(axis=0)
    df = n - 2

    def statistic(idx):
        delta = y[idx].sum(axis=0) * (1 / n0 + 1 / n1)
        v = np.maximum(total - delta**2 / (1 / n0 + 1 / n1), 0) / df
        positive = v[v > np.maximum(total, 1e-300) * 1e-12]
        if len(positive) >= 2:
            logs = np.log(positive)
            excess = max(float(logs.var(ddof=1) - polygamma(1, df / 2)), 0.0)
            prior_df = (
                1e5
                if excess <= 2.1e-5
                else 2 * brentq(lambda x: float(polygamma(1, x)) - excess, 1e-5, 1e6)
            )
            prior = np.exp(
                logs.mean()
                - digamma(df / 2)
                + np.log(df / 2)
                + digamma(prior_df / 2)
                - np.log(prior_df / 2)
            )
        else:
            prior_df = 4.0
            prior = float(np.mean(positive)) if len(positive) else 0.0
        scale = np.sqrt((df * v + prior_df * prior) / (df + prior_df) * (1 / n0 + 1 / n1))
        return delta / np.maximum(scale, np.maximum(np.sqrt(total), 1e-150) * 1e-12)

    observed = statistic(np.flatnonzero(inside))
    M = comb(n, n1)
    exact = M <= max_exact
    rng = np.random.default_rng(seed)
    permutations = (
        combinations(range(n), n1)
        if exact
        else (rng.choice(n, n1, replace=False) for _ in range(draws))
    )
    upper = np.full(y.shape[1], 0 if exact else 1)
    lower = upper.copy()
    tol = 1e-10 * np.maximum(1, abs(observed))
    for idx in permutations:
        t = statistic(list(idx))
        upper += t >= observed - tol
        lower += t <= observed + tol
    p = np.minimum(1.0, 2 * np.minimum(upper, lower) / (M if exact else draws + 1))
    return {
        "method": "exact" if exact else "monte_carlo_plus_one",
        "assignments": str(M),
        "evaluations": M if exact else draws,
        "seed": None if exact else seed,
        "variance_prior_scale": "label-blind featurewise marginal SD",
        "p": p.tolist(),
        "q": bh(p).tolist(),
        "statistic": observed.tolist(),
    }


def pseudorep_audit(cohort, max_exact=50000, draws=4095):
    sections = cohort.sections
    groups = sorted({s.group for s in sections})
    by_patient: dict = defaultdict(list)
    for s in sections:
        by_patient[s.patient_id, s.group].append(np.asarray(s.X.mean(axis=0)).ravel())
    values = [np.mean(rows, axis=0) for rows in by_patient.values()]
    test = patient_test(values, [g for _, g in by_patient], max_exact, draws)
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
    spot = float(np.mean(q < 0.05))
    patient = float(np.mean(np.asarray(test["q"]) < 0.05))
    gap = spot - patient
    return {
        "status": "red" if gap >= 0.25 else ("yellow" if gap >= 0.1 else "green"),
        "message": "Difference in BH call rates across sampling units; true discovery status is unknown for public cohorts.",
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
