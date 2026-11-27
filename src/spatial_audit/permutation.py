"""Patient-level exchangeability blocks and allocation counts."""

from itertools import combinations, product
from math import comb, prod

import numpy as np


def exchangeability_blocks(labels, design, patient_ids=None, strata=None):
    inside = np.asarray(labels) == sorted(set(labels))[1]
    if design == "independent":
        return [(np.arange(len(labels)), int(inside.sum()))]
    if design not in {"paired", "stratified"}:
        raise ValueError("Unknown permutation design")
    keys = patient_ids if design == "paired" else strata
    if keys is None or len(keys) != len(labels) or any(k is None for k in keys):
        raise ValueError("Permutation block labels must align with patient summaries")
    groups: dict = {}
    for i, key in enumerate(keys):
        groups.setdefault(key, []).append(i)
    blocks = [(np.asarray(rows), int(inside[rows].sum())) for rows in groups.values()]
    if design == "paired" and any(len(rows) != 2 or count != 1 for rows, count in blocks):
        raise ValueError("Each pair requires one summary from each condition")
    return blocks


def assignment_count(blocks):
    return prod(comb(len(rows), count) for rows, count in blocks)


def allocations(blocks, exact, draws, rng):
    if exact:
        pools = [combinations(rows, count) for rows, count in blocks]
        if len(pools) == 1:
            yield from pools[0]
        else:
            for selected in product(*pools):
                yield tuple(i for part in selected for i in part)
    else:
        for _ in range(draws):
            yield np.concatenate([rng.choice(rows, count, replace=False) for rows, count in blocks])


def cohort_test_options(cohort, patient_ids):
    levels = cohort.patient_strata() if cohort.design == "stratified" else {}
    return {
        "design": cohort.design,
        "patient_ids": patient_ids,
        "strata": [levels.get(pid) for pid in patient_ids],
    }
