"""Deterministic synthetic positive and negative controls with specified mechanisms."""

import numpy as np

from .schema import Cohort, Section


def trap_cohort(trap, seed=36):
    if trap not in {"null", "pseudorep", "composition", "stage"}:
        raise ValueError("unknown trap")
    y, x = np.mgrid[:28, :28]
    coords = np.column_stack((x.ravel(), y.ravel())) * 55.0
    groups = ("A", "B")
    sections = []
    offsets = ([-5, -4, -2, 0, 2, 3], [-3, -2, 0, 2, 4, 5])
    for i in range(6):
        rng = np.random.default_rng([seed, i])
        component = (
            0.5
            + 0.25 * np.cos((x.ravel() + 0.5) * np.pi / 8)
            + 0.20 * np.cos((y.ravel() + 0.5) * np.pi / 4)
        )
        smooth = np.cos((y.ravel()[:, None] + 0.5) * np.pi / 8) * np.linspace(1.2, 3.0, 8)
        # Common random numbers make the null an exact negative control for
        # group differences. These fixtures measure behavior on four named
        # constructions; random-null calibration is evaluated separately.
        background = rng.normal(0, 0.4, size=(len(coords), 8))
        for g, label in enumerate(groups):
            matrix = 45.0 + i / 3 + smooth + background
            c = np.cos((x.ravel() + 0.5) * np.pi / 8)[:, None]
            if trap == "pseudorep":
                matrix = matrix + offsets[g][i]
            if trap == "composition":
                level = (25.0, 100.0)[g]
                c = component[:, None]
                # Linear composition plus a group-dependent count mean; the
                # shot noise has conditional variance equal to the local mean.
                expectation = level + (component[:, None] * (22.0 + g * 36.0)) * np.linspace(
                    0.8, 1.2, 8
                )
                matrix = rng.poisson(expectation).astype(float)
            stage = (float(i % 3) * 600.0, float(i // 3) * 600.0)
            if trap == "stage":
                stage = (stage[0] + g * 4000, stage[1])
            patient = f"{label}{i + 1}"
            sections.append(
                Section(
                    patient,
                    label,
                    coords,
                    matrix,
                    C=c,
                    slide_id=f"s{i % 2}",
                    batch=f"b{i % 3}",
                    stage=stage,
                    section_id=patient,
                    domain=np.where(x.ravel() < 14, "left", "right"),
                )
            )
    return Cohort(
        f"synthetic-{trap}",
        "synthetic",
        [f"F{i + 1}" for i in range(8)],
        sections,
        platform="Synthetic lattice",
        source="Synthetic Benchmark Fixtures",
        provenance={
            "measurement": (
                "synthetic Poisson counts"
                if trap == "composition"
                else "synthetic continuous intensities"
            ),
            "selection": "All eight generated features and twelve constructed independent-unit labels",
            "fixture_seed": int(seed),
            "construction": {
                "null": "Matched fields and common noise in each group",
                "pseudorep": "Matched fields with patient-specific constant offsets",
                "composition": "Smooth covariate and group-specific conditional Poisson means",
                "stage": "Matched fields with a group-specific acquisition-coordinate offset",
            }[trap],
        },
    ).validate()


def fixtures():
    return [trap_cohort(name) for name in ("null", "pseudorep", "composition", "stage")]
