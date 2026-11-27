"""Compute auditable records shared by CLI, dashboard, benchmark and MCP."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from .design_audit import design_audit, power_audit
from .permutation import cohort_test_options
from .pseudorep import patient_test, pseudorep_audit
from .schema import clean, write_json
from .spectral_audit import spectral_audit

DIMENSIONS = ("power", "design", "pseudorep", "composition", "nugget")
STATEMENT = (
    "基于公开存档数据，在患者/生物学重复级去混杂框架下的可复现性、统计功效与混杂敏感性审计。"
)


def describe_noise_scale(report):
    """Attach the scope of a sensitivity screen on the input's measurement scale."""
    provenance = report["provenance"]
    nugget = report["dimensions"]["nugget"]
    nugget["measurement_scope"] = {
        "measurement": provenance.get("measurement"),
        "layer": provenance.get("layer"),
        "interpretation": (
            "Energy, semivariance and attenuation refer to the supplied measurement scale. "
            "For log-normalized inputs they describe variability after transformation. "
            "The procedure does not fit a negative-binomial count model or identify a noise cause."
        ),
    }
    if nugget.get("nugget_mode") == "vst":
        nugget["measurement_scope"]["interpretation"] = (
            "Fine-band residual energies and noise estimates use the Anscombe transform of "
            "nonnegative volume means. Composition shares use the supplied measurement scale."
        )
    return report


def audit(
    cohort,
    wavelengths=(250.0, 600.0, 1400.0),
    max_exact=50000,
    draws=4095,
    *,
    design=None,
    stratum=None,
    alpha=0.05,
    nugget_mode="raw",
    graph_mode="grid",
    leave_one_out_tic=False,
):
    cohort = replace(cohort, design=design or cohort.design, stratum=stratum or cohort.stratum)
    cohort.validate()
    power = power_audit(cohort, max_exact, draws, alpha)
    composition, nugget = spectral_audit(
        cohort,
        wavelengths=wavelengths,
        nugget_mode=nugget_mode,
        graph_mode=graph_mode,
        leave_one_out_tic=leave_one_out_tic,
    )
    energies = nugget.get("patient_energies", [])
    for endpoint in ("raw", "corrected"):
        if not energies:
            continue
        values = np.asarray([row[endpoint] for row in energies], float)
        valid = np.isfinite(values).all(axis=0)
        if valid.any():
            nugget[endpoint + "_test"] = dict(
                features=np.asarray(cohort.features)[valid].tolist(),
                **patient_test(
                    values[:, valid],
                    [r["group"] for r in energies],
                    max_exact,
                    draws,
                    **cohort_test_options(cohort, [r["patient_id"] for r in energies]),
                ),
            )
    dimensions = {
        "power": power,
        "design": design_audit(cohort),
        "pseudorep": pseudorep_audit(cohort, max_exact, draws, alpha),
        "composition": composition,
        "nugget": nugget,
    }
    return describe_noise_scale(
        clean(
            {
                "schema_version": 1,
                "id": cohort.id,
                "kind": cohort.kind,
                "platform": cohort.platform,
                "source": cohort.source,
                "accession": cohort.provenance.get("accession"),
                "provenance": cohort.provenance,
                **(
                    {"preprocessing_evidence": cohort.provenance["preprocessing_evidence"]}
                    if "preprocessing_evidence" in cohort.provenance
                    else {}
                ),
                "n_patients": len({s.patient_id for s in cohort.sections}),
                "n_sections": len(cohort.sections),
                "n_observations": sum(s.X.shape[0] for s in cohort.sections),
                "n_features": len(cohort.features),
                "dimensions": dimensions,
                "statement": STATEMENT
                if cohort.kind == "public"
                else "Local spatial cohort"
                if cohort.kind == "local"
                else "Synthetic Benchmark Fixture · 合成模拟夹具",
            }
        )
    )


def save_report(report, out):
    folder = Path(out) / report["kind"]
    folder.mkdir(parents=True, exist_ok=True)
    name = report["id"]
    write_json(folder / (name + ".json"), report)
    d = report["dimensions"]
    lines = [
        f"# {name}",
        "",
        report["statement"],
        "",
        (
            f"Biological units {report['n_patients']} · "
            f"sections {report['n_sections']} · observations {report['n_observations']}"
        ),
        "",
        "| Dimension | Status | Interpretation |",
        "|---|---|---|",
    ]
    lines += [f"| {key} | {d[key]['status']} | {d[key]['message']} |" for key in DIMENSIONS]
    lines += [
        "",
        (
            f"Two-sided floor {d['power']['p_floor']:.6g}. "
            f"Spot BH calls {d['pseudorep']['spot_rate']:.1%}, "
            f"patient BH calls {d['pseudorep']['patient_rate']:.1%}."
        ),
        "",
        f"Permutation design: {d['power'].get('design', 'independent')}. BH threshold: {d['power']['alpha']:g}.",
        d["nugget"].get("measurement_scope", {}).get("interpretation", ""),
        "",
        report["source"],
    ]

    def cell(value):
        return str(value).replace("|", r"\|").replace("\n", " ")

    def number(value):
        return "Unassessed" if value is None else f"{value:.4g}"

    power, pseudo, noise = d["power"], d["pseudorep"], d["nugget"]
    lines += ["", "## Group counts", "", "| Group | Biological units |", "|---|---:|"]
    lines += [f"| {cell(g)} | {n} |" for g, n in zip(power["groups"], [power["n0"], power["n1"]])]
    lines += [
        "",
        "## Physical wavelength shares",
        "",
        "| Wavelength (µm) | Covariate share | Supported sections |",
        "|---:|---:|---:|",
    ]
    for band in d["composition"]["bands"]:
        lines.append(
            f"| {number(band['wavelength_um'])} | {number(band['share'])} | {band.get('supported_sections', 'Unassessed')} |"
        )
    lines += [
        "",
        "## Acquisition design flags",
        "",
        ", ".join(map(cell, d["design"]["reasons"])) or "No acquisition design flags recorded.",
    ]
    rows = sorted(
        zip(pseudo["features"], pseudo["spot_q"], pseudo["patient"]["q"], strict=True),
        key=lambda row: row[2] - row[1],
        reverse=True,
    )
    lines += [
        "",
        "## Top 10 pseudoreplication gaps",
        "",
        "Ranked by patient q minus pooled q.",
        "",
        "| Feature | Pooled q | Patient q | Difference |",
        "|---|---:|---:|---:|",
    ]
    lines += [
        f"| {cell(f)} | {number(q)} | {number(p)} | {number(p - q)} |" for f, q, p in rows[:10]
    ]
    flagged = sorted(
        (r for r in noise["features"] if r["mean_noise_sensitive"]),
        key=lambda r: abs(r["raw_delta"]) - abs(r["debiased_delta"]),
        reverse=True,
    )
    lines += [
        "",
        "## Top 10 nugget-sensitive features",
        "",
        f"Calibration mode: {noise.get('nugget_mode', 'raw')}. Ranked by absolute contrast attenuation.",
        "",
    ]
    if flagged:
        lines += ["| Feature | Raw energy contrast | Debiased energy contrast |", "|---|---:|---:|"]
        lines += [
            f"| {cell(r['feature'])} | {number(r['raw_delta'])} | {number(r['debiased_delta'])} |"
            for r in flagged[:10]
        ]
    else:
        lines.append(f"No features flagged among {noise['evaluated']} evaluable features.")
    (folder / (name + ".md")).write_text("\n".join(lines) + "\n", encoding="utf-8")
