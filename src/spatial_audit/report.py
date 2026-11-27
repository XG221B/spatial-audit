"""Compute auditable records shared by CLI, dashboard, benchmark and MCP."""

from pathlib import Path

import numpy as np

from .design_audit import design_audit, power_audit
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
    return report


def audit(cohort, wavelengths=(250.0, 600.0, 1400.0), max_exact=50000, draws=4095):
    cohort.validate()
    composition, nugget = spectral_audit(cohort, wavelengths=wavelengths)
    energies = nugget.get("patient_energies", [])
    for endpoint in ("raw", "corrected"):
        if not energies:
            continue
        values = np.asarray([row[endpoint] for row in energies], float)
        valid = np.isfinite(values).all(axis=0)
        if valid.any():
            nugget[endpoint + "_test"] = dict(
                features=np.asarray(cohort.features)[valid].tolist(),
                **patient_test(values[:, valid], [r["group"] for r in energies], max_exact, draws),
            )
    dimensions = {
        "power": power_audit(cohort, max_exact, draws),
        "design": design_audit(cohort),
        "pseudorep": pseudorep_audit(cohort, max_exact, draws),
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
        "Group tests assume exchangeable patients. Inspect acquisition support before interpretation.",
        d["nugget"].get("measurement_scope", {}).get("interpretation", ""),
        "",
        report["source"],
    ]
    (folder / (name + ".md")).write_text("\n".join(lines) + "\n", encoding="utf-8")
