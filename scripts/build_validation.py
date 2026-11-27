"""Render the validation results from frozen public reports and synthetic measurements."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

from spatial_audit.schema import strict_read

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"


def build_validation(save):
    def read(kind, name):
        return strict_read(ROOT / f"results/{kind}/{name}.json")

    vst = read("synthetic", "vst_calibration")
    scale = read("synthetic", "scaling_curve")["rows"]
    cases = read("public", "case_studies")
    panel = read("public", "panel_sensitivity")
    null = read("public", "empirical_null")
    ad, tnbc, msi = [
        read("public", k)["dimensions"] for k in ("GSE220442", "GSE348499", "MTBLS12782")
    ]
    rows = vst["rows"]
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.6))
    for ax, state in zip(axes[0], ("null", "alternative"), strict=True):
        for mode, color in (("raw", "#235f65"), ("vst", "#bd7048")):
            sub = sorted(
                (
                    r
                    for r in rows
                    if r["state"] == state
                    and r["nugget_mode"] == mode
                    and r["case_mean_multiplier"] == 2
                    and r["n_per_arm"] == 7
                ),
                key=lambda r: r["nb_dispersion"],
            )
            rates = np.array([r["rate"] for r in sub])
            bounds = np.array([r["wilson95"] for r in sub])
            ax.errorbar(
                [r["nb_dispersion"] for r in sub],
                rates,
                yerr=np.maximum(0, np.array([rates - bounds[:, 0], bounds[:, 1] - rates])),
                color=color,
                marker="o",
                label=mode.upper(),
                capsize=3,
            )
        if state == "null":
            ax.axhline(0.05, ls=":", color="#555555")
        ax.set(
            title="Latent-amplitude " + state,
            xlabel="Negative-binomial dispersion",
            ylabel="Rejection fraction",
            ylim=(0, 1.03),
        )
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.legend(fontsize=8)
    for ax, key, label, divisor in zip(
        axes[1],
        ("wall_seconds", "peak_rss_bytes"),
        ("Audit wall time (seconds)", "Peak RSS (MiB)"),
        (1, 2**20),
        strict=True,
    ):
        x = [r["n_observations"] for r in scale]
        ax.plot(x, [r[key] / divisor for r in scale], color="#235f65", marker="o")
        ax.set(xscale="log", xlabel="Observations", ylabel=label, ylim=(0, None))
        ax.set_xticks(x, [f"{n // 1000}k" if n < 1e6 else f"{n / 1e6:g}M" for n in x])
    save(fig, "vst-scaling")

    def feature(d, name):
        p = d["pseudorep"]
        j = p["features"].index(name)
        return p["spot_q"][j], p["patient"]["q"][j], d["nugget"]["features"][j]

    pooled_ad, patient_ad, _ = feature(ad, "SDF4")
    _, _, abc = feature(tnbc, "ABCB4")
    _, _, ion = feature(msi, "mz_178.00810")
    pooled_abc, patient_abc, _ = feature(tnbc, "ABAT")
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8))
    for ax, values, title, unit, labels in zip(
        axes,
        (
            [pooled_ad, patient_ad],
            [abc["raw_delta"] * 1e5, abc["debiased_delta"] * 1e5],
            [ion["raw_delta"], ion["debiased_delta"]],
        ),
        ("Cortex · SDF4", "TNBC · ABCB4", "PDAC · m/z 178.00810"),
        ("BH-adjusted p", "Energy difference (× 10⁻⁵)", "Energy difference"),
        (["Pooled", "Donor"], ["Before", "After"], ["Before", "After"]),
        strict=True,
    ):
        ax.set(title=title, ylabel=unit)
        if ax == axes[0]:
            ax.scatter(labels, values, color=["#bd7048", "#235f65"], s=40)
            ax.set_yscale("log")
            ax.axhline(0.05, color="#888888", ls=":")
            ax.set_ylim(min(values) * 0.1, 3)
        else:
            ax.bar(labels, values, color=["#bd7048", "#235f65"])
            ax.axhline(0, color="#555555", lw=0.6)
    save(fig, "case-studies")

    small, large = panel["panels"]
    worst = {
        m: max(
            (r for r in rows if r["state"] == "null" and r["nugget_mode"] == m),
            key=lambda r: r["rate"],
        )
        for m in ("raw", "vst")
    }

    def pair(state):
        return [
            r["rate"]
            for r in rows
            if r["state"] == state
            and r["n_per_arm"] == 7
            and r["case_mean_multiplier"] == 2
            and r["nb_dispersion"] == 0.35
        ]

    raw, vst_alt = pair("alternative")
    fibrosis = cases["fibrosis"]
    baseline, loo = cases["msi"]["records"]
    paragraphs = [
        (
            f"In {len(rows) // 4} simulation settings with no batch alignment, each mode used the same count "
            f"realizations ({vst['replicates_per_state']} cohorts per effect state and setting). Across settings, "
            f"the largest null rejection fractions were {worst['raw']['rate']:.1%} on the raw scale and "
            f"{worst['vst']['rate']:.1%} after Anscombe transformation. The latter was "
            f"{worst['vst']['rejections']}/{vst['replicates_per_state']} cohorts (pointwise 95% Wilson "
            f"interval {worst['vst']['wilson95'][0]:.1%} to {worst['vst']['wilson95'][1]:.1%}) under Poisson "
            "sampling with seven donors per arm and a twofold case baseline. At the same donor count and mean "
            f"shift with dispersion 0.35, alternative detection was {raw:.1%} and {vst_alt:.1%}, respectively. "
            "A fixed additive pattern in counts changes amplitude after a square-root transformation as baseline "
            "abundance changes. These rates evaluate the original latent-amplitude null."
        ),
        f"GSE348499 supplied {null['splits_per_arm']} random donor splits within each pretreatment breast cancer "
        "response group. Across splits and features, the BH call fractions for pooled cells, donor means, "
        "raw-scale corrected energy and VST corrected energy, respectively, were "
        + " and ".join(
            ", ".join(
                f"{a['summary'][k]['bh_fraction']:.2%}"
                for k in ("pooled", "patient_mean", "raw_corrected", "vst_corrected")
            )
            + f" within {a['group']} ({len(a['donors'])} donors)"
            for a in null["arms"]
        )
        + ". Splits reuse donor material. The fractions describe random-label sensitivity conditional on these donors.",
        (
            f"In GSE220442, expanding the same {panel['n_patients']} cortical donors from {small['features']} to "
            f"{large['features']:,} archived features preserved all five screening categories. Pooled BH calls "
            f"rose from {small['pooled_rate']:.1%} to {large['pooled_rate']:.1%}. Both panels had zero donor-level "
            f"calls. Cortical-layer covariates accounted for {small['bands'][-1]['share']:.1%} and "
            f"{large['bands'][-1]['share']:.1%} of energy at 1,400 micrometres. Neither panel supported the "
            "250-micrometre nugget endpoint."
        ),
        (
            f"For cortical SDF4, the pooled adjusted p was {pooled_ad:.2g} and the donor-adjusted p was "
            f"{patient_ad:.2g}. Three donors per arm limit the smallest two-sided donor p to 0.10. In the "
            f"breast cancer panel, ABAT yielded pooled and donor-adjusted p values of {pooled_abc:.2g} and "
            f"{patient_abc:.3f}. For ABCB4, the pCR-minus-RD energy difference changed from "
            f"{abc['raw_delta']:.3g} to {abc['debiased_delta']:.3g} after short-lag subtraction. "
            "The ABCB4 contrast compares pretreatment spatial variation between eventual response groups. "
            "Without measured cell-type covariates in these sections, the contribution of cellular "
            "composition to this difference remains unresolved."
        ),
        (
            f"In pancreatic cancer, the without-minus-with-diabetes energy difference for m/z 178.00810 "
            f"changed from {ion['raw_delta']:.3f} to {ion['debiased_delta']:.3f}. The archive identifies "
            "this feature by mass window. Across the panel, TIC and raster coordinates together accounted for "
            f"{baseline['shares'][-1]:.2%} of energy at 1,400 micrometres, versus {loo['shares'][-1]:.2%} "
            "after subtracting each target window from its own TIC covariate. The resulting share still "
            "includes spatial gradients shared by tissue and acquisition."
        ),
        (
            f"The {fibrosis['n_patients']}-donor lung fibrosis Xenium contrast GSE250346 yielded "
            f"{len(fibrosis['called']['patient'])}/{fibrosis['n_features']} donor-level BH calls, "
            f"{len(fibrosis['called']['edgeR'])} edgeR calls and {len(fibrosis['called']['DESeq2'])} "
            f"DESeq2 calls. All three procedures called {len(fibrosis['shared_calls'])} shared features "
            f"({', '.join(fibrosis['shared_calls'])}). AGER and CA4 link the panel to alveolar and vascular "
            "compartments, respectively. Whole-tissue summaries combine changes in compartment "
            "abundance and expression within cells. All three analyses share donors and assay data."
        ),
        (
            f"Across {len(scale)} independently launched CPU audits with 10 donors and 64 sparse features, "
            f"wall time increased from {scale[0]['wall_seconds']:.1f} to {scale[-1]['wall_seconds']:.1f} "
            f"seconds and peak RSS from {scale[0]['peak_rss_bytes'] / 2**20:.1f} to "
            f"{scale[-1]['peak_rss_bytes'] / 2**20:.1f} MiB between {scale[0]['n_observations']:,} and "
            f"{scale[-1]['n_observations']:,} observations. Each size has one host-specific measurement."
        ),
    ]

    def section(path, blocks):
        start, end = "<!-- BEGIN measured validation -->", "<!-- END measured validation -->"
        before, rest = path.read_text(encoding="utf-8").split(start)
        _, after = rest.split(end)
        path.write_text(
            before + start + "\n\n" + "\n\n".join(blocks) + "\n\n" + end + after, encoding="utf-8"
        )

    section(ROOT / "README.md", paragraphs)
    section(PAPER / "paper.md", paragraphs)
    tex = "\n\n".join(paragraphs).replace("%", r"\%")
    tex = tex.replace("GSE250346 yielded", r"GSE250346 \cite{geo250346} yielded")
    # Set scientific notation as math so exponents cannot wrap away from their values.
    import re

    tex = re.sub(
        r"(?<![A-Za-z])(-?\d+(?:\.\d+)?)e([+-]\d+)",
        lambda m: "$" + m[1] + r"\times 10^{" + str(int(m[2])) + "}$",
        tex,
    )
    (PAPER / "generated/validation.tex").write_text(tex + "\n", encoding="utf-8")

    table = [
        r"\begin{center}\small\begin{tabular}{llrrr}\toprule",
        r"Arm & Endpoint & Nominal (\%) & BH calls/tests & Splits with calls\\\midrule",
    ]
    endpoints = {
        "pooled": "Pooled cells",
        "patient_mean": "Donor mean",
        "raw_corrected": "Raw-scale energy",
        "vst_corrected": "VST energy",
    }
    for arm in null["arms"]:
        for key, label in endpoints.items():
            data = arm["summary"][key]
            table.append(
                f"{arm['group']} & {label} & {100 * data['nominal_fraction']:.2f} & "
                f"{data['bh_calls']}/{data['tests']} & "
                f"{data['splits_with_calls']}/{null['splits_per_arm']}" + r"\\"
            )
    table.extend(
        [
            r"\bottomrule\end{tabular}\end{center}",
            r"\begin{center}\small\begin{tabular}{rrr}\toprule",
            r"Observations & Audit seconds & Peak RSS (MiB)\\\midrule",
        ]
    )
    for point in scale:
        table.append(
            f"{point['n_observations']:,} & {point['wall_seconds']:.2f} & "
            f"{point['peak_rss_bytes'] / 2**20:.1f}" + r"\\"
        )
    table.append(r"\bottomrule\end{tabular}\end{center}")
    (PAPER / "generated/validation-tables.tex").write_text(
        "\n".join(table) + "\n", encoding="utf-8"
    )
