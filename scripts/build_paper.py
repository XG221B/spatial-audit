"""Derive manuscript statistics and vector figures exclusively from result JSON."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spatial_audit.schema import strict_read

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PAPER = ROOT / "paper"
INK = "#235f65"
ACCENT = "#bd7048"


def save(fig, stem):
    fig.savefig(PAPER / "figures" / (stem + ".pdf"), dpi=300, metadata={"Creator": "spatial-audit"})
    plt.close(fig)


def audits(kind):
    values = [strict_read(p) for p in sorted((RESULTS / kind).glob("*.json"))]
    return [
        r for r in values if isinstance(r.get("dimensions"), dict) and r.get("schema_version") == 1
    ]


def refinement_table():
    """Render before/after measurements in the review from versioned scientific receipts."""
    rank = strict_read(RESULTS / "public/STT0000036_refinement.json")
    noise = strict_read(RESULTS / "synthetic/nugget_refinement.json")
    lattice = strict_read(RESULTS / "public/GSE220442_reconstruction.json")["sections"]
    doc = ROOT / "docs/SELF_CRITIQUE_AND_REFINEMENT.md"
    opening, closing = "<!-- BEGIN measured comparisons -->", "<!-- END measured comparisons -->"
    before, rest = doc.read_text(encoding="utf-8").split(opening)
    _, after = rest.split(closing)

    def rates(record, key):
        return " | ".join(f"{record[stage][key]:.2%}" for stage in ("before", "after"))

    pitch = [r["legacy_image_derived_pitch_um"] for r in lattice]
    nominal = {r["nominal_array_pitch_um"] for r in lattice}
    if len(nominal) != 1:
        raise ValueError("Expected a common GSE220442 lattice pitch")
    rows = [
        "| 实测指标 | 修复前 | 修复后 |",
        "|---|---:|---:|",
        f"| STT0000036 池化秩检验 BH 调用率 | {rates(rank, 'spot_rate')} |",
        f"| GSE220442 切片内相邻点中心距（µm） | {min(pitch):.2f}–{max(pitch):.2f} | {nominal.pop():g} |",
        f"| 模拟，无批次关联且每组至少五人的最大修正后零效应拒绝率 | {rates(noise, 'zero_alignment_n_ge_5_max_null_corrected')} |",
        f"| 模拟，论文预设切片的最大修正后零效应拒绝率 | {rates(noise, 'figure_slice_max_null_corrected')} |",
    ]
    doc.write_text(
        before + opening + "\n\n" + "\n".join(rows) + "\n\n" + closing + after, encoding="utf-8"
    )


def main():
    refinement_table()
    for sub in ("figures", "generated"):
        (PAPER / sub).mkdir(parents=True, exist_ok=True)
    public = audits("public")
    demo = audits("synthetic")
    phase = strict_read(RESULTS / "synthetic/phase_diagram.json")
    memory = strict_read(RESULTS / "synthetic/memory.json")
    sota = strict_read(RESULTS / "public/sota.json")["benchmarks"]
    coverage = strict_read(RESULTS / "public/coverage.json")
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "figure.constrained_layout.use": True,
        }
    )
    from build_validation import build_validation

    build_validation(save)
    labels = [r["id"] for r in public]
    x = np.arange(len(public))
    height = max(3.0, len(public) * 0.26 + 0.65)
    fig, ax = plt.subplots(figsize=(6.5, height))
    floors = [r["dimensions"]["power"]["p_floor"] for r in public]
    ax.scatter(floors, x, s=38, color=INK)
    alpha = public[0]["dimensions"]["power"]["alpha"]
    ax.axvline(alpha, color=ACCENT, linestyle=":", label="Nominal threshold")
    ax.set(yticks=x, yticklabels=labels, xscale="log", xlabel="Minimum two-sided p")
    ax.invert_yaxis()
    ax.legend()
    save(fig, "exact-floor")
    fig, ax = plt.subplots(figsize=(6.5, height))
    for i, (key, label, color) in enumerate(
        [("spot_rate", "Pooled observations", ACCENT), ("patient_rate", "Patients", INK)]
    ):
        ax.barh(
            x + (i - 0.5) * 0.3,
            [r["dimensions"]["pseudorep"][key] for r in public],
            0.29,
            color=color,
            label=label,
        )
    ax.set(yticks=x, yticklabels=labels, xlim=(0, 1.08), xlabel="BH call fraction")
    ax.invert_yaxis()
    ax.legend()
    save(fig, "replication")
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    for r in public:
        bands = r["dimensions"]["composition"]["bands"]
        if any(b["share"] is not None for b in bands):
            ax.plot(
                [b["wavelength_um"] for b in bands],
                [b["share"] if b["share"] is not None else np.nan for b in bands],
                marker="o",
                label=r["id"],
            )
    ax.set(xlabel="Wavelength (µm)", ylabel="Explained energy fraction", ylim=(0, 1))
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "composition")
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    for r in public:
        d = r["dimensions"]["nugget"]
        f = [s for s in d["features"] if s["debiased_delta"] is not None and any(s["raw_energy"])]
        if not f:
            continue
        # Divide by each feature's larger group energy to compare across modalities.
        ax.scatter(
            [s["raw_delta"] / max(s["raw_energy"]) for s in f],
            [s["debiased_delta"] / max(s["raw_energy"]) for s in f],
            alpha=0.6,
            s=13,
            label=r["id"],
        )
    ax.axline((0, 0), slope=1, ls=":", color="#777777")
    ax.set(xlabel="Raw contrast / larger group energy", ylabel="Corrected contrast / same scale")
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "nugget")
    sl = phase["figure_slice"]
    rows = [r for r in phase["rows"] if all(r[k] == v for k, v in sl.items())]
    xs = sorted({r["n_per_group"] for r in rows})
    ys = sorted({r["composition_setting"] for r in rows})
    by = {(r["n_per_group"], r["composition_setting"]): r for r in rows}
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.1), sharey=True)
    for ax, target, title in zip(
        axes,
        ("null_corrected", "alternative_corrected"),
        ("Adjusted null rejection", "Adjusted detection"),
    ):
        mesh = ax.contourf(
            xs,
            ys,
            [[by[x, y][target] for x in xs] for y in ys],
            levels=np.linspace(0, 1, 11),
            cmap="viridis",
            vmin=0,
            vmax=1,
        )
        ax.set(xlabel="Patients per group", title=title)
    axes[0].set_ylabel("Composition parameter")
    fig.colorbar(mesh, ax=axes.tolist(), label="Fraction of synthetic cohorts")
    save(fig, "phase-diagram")
    fig, ax = plt.subplots(figsize=(6.5, 3.1))
    names_sorted = sorted({r["method"] for r in sota})
    colors = {name: plt.get_cmap("tab10")(i) for i, name in enumerate(names_sorted)}
    for r in sota:
        if r["status"] == "ok" and r.get("peak_rss_bytes") and r.get("wall_seconds"):
            ax.scatter(
                r["wall_seconds"],
                r["peak_rss_bytes"] / 2**20,
                label=r["method"],
                color=colors[r["method"]],
                s=25,
            )
    handles, names = ax.get_legend_handles_labels()
    legend = dict(zip(names, handles))
    if legend:
        ax.legend(legend.values(), legend.keys(), fontsize=7, ncol=2)
    ax.set(
        xlabel="Elapsed seconds per method run",
        ylabel="Sampled peak RSS (MiB)",
        xscale="log",
        yscale="log",
    )
    save(fig, "benchmark")
    ref = next(r for r in public if r["id"] == "GSE220442")
    macro = {
        "ArchiveCount": coverage["human_records"],
        "MeasuredArchives": len(coverage["audited_accessions"]),
        "MeasuredPlatforms": len({r["platform"] for r in public}),
        "MaxNullSlice": f"{100 * max(r['null_corrected'] for r in rows):.1f}",
        "MaxUnalignedNull": f"{100 * max(r['null_corrected'] for r in phase['rows'] if r['batch_alignment'] == 0 and r['n_per_group'] >= 5):.1f}",
        "PublicRuns": len(public),
        "FixtureCount": len(demo),
        "PhaseCells": phase["scenarios"],
        "SimCohorts": phase["simulated_cohorts"],
        "RepeatCount": phase["replicates_per_state"],
        "MemoryCells": memory["n_observations"],
        "MemoryFeatures": memory["n_features"],
        "MemoryPatients": memory["n_patients"],
        "RSSGiB": f"{memory['peak_rss_bytes'] / 2**30:.3f}",
        "DensityPct": f"{100 * memory['density']:.2f}",
        "ADUnits": ref["n_patients"],
        "ADFeatures": ref["n_features"],
        "ADFloor": f"{ref['dimensions']['power']['p_floor']:.2f}",
        "ADSpot": f"{100 * ref['dimensions']['pseudorep']['spot_rate']:.1f}",
        "ADPatient": f"{100 * ref['dimensions']['pseudorep']['patient_rate']:.1f}",
        "SOTARuns": sum(r["status"] == "ok" for r in sota),
        "SOTAMethods": len({r["method"] for r in sota if r["status"] == "ok"}),
    }
    (PAPER / "generated/metrics.tex").write_text(
        "\n".join(rf"\newcommand{{\{key}}}{{{value}}}" for key, value in macro.items()) + "\n",
        encoding="utf-8",
    )
    table = []
    citations = {
        rid: entry["citation_key"]
        for entry in coverage["by_accession"].values()
        for rid in entry["reports"]
    }
    for r in public:
        d = r["dimensions"]
        table.append(
            f"{r['id']}~\\cite{{{citations[r['id']]}}} & {r['n_patients']} & {r['n_observations']:,} & {r['n_features']} & "
            f"{d['power']['p_floor']:.3g} & {100 * d['pseudorep']['spot_rate']:.1f} & "
            f"{100 * d['pseudorep']['patient_rate']:.1f} " + r"\\"
        )
    (PAPER / "generated/cohorts.tex").write_text(
        "\n".join([*table, r"\bottomrule"]) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
