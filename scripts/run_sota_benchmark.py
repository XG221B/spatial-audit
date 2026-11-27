"""Run external statistical packages on the same registered cohorts and feature panel.

Example: python scripts/run_sota_benchmark.py --manifest cohort.json --out results/public/sota.json
Set R_LIBS to add Bioconductor libraries, or pass --rscript if it is outside PATH.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from dataclasses import replace
from pathlib import Path

# The R model processes are single-threaded. Bound the Python staging libraries
# before NumPy loads its BLAS runtime, respecting an explicit user override.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.io import mmwrite

from spatial_audit.h5_backed import H5Matrix
from spatial_audit.pseudorep import pseudorep_audit
from spatial_audit.schema import (
    Cohort,
    load_manifest,
    sha256,
    strict_read,
    value_blocks,
    write_json,
)

ROOT = Path(__file__).resolve().parent.parent
WORKER = Path(__file__).resolve().with_suffix(".R")
ESTIMANDS = {
    "edgeR": "Patient pseudobulk relative abundance, TMM quasi-likelihood group contrast",
    "DESeq2": "Patient pseudobulk relative abundance, poscounts Wald group contrast",
    "DESpace2": "Differential spatial patterns across annotated domains and patient groups",
    "nnSVG": "Within-section spatial variability of log-normalized expression",
    "SPARK-X": "Within-section spatial variability, mixture of spatial kernels",
    "Cardinal": "Group differences in patient-averaged pixel intensities",
}
SPATIAL = {"nnSVG", "SPARK-X"}
MARGINAL = {"edgeR", "DESeq2", "Cardinal"}
UNASSIGNED = {
    "",
    "na",
    "nan",
    "none",
    "unknown",
    "noise",
    "unassigned",
    "unlabeled",
    "unlabelled",
    "not assigned",
}
NORMALIZATION = {
    "edgeR": "Patient pseudobulk sums, TMM library sizes estimated on the registered panel.",
    "DESeq2": "Patient pseudobulk sums, positive-count size factors estimated on the registered panel.",
    "DESpace2": "Raw counts after removal of unassigned domain observations; DESpace within-domain count model on the registered panel.",
    "SPARK-X": "Raw registered panel counts passed to SPARK-X.",
    "Cardinal": "Archived target-window intensities, patient means with no additional normalization.",
}


def domain_values(section):
    if section.domain is None:
        return np.full(len(section.coords), "", dtype=str)
    labels = pd.Series(section.domain, dtype="string").fillna("").str.strip()
    return labels.mask(labels.str.lower().isin(UNASSIGNED), "").to_numpy(dtype=str)


def depth_source(section):
    if section.library_size is not None:
        return "archived per-observation full library size"
    if isinstance(section.X, H5Matrix):
        return "all columns of the declared H5AD count layer"
    return "registered panel only (panel-closed)"


def normalization(cohort, method):
    if method != "nnSVG":
        return NORMALIZATION[method]
    if cohort.provenance.get("measurement") == "log_normalized":
        return "Deposited log-normalized expression passed directly to nnSVG."
    sources = sorted({depth_source(s) for s in cohort.sections})
    return "log1p(count * 10000 / depth); depth from " + ", ".join(sources) + "."


def depths(section, rows, panel_totals):
    if section.library_size is not None:
        values = np.asarray(section.library_size, dtype=float)
        if values.shape != (len(section.coords),):
            raise ValueError("Library size must have one value per observation")
        result = values[rows]
    elif isinstance(section.X, H5Matrix):
        full = replace(section.X, cols=np.arange(section.X.full_shape[1]))
        chunk = max(1, min(2048, 1_048_576 // full.shape[1]))
        result = np.concatenate(
            [
                np.asarray(full[rows[j : j + chunk], :].sum(axis=1, dtype=float)).ravel()
                for j in range(0, len(rows), chunk)
            ]
        )
    else:
        result = panel_totals
    if not np.isfinite(result).all() or np.any(result < panel_totals):
        raise ValueError("Full library depths must be finite and at least the panel totals")
    return result


def choose_rows(size: int, cap: int | None, index: int, seed: int) -> np.ndarray:
    if cap is None or size <= cap:
        return np.arange(size)
    rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
    return np.sort(rng.choice(size, cap, replace=False))


def export_input(
    cohort: Cohort, directory: Path, cap: int | None, seed: int, use_depth: bool = False
) -> list[int]:
    directory.mkdir(parents=True, exist_ok=True)
    patients = sorted({(s.patient_id, s.group) for s in cohort.sections})
    pd.DataFrame(patients, columns=["patient_id", "group"]).to_csv(
        directory / "patients.csv", index=False
    )
    sizes = []
    for index, section in enumerate(cohort.sections, 1):
        rows = choose_rows(section.X.shape[0], cap, index, seed)
        measurements = sparse.coo_matrix(section.X[rows, :].T)
        mmwrite(directory / f"s{index:05d}.mtx", measurements)
        xy = section.coords[rows, :]
        domain = domain_values(section)[rows]
        panel_totals = np.asarray(measurements.sum(axis=0, dtype=float)).ravel()
        depth = depths(section, rows, panel_totals) if use_depth else panel_totals
        pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "domain": domain, "depth": depth}).to_csv(
            directory / f"s{index:05d}.csv", index=False
        )
        sizes.append(len(rows))
    write_json(
        directory / "input.json",
        {
            "features": cohort.features,
            "n_sections": len(cohort.sections),
            "patients": [s.patient_id for s in cohort.sections],
            "groups": [s.group for s in cohort.sections],
            "log_normalized": cohort.provenance.get("measurement") == "log_normalized",
        },
    )
    return sizes


def inapplicable(cohort: Cohort, method: str, count_scale: bool) -> str | None:
    platform = cohort.platform.upper()
    msi = "MALDI" in platform or "MSI" in platform or "MASS SPECTROM" in platform
    if method == "Cardinal":
        return None if msi else "Cardinal comparator requires an MSI intensity cohort."
    if msi:
        return "RNA comparator targets transcript counts; this cohort contains MSI intensities."
    if not count_scale:
        if method == "nnSVG" and cohort.provenance.get("measurement") == "log_normalized":
            return None
        return (
            "RNA comparator requires an explicitly declared raw-count layer or measurement=counts."
        )
    for section in cohort.sections:
        for arr in value_blocks(section.X):
            if np.any(arr < 0) or np.any(arr != np.floor(arr)):
                return "Count input contains negative or noninteger observations."
    if method == "DESpace2":
        if any(s.domain is None for s in cohort.sections):
            return "Curated domain labels are required for each observation."
        by_group: dict[str, set[str]] = {}
        for section in cohort.sections:
            assigned = set(domain_values(section)) - {""}
            if not assigned:
                return "A section has no assigned domains."
            by_group.setdefault(section.group, set()).update(assigned)
        if len(set.intersection(*by_group.values())) < 2:
            return "The condition-by-domain contrast requires at least two domains shared across groups."
    return None


def sample_rss(process: psutil.Process) -> int:
    try:
        processes = [process, *process.children(recursive=True)]
    except psutil.NoSuchProcess:
        return 0
    rss = 0
    for proc in processes:
        try:
            rss += proc.memory_info().rss
        except (psutil.AccessDenied, psutil.NoSuchProcess):
            continue
    return rss


def terminate_tree(process: psutil.Process) -> None:
    try:
        descendants = process.children(recursive=True)
    except psutil.NoSuchProcess:
        return
    for proc in [*descendants, process]:
        try:
            proc.kill()
        except psutil.NoSuchProcess:
            continue


def invoke(rscript: str, method: str, directory: Path, timeout: float) -> dict:
    output = directory / "worker.json"
    output.unlink(missing_ok=True)
    log = directory / "worker.log"
    worker = directory / "run_sota_benchmark.R"
    shutil.copyfile(WORKER, worker)
    environment = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
    )
    if os.name == "nt":
        for key in ("LC_ALL", "LC_CTYPE"):
            if environment.get(key, "").lower() in {"c.utf8", "c.utf-8"}:
                environment[key] = ".UTF-8"
    started = time.perf_counter()
    peak = 0
    timed_out = False
    with log.open("wb") as handle:
        process = subprocess.Popen(
            [rscript, "--vanilla", str(worker), method, str(directory), str(output)],
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=environment,
        )
        monitor = psutil.Process(process.pid)
        try:
            while process.poll() is None:
                peak = max(peak, sample_rss(monitor))
                if time.perf_counter() - started > timeout:
                    timed_out = True
                    break
                time.sleep(0.03)
        finally:
            if process.poll() is None:
                terminate_tree(monitor)
            process.wait()
    elapsed = time.perf_counter() - started
    result: dict
    if timed_out:
        result = {"status": "timeout", "reason": "Configured per-method wall-time limit reached."}
    elif process.returncode or not output.is_file():
        result = {
            "status": "error",
            "reason": log.read_text(encoding="utf-8", errors="replace")[-900:],
        }
    else:
        result = strict_read(output)
    result.update(wall_seconds=elapsed, peak_rss_bytes=peak or None)
    return result


def summarize(row: dict, features: list[str], spot_rate: float) -> dict:
    if row["status"] != "ok":
        return row
    panels = row["panels"]
    for panel in panels:
        if panel["n_family"] != len(features) or len(panel["p"]) != len(features):
            raise ValueError("Worker returned a different hypothesis family")
        if len(panel["q"]) != len(features):
            raise ValueError("Incomplete adjusted p-value vector")
        raw = np.asarray(panel["p"], dtype=float)
        adjusted = np.asarray(panel["q"], dtype=float)
        available = np.isfinite(raw)
        if not np.array_equal(np.isfinite(adjusted), available) or np.any(
            (raw[available] < 0) | (raw[available] > 1)
        ):
            raise ValueError("Invalid raw or adjusted probabilities")
        observed = raw[available]
        order = np.argsort(observed, kind="stable")
        ranked = observed[order] * len(features) / np.arange(1, len(order) + 1)
        expected = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
        if (
            not np.allclose(adjusted[available][order], expected, rtol=1e-10, atol=1e-14)
            or available.sum() != panel["n_tests"]
            or np.count_nonzero(adjusted < 0.05) != panel["n_called"]
        ):
            raise ValueError("Package report does not reproduce its BH family and calls")
    tests = sum(panel["n_tests"] for panel in panels)
    called = sum(panel["n_called"] for panel in panels)
    opportunities = sum(panel["n_family"] for panel in panels)
    rate = called / opportunities if opportunities and tests else None
    row.update(
        n_tests=tests,
        n_family=opportunities,
        n_called=called,
        call_rate=rate,
        tested_fraction=tests / opportunities if opportunities else None,
    )
    if rate is not None and row["method"] in MARGINAL:
        row["spot_minus_method"] = spot_rate - rate
    if not tests:
        row.update(status="no_estimable_features", reason="Every feature test was unestimable.")
    return row


def benchmark(
    cohort: Cohort,
    count_scale: bool,
    rscript: str | None,
    methods: list[str],
    timeout: float,
    cap: int | None,
    seed: int,
) -> tuple[dict, list[dict]]:
    reference = pseudorep_audit(cohort)
    cohort_info = {
        "id": cohort.id,
        "kind": cohort.kind,
        "platform": cohort.platform,
        "source": cohort.source,
        "features": cohort.features,
        "n_patients": len({s.patient_id for s in cohort.sections}),
        "n_sections": len(cohort.sections),
        "n_observations": sum(len(s.coords) for s in cohort.sections),
        "sections": [
            {
                "section": s.section_id,
                "patient": s.patient_id,
                "group": s.group,
                "observations": len(s.coords),
            }
            for s in cohort.sections
        ],
        "spot_rate": reference["spot_rate"],
        "patient_rate": reference["patient_rate"],
        "spot_minus_patient": reference["gap"],
        "patient_test": reference["patient"]["method"],
    }
    rows = []
    with tempfile.TemporaryDirectory(prefix="spatial-audit-r-") as tmp:
        directory = Path(tmp)
        inputs: dict[tuple, tuple[Path, list[int]]] = {}
        for method in methods:
            row: dict = {
                "cohort_id": cohort.id,
                "kind": cohort.kind,
                "method": method,
                "estimand": ESTIMANDS[method],
                "normalization": normalization(cohort, method),
                "unassigned_domain_observations": sum(
                    int((domain_values(s) == "").sum()) for s in cohort.sections
                )
                if method == "DESpace2"
                else None,
                "version": None,
                "reason": None,
                "n_tests": None,
                "n_called": None,
                "call_rate": None,
                "spot_minus_method": None,
                "peak_rss_bytes": None,
                "wall_seconds": None,
                "panels": [],
            }
            reason = inapplicable(cohort, method, count_scale)
            if reason:
                row.update(status="not_applicable", reason=reason)
            elif rscript is None:
                row.update(status="unavailable", reason="Rscript executable unavailable.")
            else:
                sample_cap = cap if method in SPATIAL else None
                key = (sample_cap, method == "nnSVG")
                if key not in inputs:
                    root = directory / f"input-{len(inputs)}"
                    inputs[key] = root, export_input(cohort, root, sample_cap, seed, key[1])
                root, sizes = inputs[key]
                row.update(
                    input_observations=sum(sizes),
                    per_section_input=sizes,
                    spatial_cap=sample_cap,
                    sampling_seed=seed if sample_cap else None,
                )
                row.update(invoke(rscript, method, root, timeout))
            summarize(row, cohort.features, reference["spot_rate"])
            print(f"{cohort.id} {method}: {row['status']} {row.get('reason') or ''}", flush=True)
            rows.append(row)
    return cohort_info, rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--manifest", type=Path, action="append")
    inputs.add_argument(
        "--release-inputs",
        type=Path,
        metavar="DIR",
        help="Read the published manifest names and spatial caps from the suite registry",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--markdown", type=Path, help="Render the same measured result as a cohort-by-method table"
    )
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--methods", choices=list(ESTIMANDS), nargs="+", default=list(ESTIMANDS))
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument(
        "--spatial-cap",
        type=int,
        help="Optional spot limit per section for within-section SVG only",
    )
    parser.add_argument(
        "--cap-for",
        metavar="COHORT=N",
        action="append",
        default=[],
        help="Override the within-section SVG cap for one registered cohort",
    )
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    if args.timeout <= 0 or (args.spatial_cap is not None and args.spatial_cap < 12):
        parser.error("Timeout must be positive. A spatial cap must be at least 12.")
    if len(set(args.methods)) != len(args.methods):
        parser.error("Each method should appear once.")
    paths = args.manifest
    if args.release_inputs:
        registry = strict_read(ROOT / "docs/data/sota_registry.json")["cohorts"]
        paths = [args.release_inputs / row["manifest"] for row in registry]
        for path, row in zip(paths, registry, strict=True):
            if [s["id"] for s in strict_read(path)["cohorts"]] != [row["id"]]:
                parser.error(f"Registered cohort ID differs from {path.name}")
            if row["spatial_cap"] is not None:
                args.cap_for.insert(0, f"{row['id']}={row['spatial_cap']}")
    jobs = [s for path in paths for s in strict_read(path)["cohorts"]]
    kinds = {s["kind"] for s in jobs}
    ids = {s["id"] for s in jobs}
    if len(kinds) != 1 or len(ids) != len(jobs):
        parser.error("Each output requires one provenance kind and distinct cohort IDs.")
    caps = dict.fromkeys(ids, args.spatial_cap)
    for item in args.cap_for:
        name, separator, number = item.partition("=")
        if not separator or name not in ids or not number.isdigit() or int(number) < 12:
            parser.error("--cap-for requires a registered cohort ID and an integer >= 12")
        caps[name] = int(number)
    kind = kinds.pop()
    if ({"public", "synthetic"} - {kind}).intersection(args.out.parts):
        parser.error("Output path conflicts with the dataset provenance.")
    rscript = shutil.which(args.rscript)
    summary: dict = {
        "kind": kind,
        "schema_version": 1,
        "alpha": 0.05,
        "cohorts": [],
        "benchmarks": [],
        "environment": {
            "os": platform.system(),
            "os_release": platform.release(),
            "architecture": platform.machine(),
            "python": platform.python_version(),
        },
        "multiple_testing": "BH within each prespecified feature panel; untested members remain in family size.",
        "rates": "Fraction called out of the registered family. SVG denominators count feature-section pairs.",
        "timing": "Per-method R process launch through exit, including package loading and input read; "
        "Python staging and reference audit excluded. Peak process-tree summed RSS sampled every 30 ms.",
        "gap_definition": "Spot Wilcoxon call rate minus package call rate for marginal contrasts. "
        "Count methods use patient pseudobulk library normalization on the registered panel.",
    }

    def loaded_jobs():
        for path in paths:
            for spec, c in zip(strict_read(path)["cohorts"], load_manifest(path), strict=True):
                yield path, spec, c

    for path, spec, cohort in loaded_jobs():
        count_scale = spec.get("layer") == "counts" or spec.get("measurement") == "counts"
        info, rows = benchmark(
            cohort, count_scale, rscript, args.methods, args.timeout, caps[cohort.id], args.seed
        )
        info["inputs"] = [
            {"name": Path(f["path"]).name, "sha256": sha256(path.parent / f["path"])}
            for f in spec["files"]
        ]
        info.update({key: spec.get(key) for key in ("measurement", "layer", "selection")})
        info["depth_column"] = spec.get("library_size_column")
        info["depth_inputs"] = [
            {"name": Path(f["library_size_path"]).name, "sha256": f["library_size_sha256"]}
            for f in spec["files"]
            if f.get("library_size_path")
        ]
        summary["cohorts"].append(info)
        summary["benchmarks"].extend(rows)
    write_json(args.out, summary)
    if args.markdown:
        render_table(summary, args.markdown)
    if any(row["status"] in {"error", "timeout"} for row in summary["benchmarks"]):
        raise SystemExit(1)


def render_table(summary: dict, path: Path) -> None:
    """Present recorded calls, family sizes and estimability without comparing distinct nulls."""
    methods = list(ESTIMANDS)
    groups = Counter(row["status"] for row in summary["benchmarks"])
    lookup = {(r["cohort_id"], r["method"]): r for r in summary["benchmarks"]}
    body = [
        "# Cross-platform package results",
        "",
        f"The benchmark covers {len(summary['cohorts'])} {summary['kind']} cohort contrasts. Method outcomes include "
        + ", ".join(
            f"{value} {status.replace('_', ' ')}" for status, value in sorted(groups.items())
        )
        + ".",
        "",
        "Cells show BH calls / registered hypotheses (estimable tests). The SVG columns count feature–section pairs, with separate correction in each section. Other columns count registered cohort features.",
        "",
        "| Contrast | Donors | " + " | ".join(methods) + " |",
        "|---|---:|" + "---:|" * len(methods),
    ]
    for c in summary["cohorts"]:
        cells = []
        for method in methods:
            row = lookup.get((c["id"], method))
            if row is None:
                cells.append("not requested")
            elif row["status"] == "ok":
                cells.append(f"{row['n_called']} / {row['n_family']} ({row['n_tests']})")
            else:
                cells.append(row["status"].replace("_", " "))
        body.append(f"| {c['id']} | {c['n_patients']} | " + " | ".join(cells) + " |")
    body.extend(
        [
            "",
            "Package versions, raw p-values, adjusted p-values, reasons for inapplicability, normalization, timings, peak memory and sampling depths are recorded in the companion JSON.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(body), encoding="utf-8")


if __name__ == "__main__":
    main()
