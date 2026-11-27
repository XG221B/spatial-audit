"""Release checks. A required check that cannot run makes the gate fail."""

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from spatial_audit.schema import strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {
    ".git",
    ".venv",
    ".ruff_cache",
    ".mypy_cache",
    ".pytest_cache",
    "__pycache__",
    "dist",
    "build",
    ".tox",
}


def files():
    return [
        p
        for p in ROOT.rglob("*")
        if p.is_file() and not EXCLUDED.intersection(p.relative_to(ROOT).parts)
    ]


def command(args, cwd=ROOT, env=None):
    completed = subprocess.run(
        args, check=False, cwd=cwd, env=env, capture_output=True, encoding="utf-8", errors="replace"
    )
    if completed.returncode:
        raise RuntimeError(f"{args[0]} failed\n" + (completed.stdout + completed.stderr)[-4500:])
    return completed.stdout


def verify_native_fit(row, family_size, alpha):
    """Independently reconstruct the BH families and reported package call count."""
    import numpy as np
    from scipy.stats import false_discovery_control

    totals = np.zeros(3, dtype=int)
    for panel in row["panels"]:
        p, q = (np.array(panel[k], dtype=float) for k in ("p", "q"))
        if p.shape != (family_size,) or q.shape != p.shape:
            raise ValueError("Native panel has the wrong feature axis")
        finite = np.isfinite(p)
        if not np.array_equal(finite, np.isfinite(q)) or np.any((p[finite] < 0) | (p[finite] > 1)):
            raise ValueError("Inconsistent native feature probabilities")
        # Missing tests still occupy places in the prespecified family.
        padded = np.where(finite, p, 1.0)
        recomputed = false_discovery_control(padded, method="bh")[finite]
        evaluated = int(finite.sum())
        called = int(np.count_nonzero(recomputed < alpha))
        if not np.allclose(q[finite], recomputed, atol=2e-12, rtol=2e-11):
            raise ValueError("Native adjusted probabilities fail independent BH reconstruction")
        counts = [family_size, evaluated, called]
        if counts != [panel[k] for k in ("n_family", "n_tests", "n_called")]:
            raise ValueError("Native family or call counts disagree with raw probabilities")
        totals += counts
    if (
        not np.array_equal(totals, [row[k] for k in ("n_family", "n_tests", "n_called")])
        or not np.isclose(row["call_rate"], totals[2] / totals[0])
        or not all(
            row[k] and np.isfinite(row[k]) and row[k] > 0
            for k in ("wall_seconds", "peak_rss_bytes")
        )
    ):
        raise ValueError("Native call fraction or resource measurement is inconsistent")


def validate_artifacts():
    import json

    from spatial_audit.report import DIMENSIONS
    from spatial_audit.schema import clean

    paths = files()
    json_count = 0
    for path in paths:
        if path.suffix == ".json":
            data = strict_read(path)
            if clean(data) != data:
                raise ValueError(f"Nonfinite value in {path.name}")
            json_count += 1
            parts = path.relative_to(ROOT).parts
            if parts[0] == "results" and (
                parts[1] not in {"synthetic", "public"} or data.get("kind") != parts[1]
            ):
                raise ValueError(f"Provenance mismatch {path.name}")
        if (
            path.suffix
            in {
                ".py",
                ".R",
                ".md",
                ".tex",
                ".toml",
                ".html",
                ".js",
                ".css",
                ".json",
                ".yml",
                ".txt",
            }
            or path.name == "Dockerfile"
        ) and re.search(r"(?<![a-zA-Z0-9])[A-Za-z]:[\\/]", path.read_text(encoding="utf-8-sig")):
            raise ValueError(f"Embedded drive path: {path.relative_to(ROOT)}")
    catalog = strict_read(ROOT / "docs/data/catalog.json")["cohorts"]
    registered = {r["id"] for r in catalog}
    primary = {r["id"] for r in catalog if r["scope"] == "human"}
    if len(primary) < 20 or len(registered) != len(catalog):
        raise ValueError("Insufficient human archive records")
    pages = ROOT / "docs/data"
    idx = strict_read(pages / "audit_summary.json")["reports"]
    keys = set()
    for entry in idx:
        key = (entry["kind"], entry["id"])
        if key in keys or entry["path"] != f"{key[0]}/{key[1]}.json":
            raise ValueError("Ambiguous public/fixture dashboard reference")
        keys.add(key)
        data = strict_read(pages / entry["path"])
        original = strict_read(ROOT / "results" / entry["path"])
        if (
            data != original
            or (data["kind"], data["id"]) != key
            or set(data["dimensions"]) != set(DIMENSIONS)
        ):
            raise ValueError(f"Stale report {key}")
        if key[0] == "public" and (
            not data["source"].startswith("https://")
            or key[1].startswith("synthetic-")
            or data.get("accession") not in registered
        ):
            raise ValueError("Missing public archive provenance")
        if key[0] == "synthetic" and "Synthetic" not in data["statement"]:
            raise ValueError("Missing fixture label")
    measured = {
        (p.parent.name, r["id"])
        for k in ("public", "synthetic")
        for p in (ROOT / "results" / k).glob("*.json")
        if isinstance((r := strict_read(p)).get("dimensions"), dict)
    }
    if measured != keys:
        raise ValueError("Dashboard omits or duplicates audit records")
    for kind in ("public", "synthetic"):
        js = (pages / kind / "bundle.js").read_text(encoding="utf-8")
        payload = json.loads(
            js.removeprefix("globalThis.SPATIAL_" + kind.upper() + "=").removesuffix(";\n")
        )
        if {(kind, r["id"]) for r in payload["reports"]} != {k for k in keys if k[0] == kind}:
            raise ValueError("Offline and JSON dashboard inventories differ")
        for r in payload["reports"]:
            if r != strict_read(pages / kind / (r["id"] + ".json")):
                raise ValueError("Offline report differs from measured result")
        if kind == "public":
            coverage = strict_read(ROOT / "results/public/coverage.json")
            records = {row["id"]: row for row in payload["catalog"]}
            if set(records) != {r["id"] for r in catalog}:
                raise ValueError("Offline catalog differs from registered inventory")
            linked = {
                record["id"]: sorted(record["audit_ids"])
                for record in records.values()
                if record["audit_ids"]
            }
            if set(coverage["audited_accessions"]) != set(linked):
                raise ValueError("Coverage and offline public audits disagree")
            if primary - linked.keys():
                raise ValueError(
                    f"Primary accessions missing measured audits: {sorted(primary - linked.keys())}"
                )
            if {rid for cohort in linked.values() for rid in cohort} != {
                r["id"] for r in payload["reports"]
            }:
                raise ValueError("A published public audit lacks an accession")
            if payload["sota"] != strict_read(ROOT / "results/public/sota.json")["benchmarks"]:
                raise ValueError("Stale offline comparator results")
            for report_id, name in payload["evidence"].items():
                if ("public", report_id) not in keys or Path(name).name != name:
                    raise ValueError("Invalid dashboard preprocessing evidence reference")
                if strict_read(pages / kind / name) != strict_read(ROOT / "results" / kind / name):
                    raise ValueError(f"Outdated input evidence for {report_id}")
            if set(payload["evidence"]) != {r["id"] for r in payload["reports"]}:
                raise ValueError("A measured public report lacks linked archive evidence")
        else:
            grid = strict_read(ROOT / "results/synthetic/phase_diagram.json")
            if payload["phase"] != {k: v for k, v in grid.items() if k != "rows"}:
                raise ValueError("Stale offline synthetic phase diagram metadata")
    memory = strict_read(ROOT / "results/synthetic/memory.json")
    phase = strict_read(ROOT / "results/synthetic/phase_diagram.json")
    if (
        memory["n_observations"] < 10**6
        or memory["peak_rss_bytes"] >= 1.5 * 2**30
        or phase["scenarios"] < 10000
    ):
        raise ValueError("Measured scale targets unmet")
    sota = strict_read(ROOT / "results/public/sota.json")
    from run_sota_benchmark import ESTIMANDS

    coverage = strict_read(ROOT / "results/public/coverage.json")
    owner = {
        report: accession
        for accession, entry in coverage["by_accession"].items()
        for report in entry["reports"]
    }
    comparator_cohorts = {c["id"]: c for c in sota["cohorts"]}
    registry = strict_read(pages / "sota_registry.json")["cohorts"]
    if (
        len(comparator_cohorts) != len(sota["cohorts"])
        or {r["accession"] for r in registry} != primary
        or len(registry) != len(comparator_cohorts)
        or {r["id"]: r["accession"] for r in registry}
        != {k: owner.get(k) for k in comparator_cohorts}
    ):
        raise ValueError("R benchmark inventory differs from the registered primary contrasts")
    for contrast, model in comparator_cohorts.items():
        source = strict_read(ROOT / "results/public" / (contrast + ".json"))
        if (
            model["features"] != source["dimensions"]["pseudorep"]["features"]
            or any(
                model[k] != source[k]
                for k in ("n_patients", "n_sections", "n_observations", "source")
            )
            or {r["patient"] for r in model["sections"]}
            != {s["patient_id"] for s in source["dimensions"]["composition"]["sections"]}
        ):
            raise ValueError(
                f"Registered feature panel or cohort membership differs between Python and R: {contrast}"
            )
    by_contrast: dict[str, set[str]] = {}
    for row in sota["benchmarks"]:
        contrast, method = row["cohort_id"], row["method"]
        group = by_contrast.setdefault(contrast, set())
        if method in group or contrast not in owner or contrast not in comparator_cohorts:
            raise ValueError(f"Duplicated or unregistered comparator: {contrast}, {method}")
        if row["status"] not in {"ok", "not_applicable", "no_estimable_features"}:
            raise ValueError(f"Unresolved comparator: {contrast}, {method}, {row['status']}")
        if row["status"] == "ok":
            verify_native_fit(row, len(comparator_cohorts[contrast]["features"]), sota["alpha"])
        elif not row.get("reason") or row.get("call_rate") is not None:
            raise ValueError(
                f"Unassessed comparison lacks an explicit reason: {contrast}, {method}"
            )
        group.add(method)
    if set(by_contrast) != set(comparator_cohorts) or any(
        methods != ESTIMANDS.keys() for methods in by_contrast.values()
    ):
        raise ValueError("Incomplete set of comparator methods for a registered contrast")
    if primary - {owner[c] for c in by_contrast}:
        raise ValueError("A primary accession lacks a registered R benchmark contrast")
    successful = {row["method"] for row in sota["benchmarks"] if row["status"] == "ok"}
    if len(successful) != 6 or any(r["status"] in {"timeout", "error"} for r in sota["benchmarks"]):
        raise ValueError(f"Comparator validation incomplete: {sorted(successful)}")
    return {
        "json_files": json_count,
        "portable_paths": True,
        "provenance_matched": len(keys),
        "successful_methods": sorted(successful),
        "primary_accessions_measured": len(primary),
        "comparator_contrasts": len(by_contrast),
        "comparator_rows": len(sota["benchmarks"]),
        "independently_verified_r_fits": sum(r["status"] == "ok" for r in sota["benchmarks"]),
    }


def latex():
    command([sys.executable, str(ROOT / "scripts/build_paper.py")])
    outputs = []
    with tempfile.TemporaryDirectory(prefix="spatial-paper-") as directory:
        work = Path(directory) / "paper"
        shutil.copytree(
            ROOT / "paper", work, ignore=shutil.ignore_patterns("main.pdf", "supplement.pdf")
        )
        for name in ("main", "supplement"):
            args = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"]
            command(args, work)
            if name == "main":
                command(["bibtex", name], work)
            command(args, work)
            command(args, work)
            log = (work / (name + ".log")).read_text(encoding="utf-8", errors="replace")
            bad = re.findall(
                r"^.*(?:Overfull \\[hv]box|undefined references|Citation .+ undefined|LaTeX Error).*$",
                log,
                re.MULTILINE,
            )
            if bad:
                raise ValueError(name + " typesetting problems\n" + "\n".join(bad))
            outputs.append(name + ".pdf")
        for document in outputs:
            shutil.copyfile(work / document, ROOT / "paper" / document)
    return outputs


def installed_distribution():
    from tarfile import open as open_tar

    with tempfile.TemporaryDirectory(prefix="spatial-wheel-") as tmp:
        work = Path(tmp)
        command([sys.executable, "-m", "build", "--outdir", str(work), str(ROOT)], work)
        (wheel,) = work.glob("*.whl")
        (sdist,) = work.glob("*.tar.gz")
        with open_tar(sdist) as archive:
            names = {p.name.split("/", 1)[-1] for p in archive.getmembers()}
            for path in (
                "docs/index.html",
                "scripts/run_sota_benchmark.R",
                "paper/main.tex",
                "src/spatial_audit/mcp_server.py",
                "LICENSE",
            ):
                if path not in names:
                    raise ValueError(f"Missing source distribution member: {path}")
        target = work / "installed"
        command(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-deps",
                "--target",
                str(target),
                str(wheel),
            ],
            work,
        )
        env = dict(os.environ, PYTHONPATH=str(target), PYTHONUTF8="1")
        smoke = (
            "from importlib.metadata import distribution; from pathlib import Path; "
            "import spatial_audit, sys; "
            "assert Path(spatial_audit.__file__).resolve().is_relative_to(Path('installed').resolve()); "
            "ep = {e.name: e for e in distribution('spatial-audit').entry_points}; "
            "assert callable(ep['spatial-audit-mcp'].load()); "
            "sys.argv = ['spatial-audit','demo','--out','wheel-results']; ep['spatial-audit'].load()()"
        )
        command([sys.executable, "-c", smoke], work, env)
        if len(list((work / "wheel-results/synthetic").glob("*.json"))) != 4:
            raise ValueError("Installed wheel did not finish demo")
        return {
            "wheel": wheel.name,
            "sdist": sdist.name,
            "installed_demo": "PASS",
            "installed_mcp_entrypoint": "PASS",
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-only", action="store_true", help="CI core checks without LaTeX")
    parser.add_argument(
        "--evidence",
        type=Path,
        help="Save executed check results outside the scientific result folders",
    )
    args = parser.parse_args()
    checks: dict = {}
    with tempfile.TemporaryDirectory(prefix="spatial-tests-") as d:
        junit = Path(d) / "report.xml"
        command([sys.executable, "-m", "pytest", "-q", f"--junitxml={junit}"])
        tests = ET.parse(junit).findall(".//testcase")
        skipped = [t.attrib["name"] for t in tests if t.findall("skipped")]
        if skipped and not args.core_only:
            raise ValueError(f"Full release validation requires executed tests: {skipped}")
        passed = sum(
            len(t.findall("skipped")) == 0
            and len(t.findall("error")) == 0
            and len(t.findall("failure")) == 0
            for t in tests
        )
        if passed < 100:
            raise ValueError("Release needs at least 100 passing tests")
        checks["pytest_passed"] = passed
        checks["pytest_skipped"] = skipped
    print(f"PASS pytest ({passed})", flush=True)
    for name, cmd in [
        ("ruff", ["ruff", "check", "."]),
        ("mypy", ["mypy", "src", "scripts"]),
        ("references", [str(ROOT / "scripts/verify_dois.py"), "--offline"]),
    ]:
        command([sys.executable] + (["-m"] if name != "references" else []) + cmd)
        checks[name] = "PASS"
        print("PASS", name, flush=True)
    checks.update(validate_artifacts())
    print("PASS JSON, paths, provenance and benchmark evidence", flush=True)
    if not args.core_only:
        checks["compiled"] = latex()
        print("PASS LaTeX and layout", flush=True)
        checks["distributions"] = installed_distribution()
        print("PASS sdist and installed wheel", flush=True)
    if args.evidence:
        versions: dict[str, str | None] = {}
        for package in (
            "numpy",
            "scipy",
            "pandas",
            "h5py",
            "anndata",
            "statsmodels",
            "pyarrow",
            "mcp",
            "playwright",
            "ruff",
            "mypy",
            "pytest",
            "matplotlib",
        ):
            try:
                versions[package] = version(package)
            except PackageNotFoundError:
                versions[package] = None
        write_json(
            args.evidence,
            {
                "status": "PASS",
                "scope": "core" if args.core_only else "full",
                "checks": checks,
                "runtime": {
                    "python": platform.python_version(),
                    "os": platform.platform(),
                    "packages": versions,
                },
            },
        )
    print("PASS")


if __name__ == "__main__":
    main()
