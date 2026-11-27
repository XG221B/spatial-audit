"""Run standardized cohorts, generate the Visium baseline manifest, and build offline dashboard data."""

import argparse
import json
import shutil
from pathlib import Path

from spatial_audit.report import audit, save_report
from spatial_audit.schema import load_manifest, sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "data"


EXPERIMENTS = {
    "public": ("empirical_null", "panel_sensitivity", "case_studies", "cohort_screening"),
    "synthetic": ("vst_calibration", "scaling_curve"),
}


def make_visium(folder, path):
    import os

    files = sorted(Path(folder).glob("*.h5ad"))
    if len(files) != 6:
        raise ValueError("Expected the six per-section GSE220442_ready H5ADs")
    spec = {
        "id": "GSE220442",
        "accession": "GSE220442",
        "kind": "public",
        "platform": "Visium",
        "units": "um",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442",
        "layer": "counts",
        "measurement": "counts",
        "coordinate_calibration": "classic_visium_100um",
        "preprocessing_evidence": "GSE220442_reconstruction.json",
        "library_size_column": "total_counts",
        "selection": "First 32 frozen Visium panel genes, all six donors; physical calibration from classic Visium array addresses or verified 100-um hexagonal positions.",
        "feature_limit": 32,
        "domain": "layer_clean",
        "columns": {"group": "condition"},
        "covariates": [f"ct_Layer_{i}" for i in range(2, 7)] + ["ct_White_Matter"],
    }

    def relative(f):
        try:
            return os.path.relpath(f.resolve(), path.parent.resolve())
        except ValueError:
            return str(f.resolve())

    spec["files"] = [{"path": relative(f), "sha256": sha256(f)} for f in files]
    write_json(path, {"cohorts": [spec]})


def make_xenium(folder, path, all_patients=False):
    """Select the first three GSM files per response group in lexical order."""
    import os
    from collections import Counter

    import h5py
    import numpy as np
    from anndata.io import read_elem

    files = []
    counts: Counter = Counter()
    for f in sorted(Path(folder).glob("GSM*.h5ad")):
        with h5py.File(f) as handle:
            labels = set(np.asarray(read_elem(handle["obs/condition"])).tolist())
        if len(labels) != 1:
            raise ValueError("Expected one patient response per file")
        label = str(labels.pop())
        if not all_patients and counts[label] >= 3:
            continue
        try:
            relative = os.path.relpath(f.resolve(), path.parent.resolve())
        except ValueError:
            relative = str(f.resolve())
        files.append({"path": relative, "sha256": sha256(f)})
        counts[label] += 1
    if sorted(counts.values()) != ([19, 25] if all_patients else [3, 3]):
        raise ValueError("Expected three independent patients in each response group")
    write_json(
        path,
        {
            "cohorts": [
                {
                    "id": "GSE348499" if all_patients else "GSE348499-subset",
                    "kind": "public",
                    "platform": "Xenium",
                    "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348499",
                    "units": "um",
                    "columns": {"group": "condition"},
                    "layer": "counts",
                    "feature_limit": 32,
                    "files": files,
                }
            ]
        },
    )


def publish(results, docs):
    index: dict = {"schema_version": 1, "reports": []}
    embedded: dict = {"reports": [], "catalog": [], "sota": [], "phase": None}
    evidence_links: dict[str, str] = {}
    evidence_by_accession = {}
    experiments: dict[str, dict] = {"public": {}, "synthetic": {}}
    for kind in ("public", "synthetic"):
        for f in sorted((results / kind).glob("*.json")):
            data = strict_read(f)
            if data.get("kind") != kind:
                raise ValueError(f"Provenance mismatch in {f.name}")
            if f.stem in EXPERIMENTS[kind]:
                folder = docs / kind
                folder.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(f, folder / f.name)
                experiments[kind][f.stem] = data
            if kind == "public" and f.stem.endswith(("_inputs", "_reconstruction")):
                dest = docs / kind / f.name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(f, dest)
                evidence_by_accession[data["accession"]] = f.name
            if isinstance(data.get("dimensions"), dict) and data.get("schema_version") == 1:
                dest = docs / kind / f.name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(f, dest)
                index["reports"].append(
                    {"id": data["id"], "kind": kind, "path": f"{kind}/{f.name}"}
                )
                embedded["reports"].append(data)
            if "benchmarks" in data:
                embedded["sota"].extend(data["benchmarks"])
            if f.name == "phase_diagram.json":
                embedded["phase"] = {k: v for k, v in data.items() if k != "rows"}
    if (docs / "catalog.json").exists():
        embedded["catalog"] = strict_read(docs / "catalog.json")["cohorts"]
        acquisitions: dict[str, list[str]] = {}
        for r in embedded["reports"]:
            if r["kind"] == "public":
                accession = match_accession(r, embedded["catalog"])
                if accession:
                    acquisitions.setdefault(accession, []).append(r["id"])
                    if accession in evidence_by_accession:
                        evidence_links[r["id"]] = evidence_by_accession[accession]
        readiness = {}
        for f in sorted((results / "public").glob("*_readiness.json")):
            item = strict_read(f)
            readiness[item["accession"]] = item["preparation_status"]
        embedded["catalog"] = [
            dict(
                row,
                audit_ids=acquisitions.get(row["id"], []),
                preparation_status=readiness.get(row["id"], row["unit_note"]),
            )
            for row in embedded["catalog"]
        ]
        write_json(
            results / "public/coverage.json",
            {
                "kind": "public",
                "record_type": "catalog_coverage",
                "human_records": sum(r["scope"] == "human" for r in embedded["catalog"]),
                "records": len(embedded["catalog"]),
                "audited_accessions": sorted(acquisitions),
                "by_accession": {
                    r["id"]: {
                        "reports": r["audit_ids"],
                        "citation_key": r.get("citation_key"),
                        "status": "measured" if r["audit_ids"] else "catalog_only",
                        "preparation_note": r["preparation_status"],
                    }
                    for r in embedded["catalog"]
                },
            },
        )
    write_json(docs / "audit_summary.json", index)
    # Two separate bundles maintain public/fixture provenance, including on file://.
    for kind in ("public", "synthetic"):
        payload = {
            "reports": [r for r in embedded["reports"] if r["kind"] == kind],
            "experiments": experiments[kind],
            "catalog": embedded["catalog"] if kind == "public" else [],
            "sota": embedded["sota"] if kind == "public" else [],
            "evidence": evidence_links if kind == "public" else {},
            "phase": embedded["phase"] if kind == "synthetic" else None,
        }
        folder = docs / kind
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "bundle.js").write_text(
            "globalThis.SPATIAL_"
            + kind.upper()
            + "="
            + json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
            + ";\n",
            encoding="utf-8",
        )


def match_accession(report, catalog):
    for row in catalog:
        if report.get("accession") == row["id"] or report["id"] == row["id"]:
            return row["id"]
    matches = [r["id"] for r in catalog if r["source"] == report.get("source")]
    if len(matches) > 1:
        raise ValueError("Multiple accessions share an archive URL; supply accession explicitly")
    return matches[0] if matches else None


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path)
    p.add_argument("--visium-dir", type=Path)
    p.add_argument("--xenium-dir", type=Path)
    p.add_argument(
        "--all-patients", action="store_true", help="Use the complete 44-patient Xenium cohort"
    )
    p.add_argument("--results", type=Path, default=ROOT / "results")
    p.add_argument("--docs-data", type=Path, default=DOCS)
    args = p.parse_args()
    if args.xenium_dir:
        if not args.manifest or args.visium_dir:
            p.error("--xenium-dir requires its own --manifest output")
        make_xenium(args.xenium_dir, args.manifest, args.all_patients)
    if args.visium_dir:
        if not args.manifest:
            p.error("--visium-dir requires --manifest as the output manifest path")
        make_visium(args.visium_dir, args.manifest)
    if args.manifest:
        for cohort in load_manifest(args.manifest):
            save_report(audit(cohort), args.results)
    publish(args.results, args.docs_data)
