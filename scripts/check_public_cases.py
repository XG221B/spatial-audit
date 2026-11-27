"""Compare feature-level calls and TIC definitions on the released public panels."""

import argparse
from pathlib import Path

import numpy as np

from spatial_audit.schema import load_manifest, sha256, strict_read, write_json
from spatial_audit.spectral_audit import spectral_audit

PUBLIC = Path(__file__).resolve().parents[1] / "results/public"


def positive_check():
    path = PUBLIC / "GSE250346.json"
    report = strict_read(path)
    pseudo = report["dimensions"]["pseudorep"]
    sota = strict_read(PUBLIC / "sota.json")
    cohort = next(r for r in sota["cohorts"] if r["id"] == report["id"])
    features = pseudo["features"]
    if features != cohort["features"] or report["n_patients"] != cohort["n_patients"]:
        raise ValueError("R contrast and audit must match")
    q = {"patient": pseudo["patient"]["q"]}
    for r in sota["benchmarks"]:
        if r["cohort_id"] == report["id"] and r["method"] in {"edgeR", "DESeq2"}:
            (panel,) = r["panels"]
            q[r["method"]] = panel["q"]
    called = {
        k: [f for f, v in zip(features, values, strict=True) if v is not None and v < 0.05]
        for k, values in q.items()
    }
    overlap = sorted(set.intersection(*(set(called[k]) for k in ("patient", "edgeR", "DESeq2"))))
    return {
        "id": report["id"],
        "n_patients": report["n_patients"],
        "groups": pseudo["groups"],
        "n_features": len(features),
        "alpha": 0.05,
        "called": called,
        "shared_calls": overlap,
        "source": report["source"],
        "report_sha256": sha256(path),
        "sota_sha256": sha256(PUBLIC / "sota.json"),
        "interpretation": "Descriptive overlap on the same donors. Marginal mean-count and "
        "library-normalized pseudobulk tests have different estimands; no ground truth is available.",
    }


def main(manifest):
    path = PUBLIC / "MTBLS12782.json"
    report = strict_read(path)
    (cohort,) = load_manifest(manifest)
    if (
        cohort.id != report["id"]
        or cohort.provenance["inputs"] != report["provenance"]["inputs"]
        or cohort.features != report["dimensions"]["pseudorep"]["features"]
    ):
        raise ValueError("MSI manifest does not match measured inputs")
    records = []
    for loo in (False, True):
        composition, nugget = spectral_audit(cohort, leave_one_out_tic=loo)
        records.append(
            {
                "leave_one_out_tic": loo,
                "shares": [b["share"] for b in composition["bands"]],
                "wavelengths_um": [b["wavelength_um"] for b in composition["bands"]],
                "flagged": nugget["flagged"],
                "evaluated": nugget["evaluated"],
                "features": nugget["features"],
            }
        )
    if not np.allclose(
        records[0]["shares"], [b["share"] for b in report["dimensions"]["composition"]["bands"]]
    ):
        raise ValueError("Published TIC baseline differs")
    write_json(
        PUBLIC / "case_studies.json",
        {
            "kind": "public",
            "record_type": "case_studies",
            "fibrosis": positive_check(),
            "msi": {
                "id": cohort.id,
                "source": cohort.source,
                "report_sha256": sha256(path),
                "records": records,
            },
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--msi-manifest", type=Path, required=True)
    main(p.parse_args().msi_manifest)
