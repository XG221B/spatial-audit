"""Conservative, evidence-local claim extraction from text or text-bearing PDF."""

import hashlib
import re
from pathlib import Path

UNITS = {
    "patient": "patients",
    "participant": "patients",
    "donor": "patients",
    "individual": "patients",
    "animal": "patients",
    "mouse": "patients",
    "mice": "patients",
    "section": "sections",
    "slide": "slides",
    "sample": "samples",
    "specimen": "samples",
    "spot": "observations",
    "cell": "observations",
    "pixel": "observations",
}
STRATEGIES = {
    "batch adjustment": r"batch[- ]correct(?:ed|ion)|adjust(?:ed|ment) for batch|Harmony|ComBat",
    "composition adjustment": r"adjust(?:ed|ment) for (?:cell[- ]type )?composition|deconvol(?:ution|ved)",
    "patient-level analysis": r"pseudo[- ]?bulk|patient[- ]level|donor[- ]level|mixed[- ]effects?",
    "multiple testing": r"Benjamini[-– ]Hochberg|false discovery rate|\bFDR\b",
    "rank test": r"Wilcoxon|Mann[-– ]Whitney",
    "spatial association": r"Moran.s I|spatial autocorrelation",
}
NUMBER_WORDS = {
    w: i
    for i, w in enumerate(
        [
            "zero",
            "one",
            "two",
            "three",
            "four",
            "five",
            "six",
            "seven",
            "eight",
            "nine",
            "ten",
            "eleven",
            "twelve",
        ]
    )
}
_COUNT = re.compile(
    r"\b(?P<count>\d[\d,]*|"
    + "|".join(NUMBER_WORDS)
    + r")\s+(?:(?:unique|independent|adult|human|healthy|control|eligible|distinct|biological)\s+){0,3}(?P<unit>"
    + "|".join(UNITS)
    + r")s?\b",
    re.IGNORECASE,
)
_REVERSE = re.compile(
    r"\b(?P<unit>" + "|".join(UNITS) + r")s?\s*[(,]\s*n\s*=\s*(?P<count>\d[\d,]*)\b", re.IGNORECASE
)


def compare_methods(path: Path, report):
    path = Path(path)
    if path.suffix.lower() == ".pdf":
        from pypdf import PdfReader

        pages = [page.extract_text() or "" for page in PdfReader(path).pages]
    else:
        pages = [path.read_text(encoding="utf-8-sig")]
    if not any(p.strip() for p in pages):
        raise ValueError("No embedded text. Supply a text transcript or searchable PDF")
    if report.get("schema_version") != 1 or report.get("kind") not in {"public", "synthetic"}:
        raise ValueError("A versioned cohort audit report is required")
    measured = {k: report["n_" + k] for k in ("patients", "sections", "observations")}
    d = report["dimensions"]
    audit_evidence = {
        "composition adjustment": {"sections_with_C": d["composition"].get("covariate_sections")},
        "patient-level analysis": {
            "units": report["n_patients"],
            "test": d["pseudorep"]["patient"]["method"],
        },
        "batch adjustment": {"observed_batch_design": d["design"]["categorical"].get("batch")},
        "rank test": {"spot_level": "Mann-Whitney U"},
        "multiple testing": {"procedure": "Benjamini-Hochberg"},
        "spatial association": {
            "waves_um": [b["wavelength_um"] for b in d["composition"]["bands"]]
        },
    }
    matrix = []
    for page, content in enumerate(pages, 1):
        for sentence in re.split(r"(?<=[.!?。])\s+|\n\s*\n", content):
            words = " ".join(sentence.split())
            for match in sorted(
                [*_COUNT.finditer(words), *_REVERSE.finditer(words)], key=lambda m: m.start()
            ):
                unit = UNITS[match["unit"].lower()]
                count = match["count"].lower().replace(",", "")
                value = NUMBER_WORDS[count] if count in NUMBER_WORDS else int(count)
                observed = measured.get(unit)
                matrix.append(
                    {
                        "field": "unit_count",
                        "unit": unit,
                        "claimed": value,
                        "observed": observed,
                        "comparison": "unit_unmapped"
                        if observed is None
                        else (
                            "numeric_agreement_check_scope" if value == observed else "scope_review"
                        ),
                        "page": page,
                        "quoted": words[max(0, match.start() - 90) : match.end() + 90],
                    }
                )
            for strategy, pattern in STRATEGIES.items():
                m = re.search(pattern, words, re.IGNORECASE)
                if m:
                    matrix.append(
                        {
                            "field": "methods_mention",
                            "strategy": strategy,
                            "page": page,
                            "quoted": words[max(0, m.start() - 80) : m.end() + 80],
                            "comparison": "manual_scope_review",
                            "observed": audit_evidence[strategy],
                        }
                    )
    return {
        "schema_version": 1,
        "kind": report["kind"],
        "record_type": "claim_check",
        "cohort": report["id"],
        "methods_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_filename": path.name,
        "extraction": "English explicit numerical units and methods mentions; confirm negation, group, cohort and modality manually.",
        "matrix": matrix,
        "measured": measured,
        "audit_covariate_sections": report["dimensions"]["composition"].get("covariate_sections"),
        "design_flags": report["dimensions"]["design"]["reasons"],
    }
