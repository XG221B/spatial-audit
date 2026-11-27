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
    "batch adjustment": r"batch[- ](?:correct(?:ed|ion)|adjustment)|adjust(?:ed|ment) for batch|Harmony|ComBat",
    "composition adjustment": r"adjust(?:ed|ment) for (?:cell[- ]type )?composition|deconvol(?:ution|ved)",
    "patient-level analysis": r"pseudo[- ]?bulk|patient[- ]level|donor[- ]level|mixed[- ]effects?",
    "multiple testing": r"Benjamini[-– ]Hochberg|false discovery rate|\bFDR\b",
    "rank test": r"Wilcoxon|Mann[-– ]Whitney",
    "spatial association": r"Moran.s I|spatial autocorrelation",
}
NUMBER_WORDS = dict(
    zip(
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
            "thirteen",
            "fourteen",
            "fifteen",
            "sixteen",
            "seventeen",
            "eighteen",
            "nineteen",
        ],
        range(20),
        strict=True,
    )
)
for tens, word in enumerate(
    ["twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"], 2
):
    NUMBER_WORDS[word] = tens * 10
    for digit, unit in enumerate(
        ["one", "two", "three", "four", "five", "six", "seven", "eight", "nine"], 1
    ):
        NUMBER_WORDS[f"{word}-{unit}"] = tens * 10 + digit
_NUMBER = (
    r"(?:\d(?:,?\d)*|"
    + "|".join(w.replace("-", "[- ]") for w in sorted(NUMBER_WORDS, key=len, reverse=True))
    + ")"
)
_CN_UNITS = {
    "患者": "patients",
    "受试者": "patients",
    "参与者": "patients",
    "供体": "patients",
    "切片": "sections",
    "载玻片": "slides",
    "样本": "samples",
    "细胞": "observations",
    "像素": "observations",
    "测量点": "observations",
}
_CN_METHODS = {
    "batch adjustment": "批次校正|批次效应校正|校正批次效应",
    "composition adjustment": "细胞组成校正|去卷积|反卷积",
    "patient-level analysis": "患者水平|患者层面|混合效应|伪批量",
    "multiple testing": "多重检验|错误发现率|假发现率",
    "rank test": "秩和检验|秩检验",
    "spatial association": "空间自相关|莫兰指数",
}
_COUNT = re.compile(
    r"\b(?P<count>"
    + _NUMBER
    + r")\s+(?:(?:unique|independent|adult|human|healthy|control|eligible|distinct|biological)\s+){0,3}(?P<unit>"
    + "|".join(UNITS)
    + r")s?\b",
    re.IGNORECASE,
)
_REVERSE = re.compile(
    r"\b(?P<unit>" + "|".join(UNITS) + r")s?\s*[(,]\s*n\s*=\s*(?P<count>" + _NUMBER + r")\b",
    re.IGNORECASE,
)
_CN_COUNT = re.compile(
    r"(?P<count>\d+(?:,\d{3})*|[一二两三四五六七八九]?十[一二三四五六七八九]?|[零一二两三四五六七八九])\s*[名例张个份位只]?\s*(?P<unit>"
    + "|".join(_CN_UNITS)
    + ")"
)
_VS = re.compile(
    r"\bn\s*=\s*(?P<a>"
    + _NUMBER
    + r")\s*(?:vs\.?|versus|对比)\s*(?P<b>"
    + _NUMBER
    + r")(?![A-Za-z0-9])",
    re.IGNORECASE,
)


_CN_DIGITS = {c: n for n, c in enumerate("零一二三四五六七八九")}
_CN_DIGITS["两"] = 2


def _number(text):
    if text in _CN_DIGITS:
        return _CN_DIGITS[text]
    if "十" in text:
        tens, units = text.split("十")
        return 10 * _CN_DIGITS.get(tens, 1) + _CN_DIGITS.get(units, 0)
    word = text.lower().replace(",", "").replace(" ", "-")
    return int(word) if word.isdigit() else NUMBER_WORDS[word]


def _negated(text, match):
    before = re.split(
        r"[,;.!?，；。]|\b(?:but|however)\b|但是|但", text[: match.start()], flags=re.IGNORECASE
    )[-1]
    after = text[match.end() :]
    return bool(
        re.search(
            r"\b(?:no|not|never|without|neither)\b(?:[\s-]+\w+){0,4}\s*$|(?:未|没有|并非|不采用|无)[^，；。]{0,8}$",
            before,
            re.IGNORECASE,
        )
        or re.match(
            r"\s+(?:was|were|is|are)\s+(?:not|never)\b|\s*(?:未采用|未使用)", after, re.IGNORECASE
        )
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
    if report.get("schema_version") != 1 or report.get("kind") not in {
        "public",
        "synthetic",
        "local",
    }:
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
        for sentence in re.split(r"(?<![vV][sS]\.)(?<=[.!?])\s+|(?<=[。！？])|\n\s*\n", content):
            words = " ".join(sentence.split())
            pairs = list(_VS.finditer(words))
            for pair in pairs:
                matrix.append(
                    {
                        "field": "group_counts",
                        "claimed": [_number(pair["a"]), _number(pair["b"])],
                        "observed": [d["power"]["n0"], d["power"]["n1"]],
                        "comparison": "manual_scope_review",
                        "page": page,
                        "quoted": words,
                        "negated": _negated(words, pair),
                    }
                )
            for match in sorted(
                [*_COUNT.finditer(words), *_REVERSE.finditer(words), *_CN_COUNT.finditer(words)],
                key=lambda m: m.start(),
            ):
                if any(match.start() < pair.end() and pair.start() < match.end() for pair in pairs):
                    continue
                unit = {**UNITS, **_CN_UNITS}[match["unit"].lower()]
                value = _number(match["count"])
                negated = _negated(words, match)
                observed = measured.get(unit)
                matrix.append(
                    {
                        "field": "unit_count",
                        "unit": unit,
                        "claimed": value,
                        "negated": negated,
                        "observed": observed,
                        "comparison": "unit_unmapped"
                        if observed is None
                        else (
                            "numeric_agreement_check_scope"
                            if value == observed and not negated
                            else "scope_review"
                        ),
                        "page": page,
                        "quoted": words[max(0, match.start() - 90) : match.end() + 90],
                    }
                )
            for strategy, pattern in STRATEGIES.items():
                for m in re.finditer(pattern + "|" + _CN_METHODS[strategy], words, re.IGNORECASE):
                    matrix.append(
                        {
                            "field": "methods_mention",
                            "strategy": strategy,
                            "negated": _negated(words, m),
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
        "extraction": "English and Chinese unit counts and methods mentions. Negation cues and group scope require manual review.",
        "matrix": matrix,
        "measured": measured,
        "audit_covariate_sections": report["dimensions"]["composition"].get("covariate_sections"),
        "design_flags": report["dimensions"]["design"]["reasons"],
    }
