"""Validate every article in the repository bibliography against PubMed or Crossref.

The checked-in bibliography uses braced fields, one field per line. Unknown entry
syntax is rejected to avoid silently omitting references. No network calls with
private citations or unpublished material are required.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import tempfile
import time
import unicodedata
from datetime import UTC, datetime
from hashlib import sha256
from html import unescape
from pathlib import Path
from typing import Any
from urllib import error, parse, request
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIB = ROOT / "paper" / "references.bib"
DEFAULT_REPORT = ROOT / "paper" / "generated" / "doi_verification.json"
PUBMED = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BIB_ENTRY = re.compile(r"@article\{([A-Za-z][A-Za-z0-9_-]*),\n((?:  \w+ = \{[^\n]*\},?\n)+)\}")
FIELD = re.compile(r"^  (\w+) = \{(.+)\},?$")
REQUIRED = {"doi", "author", "title", "journal", "year"}


def bibliography(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    records: dict[str, dict[str, str]] = {}
    dois: set[str] = set()
    cursor = 0
    for match in BIB_ENTRY.finditer(text):
        if text[cursor : match.start()].strip():
            raise ValueError("Unsupported text or malformed bibliography entry")
        key, fields = match.groups()
        parsed: dict[str, str] = {}
        for line in fields.splitlines():
            field = FIELD.fullmatch(line)
            if field is None or field[1] in parsed:
                raise ValueError(f"Malformed or repeated field in {key}")
            parsed[field[1]] = field[2]
        if REQUIRED - parsed.keys() or key in records:
            raise ValueError(f"Missing fields or repeated citation key: {key}")
        doi = parsed["doi"].casefold()
        if doi in dois or not re.fullmatch(r"10\.\d{4,9}/[^\s{}]+", doi):
            raise ValueError(f"Invalid or repeated DOI in {key}")
        if "pmid" in parsed and not parsed["pmid"].isdecimal():
            raise ValueError(f"Invalid PMID in {key}")
        records[key] = parsed
        dois.add(doi)
        cursor = match.end()
    if not records or text[cursor:].strip():
        raise ValueError("Empty bibliography or unparsed trailing content")
    return records


def fetch(url: str) -> bytes:
    last_error = ""
    for attempt in range(2):
        try:
            req = request.Request(
                url, headers={"User-Agent": "spatial-audit metadata verifier/0.1"}
            )
            with request.urlopen(req, timeout=25) as response:
                return response.read()
        except (OSError, error.URLError) as exc:
            last_error = str(exc)
            if isinstance(exc, error.HTTPError) and exc.code in {429, 503}:
                wait = exc.headers.get("Retry-After", "2")
                if not wait.isdecimal() or int(wait) > 60:
                    raise RuntimeError(f"Retry {url} after {wait}") from exc
                time.sleep(max(int(wait), 1))
            elif attempt == 0:
                time.sleep(1)
    # Windows libcurl can use the OS certificate store when Python uses
    # a different certificate bundle. Both clients verify the remote server.
    curl = shutil.which("curl.exe") or shutil.which("curl")
    if curl:
        return subprocess.check_output([curl, "-fsSL", "--max-time", "30", url])
    raise RuntimeError(f"Could not read {url}: {last_error}")


def clean(text: str) -> str:
    return " ".join(unescape(re.sub(r"<[^>]+>", "", text)).split()).strip().rstrip(".")


def norm(text: str) -> str:
    # Ignore typography and TeX protection of acronyms and special characters.
    folded = unicodedata.normalize("NFKD", clean(text)).casefold()
    return "".join(char for char in folded if char.isalnum())


def metadata_pubmed(payload: bytes) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for record in ET.fromstring(payload).findall("PubmedArticle"):
        pmid = record.findtext("MedlineCitation/PMID", "")
        article = record.find("MedlineCitation/Article")
        if article is None:
            continue
        title = article.find("ArticleTitle")
        doi = record.findtext("PubmedData/ArticleIdList/ArticleId[@IdType='doi']")
        if title is None or not doi:
            continue
        authors = []
        for author in article.findall("AuthorList/Author"):
            if author.findtext("CollectiveName"):
                authors.append("{" + str(author.findtext("CollectiveName")) + "}")
            else:
                last = author.findtext("LastName", "")
                first = author.findtext("ForeName") or author.findtext("Initials", "")
                authors.append(", ".join(name for name in (last, first) if name))
        year = article.findtext("Journal/JournalIssue/PubDate/Year")
        if not year:
            date = article.findtext("Journal/JournalIssue/PubDate/MedlineDate", "")
            match = re.search(r"\d{4}", date)
            year = match[0] if match else ""
        result[pmid] = {
            "doi": doi,
            "title": clean("".join(title.itertext())),
            "author": " and ".join(authors),
            "journal": article.findtext("Journal/Title", ""),
            "year": year,
            "pmid": pmid,
        }
    return result


def metadata_crossref(payload: bytes) -> dict[str, str]:
    row = json.loads(payload)["message"]
    names = []
    for author in row["author"]:
        if author.get("name"):
            names.append("{" + author["name"] + "}")
        else:
            names.append(", ".join(author[k] for k in ("family", "given") if author.get(k)))
    return {
        "doi": row["DOI"],
        "title": clean(row["title"][0]),
        "author": " and ".join(names),
        "journal": clean(row["container-title"][0]),
        "year": str(row["published"]["date-parts"][0][0]),
    }


def verify(bib: Path) -> dict[str, Any]:
    records = bibliography(bib)
    pmids = [row["pmid"] for row in records.values() if "pmid" in row]
    pub_records: dict[str, dict[str, str]] = {}
    for offset in range(0, len(pmids), 100):
        params = parse.urlencode(
            {"db": "pubmed", "retmode": "xml", "id": ",".join(pmids[offset : offset + 100])}
        )
        pub_records.update(metadata_pubmed(fetch(f"{PUBMED}?{params}")))
    evidence = []
    for key, cited in records.items():
        if pmid := cited.get("pmid"):
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            if pmid not in pub_records:
                raise ValueError(f"PubMed did not supply DOI and metadata for {key}")
            authoritative = pub_records[pmid]
        else:
            url = "https://api.crossref.org/works/" + parse.quote(cited["doi"], safe="/")
            authoritative = metadata_crossref(fetch(url))
            time.sleep(0.4)
        mismatches = [
            k
            for k in ("title", "author", "journal", "year")
            if norm(cited[k]) != norm(authoritative[k])
        ]
        if cited["doi"].casefold() != authoritative["doi"].casefold():
            mismatches.append("doi")
        if mismatches:
            raise ValueError(f"{key} differs from {url} in {', '.join(mismatches)}")
        evidence.append({"key": key, "registry_url": url, **authoritative})
    return {
        "schema_version": "1.0",
        "checked_at_utc": datetime.now(UTC).isoformat(),
        "bibliography_sha256": sha256(bib.read_text(encoding="utf-8").encode()).hexdigest(),
        "reference_count": len(evidence),
        "references": evidence,
    }


def check_snapshot(bib: Path, report: Path) -> int:
    """Check exact coverage and fields against a previous live retrieval."""
    records = bibliography(bib)
    evidence = json.loads(report.read_text(encoding="utf-8"))
    # Universal newlines permit the same source to be checked out with CRLF.
    normalized = bib.read_text(encoding="utf-8").encode()
    if evidence["bibliography_sha256"] != sha256(normalized).hexdigest():
        raise ValueError("Bibliography changed since the live verification")
    rows = evidence["references"]
    if len(rows) != len(records) or {r["key"] for r in rows} != set(records):
        raise ValueError("Snapshot does not cover every citation exactly once")
    for row in rows:
        reference = records[row["key"]]
        if any(norm(reference[f]) != norm(row[f]) for f in REQUIRED):
            raise ValueError(f"Metadata mismatch in cached evidence for {row['key']}")
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bib", type=Path, default=DEFAULT_BIB)
    parser.add_argument("--out", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--offline", action="store_true", help="Check against the last registry snapshot"
    )
    args = parser.parse_args()
    if args.offline:
        print(f"Verified snapshot for {check_snapshot(args.bib, args.out)} references")
        return
    result = verify(args.bib)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=args.out.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        temporary.write_text(
            json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.out)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Registry verification passed for {result['reference_count']} references")


if __name__ == "__main__":
    main()
