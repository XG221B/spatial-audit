"""Inspect official GEO records for independent two-group CosMx contrasts."""

import csv
import gzip
import hashlib
import io
import urllib.request
from pathlib import Path

from spatial_audit.schema import write_json

ROOT = Path(__file__).resolve().parents[1]


def download(url):
    with urllib.request.urlopen(url, timeout=90) as response:
        data = response.read()
    return gzip.decompress(data), {
        "url": url,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def main():
    records = []
    reasons = {
        "GSE348495": "GEO titles and the cell table identify tissue specimens. Repeated CCP prefixes and "
        "adjacent tissues require an explicit specimen-to-donor and pairing crosswalk.",
        "GSE312757": "Slides cover tumour and liver regions within colorectal liver metastasis specimens. "
        "The deposited sample characteristics do not define two independent donor groups.",
        "GSE311636": "RNA and protein sample titles identify assays on rectal tissue. The inspected GEO "
        "characteristics do not link these samples to independent donors and response groups.",
    }
    for accession, reason in reasons.items():
        base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}"
        content, receipt = download(f"{base}/soft/{accession}_family.soft.gz")
        lines = content.decode("utf-8").splitlines()

        def fields(key, lines=lines):
            prefix = f"!{key} = "
            return [line.removeprefix(prefix) for line in lines if line.startswith(prefix)]

        row = {
            "accession": accession,
            "assets": [receipt],
            "titles": fields("Sample_title"),
            "characteristics": sorted(set(fields("Sample_characteristics_ch1"))),
            "status": "unregistered",
            "reason": reason,
        }
        if accession == "GSE348495":
            (url,) = (
                s for s in fields("Series_supplementary_file") if s.endswith("_metadata.csv.gz")
            )
            content, meta = download(url.replace("ftp://", "https://", 1))
            with io.TextIOWrapper(io.BytesIO(content), encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                row["cell_metadata_columns"] = reader.fieldnames
                row["tissue_labels"] = sorted({(r["group"], r["sample_name"]) for r in reader})
            row["assets"].append(meta)
        records.append(row)
    write_json(
        ROOT / "results/public/cohort_screening.json",
        {
            "kind": "public",
            "record_type": "candidate_screening",
            "scope": "Targeted review of three human CosMx records",
            "candidates": records,
        },
    )


if __name__ == "__main__":
    main()
