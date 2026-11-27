"""Fetch listed archive files or independently rehash local copies."""

import argparse
import hashlib
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from filelock import FileLock

from spatial_audit.schema import sha256, strict_read, write_json


def fetch(plan, out, local=False):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for row in plan["files"]:
        name = row["name"]
        if Path(name).name != name or "/" in name or "\\" in name or name in (".", ".."):
            raise ValueError("File name must be a basename")
        url = row["url"]
        if urlsplit(url).scheme != "https":
            raise ValueError("An archive HTTPS URL is required")
        dest = out / name
        tmp = None
        try:
            if local:
                path = dest
            else:
                with tempfile.NamedTemporaryFile(dir=out, delete=False) as handle:
                    tmp = Path(handle.name)
                    with urlopen(
                        Request(url, headers={"User-Agent": "spatial-audit/0.1"}), timeout=120
                    ) as response:
                        shutil.copyfileobj(response, handle, 1024 * 1024)
                        expected = response.headers.get("Content-Length")
                    if expected is not None and handle.tell() != int(expected):
                        raise ValueError("Incomplete transfer")
                path = tmp
            size = path.stat().st_size
            digest = sha256(path)
            if not size or size != row.get("bytes", size) or digest != row.get("sha256", digest):
                raise ValueError(f"Checksum or size mismatch for {name}")
            upstream = row.get("archive_checksum")
            if upstream:
                algorithm, expected_hash = upstream.lower().split(":", 1)
                if algorithm not in {"sha256", "md5"}:
                    raise ValueError("Unsupported archive checksum")
                if algorithm == "sha256":
                    observed = digest
                else:
                    with path.open("rb") as source:
                        observed = hashlib.file_digest(source, "md5").hexdigest()
                if observed != expected_hash:
                    raise ValueError(f"Archive {algorithm} mismatch for {name}")
            receipt = {
                "name": name,
                "accession": row["accession"],
                "url": url,
                "bytes": size,
                "sha256": digest,
                "reference_digest_matched": bool(upstream or "sha256" in row),
                "archive_checksum": upstream,
            }
            manifest = out / "dataset_manifest.json"
            with FileLock(str(manifest) + ".lock", timeout=120):
                previous = strict_read(manifest)["files"] if manifest.is_file() else []
                by_name = {r["name"]: r for r in previous}
                if name in by_name and by_name[name]["accession"] != row["accession"]:
                    raise ValueError("Filename collision between accessions")
                if tmp is not None:
                    tmp.replace(dest)
                elif sha256(dest) != digest:
                    raise ValueError("Local file changed during verification")
                by_name[name] = receipt
                write_json(
                    manifest,
                    {
                        "kind": "public",
                        "files": list(by_name.values()),
                        "checked_utc": datetime.now(UTC).isoformat(),
                    },
                )
            results.append(receipt)
        finally:
            if tmp:
                tmp.unlink(missing_ok=True)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--plan", type=Path)
    source.add_argument("--accession", help="Read asset URLs from the bundled evidence catalog")
    p.add_argument("--role", default="metadata", help="Catalog asset role to transfer")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--local", action="store_true")
    args = p.parse_args()
    if args.plan:
        plan = strict_read(args.plan)
    else:
        catalog = strict_read(Path(__file__).resolve().parents[1] / "docs/data/catalog.json")
        plan = {
            "files": [
                {"accession": c["id"], **a}
                for c in catalog["cohorts"]
                if c["id"] == args.accession
                for a in c["assets"]
                if a["role"] == args.role
            ]
        }
        if not plan["files"]:
            p.error("No catalog assets match the accession and role")
    fetch(plan, args.out, args.local)
