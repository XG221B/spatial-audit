"""Verify an archive byte stream against the locally pinned archival fingerprint."""

from __future__ import annotations

import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen

from .schema import sha256


def verified_file(path: Path, url: str, pin: dict, download: bool) -> Path:
    def check(candidate):
        if candidate.stat().st_size != pin["bytes"] or sha256(candidate) != pin["sha256"]:
            raise ValueError(f"Archive fingerprint mismatch for {path.name}")

    if download and not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=path.parent) as temp:
            part = Path(temp) / path.name
            with (
                urlopen(
                    Request(url, headers={"User-Agent": "spatial-audit/0.1"}), timeout=120
                ) as response,
                part.open("wb") as output,
            ):
                shutil.copyfileobj(response, output, 1 << 20)
            check(part)
            part.replace(path)
    else:
        check(path)
    return path
