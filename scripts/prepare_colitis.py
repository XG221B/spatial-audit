"""Prepare active colonic UC and checkpoint-inhibitor monotherapy biopsies from GSE189184."""

import argparse
import io
import re
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from zipfile import ZipFile

from fetch_and_verify import fetch
from prepare_visium import GEO, GEOMETRY, package_section, read_panel

from spatial_audit.schema import sha256, strict_read, write_json

ACCESSION = "GSE189184"
SOURCE = GEO + ACCESSION
PANEL = [
    "ACTB",
    "GAPDH",
    "RPLP0",
    "PPIA",
    "B2M",
    "EEF1A1",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "MUC2",
    "TFF3",
    "COL1A1",
    "DCN",
    "VIM",
    "PECAM1",
    "VWF",
    "PTPRC",
    "CD3D",
    "CD68",
    "LYZ",
    "HLA-DRA",
    "S100A8",
    "S100A9",
]


class RemoteZip(io.RawIOBase):
    """Seekable HTTP archive. Verify every range against the same file version."""

    def __init__(self, url):
        super().__init__()
        self.url = url
        self.pos = 0
        data, lo, hi, size, version = self.span("bytes=-65536")
        self.size, self.version = size, version
        self.cached = (lo, hi, data)

    def span(self, header):
        for attempt in range(5):
            try:
                return self.request_span(header, attempt)
            except (TimeoutError, URLError) as error:
                if isinstance(error, HTTPError) and error.code not in (429, 500, 502, 503, 504):
                    raise
                if attempt == 4:
                    raise
                time.sleep(2**attempt)
        raise RuntimeError("Could not retrieve archive range")

    def request_span(self, header, attempt):
        req = Request(
            self.url + f"?spatial_range={header[6:]}&attempt={attempt}", headers={"Range": header}
        )
        with urlopen(req, timeout=100) as response:
            match = re.fullmatch(
                r"bytes (\d+)-(\d+)/(\d+)", str(response.headers.get("Content-Range"))
            )
            if response.status != 206 or not match:
                raise ValueError("Verified HTTP byte-range support is required")
            lo, hi, size = map(int, match.groups())
            version = response.headers.get("Last-Modified")
            data = response.read(hi - lo + 2)
            if len(data) != hi - lo + 1 or not version:
                raise ValueError("Incomplete archive range or missing version")
        return data, lo, hi, size, version

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = (0, self.pos, self.size)[whence] + offset
        return self.pos

    def seekable(self):
        return True

    def read(self, n=-1):
        n = max(0, self.size - self.pos) if n < 0 else min(n, self.size - self.pos)
        if n <= 0:
            return b""
        lo, hi, data = self.cached
        if lo > self.pos or hi < self.pos + n - 1:
            end = min(self.size, self.pos + max(65536, n)) - 1
            data, lo, hi, size, version = self.span(f"bytes={self.pos}-{end}")
            if (lo, hi, size, version) != (self.pos, end, self.size, self.version):
                raise ValueError("Source archive or requested extent changed")
            self.cached = (lo, hi, data)
        p = self.pos - lo
        self.pos += n
        return data[p : p + n]


def biopsy_design(text):
    rows = []
    for sample in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        fields = dict(re.findall(r"(?m)^!Sample_characteristics_ch1 = ([^:\n]+): (.+)$", sample))
        disease = fields["disease"]
        if fields["section cut"] != "Biopsy" or fields["activity"] != "Inflamed":
            continue
        if (
            disease == "Checkpoint inhibitor induced colitis"
            and fields["immunotherapy"] == "Monotherapy"
        ):
            group = "ICI_monotherapy"
        elif disease == "Ulcerative Colitis" and fields["immunotherapy"] == "N/A":
            group = "UC"
        else:
            continue
        gsm = sample.splitlines()[0].strip()
        title = re.search(r"!Sample_title = (.+)", sample)
        link = re.search(r"!Sample_supplementary_file_1 = (\S+\.zip)", sample)
        if not title or not link or not re.fullmatch(r"GSM\d+", gsm):
            raise ValueError("Missing GSM, biopsy identifier or ZIP")
        patient = re.search(r"(?:GI\s*\d+|PRISE\s*\d+|BB\s*\d+)(?=[ )]|$)", title[1])
        if not patient:
            raise ValueError("Could not resolve archived participant identifier")
        rows.append(
            {
                "section_id": gsm,
                "patient_id": patient[0].replace(" ", ""),
                "group": group,
                "url": link[1].replace("ftp://", "https://"),
            }
        )
    if len({r["patient_id"] for r in rows}) != len(rows):
        raise ValueError("Resolve repeated participant before independent contrast")
    return rows


def extract(row, cache):
    folder = cache / row["section_id"]
    folder.mkdir(exist_ok=True)
    receipt = folder / "dataset_manifest.json"
    if receipt.exists():
        previous = strict_read(receipt)
        if previous["url"] == row["url"] and all(
            sha256(folder / f["name"]) == f["sha256"] for f in previous["members"]
        ):
            return previous
    members = []
    with RemoteZip(row["url"]) as stream, ZipFile(stream) as archive:
        for name in ("filtered_feature_bc_matrix.h5", "tissue_positions_list.csv"):
            matches = [f for f in archive.infolist() if PurePosixPath(f.filename).name == name]
            if len(matches) != 1:
                raise ValueError("Expected one archived " + name)
            f = matches[0]
            with tempfile.NamedTemporaryFile(dir=folder, delete=False) as temp:
                path = Path(temp.name)
            try:
                with archive.open(f) as member, path.open("wb") as out:
                    while block := member.read(2**20):
                        out.write(block)
                if path.stat().st_size != f.file_size:
                    raise ValueError("Incomplete ZIP member")
                digest = sha256(path)
                path.replace(folder / name)
            finally:
                path.unlink(missing_ok=True)
            members.append(
                {
                    "name": name,
                    "member": f.filename,
                    "bytes": f.file_size,
                    "sha256": digest,
                    "crc32": f"{f.CRC:08x}",
                }
            )
        result = {
            "url": row["url"],
            "archive_bytes": stream.size,
            "last_modified": stream.version,
            "members": members,
        }
    write_json(receipt, result)
    return result


def prepare(cache, manifest, provenance):
    source = SOURCE + "&targ=gsm&form=text&view=full"
    meta = fetch(
        {"files": [{"name": "samples.soft", "url": source, "accession": ACCESSION}]},
        cache,
        local=(cache / "samples.soft").exists(),
    )
    design = biopsy_design((cache / "samples.soft").read_text(encoding="utf-8"))
    if sorted([sum(r["group"] == g for r in design) for g in ("ICI_monotherapy", "UC")]) != [4, 5]:
        raise ValueError("Unexpected eligible donor inventory")
    with ThreadPoolExecutor(3) as pool:
        members = list(pool.map(lambda r: extract(r, cache), design))
    sections = []
    out = manifest.parent / "colitis-sections"
    out.mkdir(exist_ok=True, parents=True)
    for row in design:
        folder = cache / row["section_id"]
        x, xy, library_size = read_panel(
            folder / "filtered_feature_bc_matrix.h5",
            (folder / "tissue_positions_list.csv").read_bytes(),
            PANEL,
            with_depth=True,
        )
        section = {k: v for k, v in row.items() if k != "url"}
        sections.append(package_section(out, manifest, section, x, xy, PANEL, library_size))
    selection = "All inflamed colonic biopsies classified as UC (no immunotherapy) or checkpoint-inhibitor colitis on monotherapy. Unique archived GI/PRISE/BB identifiers. Raw counts at in-tissue spots; nominal Visium lattice. Treatment exposure is part of the contrast."
    write_json(
        manifest,
        {
            "cohorts": [
                {
                    "id": ACCESSION,
                    "accession": ACCESSION,
                    "platform": "Visium",
                    "kind": "public",
                    "source": SOURCE,
                    "units": "um",
                    "measurement": "counts",
                    "files": sections,
                    "selection": selection,
                }
            ]
        },
    )
    write_json(
        provenance,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "accession": ACCESSION,
            "source": SOURCE,
            "metadata": meta,
            "archives": members,
            "donors": design,
            "features": PANEL,
            "physical_calibration": GEOMETRY,
            "selection": selection,
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for option in ("cache", "manifest", "provenance"):
        p.add_argument("--" + option, required=True, type=Path)
    args = p.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
