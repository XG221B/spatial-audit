"""Build a reproducible positive-mode MSI panel from the public MTBLS12782 assay.

Download the catalog's metadata and data roles with fetch_and_verify.py, then
pass their common cache directory here. The imzML UUID and embedded IBD MD5
are checked before any spectra are interpreted. Source and section SHA-256
digests are recorded independently.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from pathlib import Path
from uuid import UUID
from xml.etree import ElementTree as ET

import numpy as np

from spatial_audit.schema import sha256, strict_read, write_json

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "MTBLS12782"
SOURCE = "https://www.ebi.ac.uk/metabolights/MTBLS12782"
PPM = 20.0
PANEL_SIZE = 32


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def mass_panel(maf: Path) -> np.ndarray:
    """Use labeled masses only, with deterministic coverage of their sorted range."""
    masses = sorted({float(r["mass_to_charge"]) for r in read_tsv(maf) if r["mass_to_charge"]})
    if len(masses) < PANEL_SIZE or not np.isfinite(masses).all() or masses[0] <= 0:
        raise ValueError("Invalid or insufficient positive m/z targets in the archive MAF")
    indices = np.rint(np.linspace(0, len(masses) - 1, PANEL_SIZE)).astype(int)
    targets = np.array(masses)[indices]
    if np.any(targets[:-1] * (1 + PPM / 1e6) >= targets[1:] * (1 - PPM / 1e6)):
        raise ValueError("Selected m/z windows overlap")
    return targets


def cv(element: ET.Element) -> dict[str, str]:
    return {
        child.attrib["accession"]: child.get("value", "") for child in element.findall("{*}cvParam")
    }


def check_pair(header: ET.Element, binary_path: Path) -> None:
    columns = cv(header)
    uid = UUID(columns["IMS:1000080"])
    checksum = columns["IMS:1000090"].lower()
    with binary_path.open("rb") as stream:
        if stream.read(16) != uid.bytes:
            raise ValueError(f"UUID mismatch for {binary_path.name}")
        stream.seek(0)
        if hashlib.file_digest(stream, "md5").hexdigest() != checksum:
            raise ValueError(f"Embedded MD5 mismatch for {binary_path.name}")


def binary_array(element: ET.Element, binary: np.memmap, types: dict) -> tuple[str, np.ndarray]:
    ref = element.find("{*}referenceableParamGroupRef")
    if ref is None or ref.attrib["ref"] not in types:
        raise ValueError("Unknown imzML external array type")
    kind, dtype = types[ref.attrib["ref"]]
    values = cv(element)
    count = int(values["IMS:1000103"])
    offset = int(values["IMS:1000102"])
    size = int(values["IMS:1000104"])
    if count < 0 or size != dtype.itemsize * count or offset < 16 or offset + size > binary.nbytes:
        raise ValueError("Invalid imzML external array bounds")
    return kind, np.frombuffer(binary, dtype=dtype, count=count, offset=offset)


def pixel(spectrum, binary, types, lower, upper):
    scan = spectrum.find("{*}scanList/{*}scan")
    if scan is None or "MS:1000130" not in cv(spectrum):
        raise ValueError("Expected a positive-mode pixel")
    position = cv(scan)
    x, y = int(position["IMS:1000050"]), int(position["IMS:1000051"])
    if x < 1 or y < 1:
        raise ValueError("Pixel indices must be positive")
    arrays = dict(
        binary_array(child, binary, types)
        for child in spectrum.findall("{*}binaryDataArrayList/{*}binaryDataArray")
    )
    mz, intensity = arrays["mz"], arrays["intensity"]
    if len(mz) != len(intensity) or not np.isfinite(mz).all() or not np.isfinite(intensity).all():
        raise ValueError("Invalid centroid spectrum")
    if np.any(intensity < 0) or np.any(mz[:-1] > mz[1:]):
        raise ValueError("Expected positive-mode sorted centroids and nonnegative intensities")
    # Sum all centroids within each closed ppm window. Prefix differences
    # use float64 so small peaks retain precision beside the higher intensities.
    prefix = np.empty(len(intensity) + 1, dtype=np.float64)
    prefix[0] = 0
    np.cumsum(intensity, dtype=np.float64, out=prefix[1:])
    selected = prefix[mz.searchsorted(upper, side="right")] - prefix[mz.searchsorted(lower)]
    return (x, y), selected, prefix[-1]


def extract(xml: Path, ibd: Path, targets: np.ndarray, spacing: tuple[float, float]):
    binary = np.memmap(ibd, dtype="uint8", mode="r")
    lower, upper = targets * (1 - PPM / 1e6), targets * (1 + PPM / 1e6)
    types: dict = {}
    expected = None
    validated_pair = False
    calibrated = False
    xy, spectra, tic = [], [], []
    for event, node in ET.iterparse(xml, events=("start", "end")):
        tag = node.tag.rsplit("}", 1)[-1]
        if event == "start":
            if tag == "spectrumList":
                expected = int(node.attrib["count"])
            continue
        if tag == "fileContent":
            check_pair(node, ibd)
            validated_pair = True
        elif tag == "referenceableParamGroup":
            codes = cv(node)
            kind = "mz" if "MS:1000514" in codes else "intensity" if "MS:1000515" in codes else None
            if kind is not None:
                dtype = "<f4" if "MS:1000521" in codes else "<f8" if "MS:1000523" in codes else None
                if not dtype or "MS:1000576" not in codes or "IMS:1000101" not in codes:
                    raise ValueError("Unsupported encoding for this MTBLS12782 adapter")
                types[node.attrib["id"]] = (kind, np.dtype(dtype))
        elif tag == "scanSettings":
            codes = cv(node)
            calibrated = np.allclose(
                [float(codes["IMS:1000046"]), float(codes["IMS:1000047"])], spacing
            )
        elif tag == "spectrum":
            if not validated_pair or not calibrated:
                raise ValueError("Pair identity or pixel calibration unavailable")
            coords, selected, total = pixel(node, binary, types, lower, upper)
            xy.append(coords)
            spectra.append(selected)
            tic.append(total)
            node.clear()
    if expected != len(xy) or not expected or len(set(xy)) != len(xy):
        raise ValueError("Spectrum count mismatch or duplicate pixel coordinates")
    um = (np.asarray(xy, dtype=float) - 1) * np.asarray(spacing)
    design = np.column_stack([np.asarray(tic), um - um.mean(axis=0)])
    return um, np.asarray(spectra, dtype=np.float32), design


def read_assay(path: Path) -> list[dict[str, str]]:
    # ISA repeats column names (including 'Unit'). Pair each pixel-size
    # column with its own unit while retaining globally unique assay fields.
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        ix = [header.index(f"Parameter Value[Pixel size {axis}]") for axis in "xy"]
        rows = []
        for line in reader:
            if len(line) != len(header) or any(line[i + 1] != "micrometer" for i in ix):
                raise ValueError("Invalid ISA row or unsupported pixel-size units")
            rows.append(dict(zip(header, line, strict=True)))
    return rows


def prepare(cache: Path, manifest: Path, evidence: Path) -> None:
    cat = strict_read(ROOT / "docs/data/catalog.json")
    cohort = next(c for c in cat["cohorts"] if c["id"] == ACCESSION)
    assets = {a["name"]: a for a in cohort["assets"]}
    inputs = {}
    for name in assets:
        path = cache / name
        inputs[name] = {
            "name": name,
            "url": assets[name]["url"],
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "accession": ACCESSION,
            "role": assets[name]["role"],
        }
    sample = cache / f"s_{ACCESSION}.txt"
    samples = read_tsv(sample)
    groups = {r["Sample Name"]: r["Factor Value[Cohort]"] for r in samples}
    if len(groups) != len(samples) or set(groups.values()) != {"Normal", "diabetes mellitus"}:
        raise ValueError("Missing, repeated, or unexpected sample identifiers/groups")
    assay = read_assay(cache / f"a_{ACCESSION}_MSImaging___metabolite_profiling.txt")
    if len(assay) != len(groups) or {r["Sample Name"] for r in assay} != set(groups):
        raise ValueError("Imaging assay and sample table differ; review the assay-to-donor mapping")
    out = manifest.parent / (manifest.stem + ".sections")
    out.mkdir(parents=True, exist_ok=True)
    sections, measurements = [], []
    targets = mass_panel(cache / f"m_{ACCESSION}_MSImaging___metabolite_profiling_v2_maf.tsv")
    features = [f"mz_{v:.5f}" for v in targets]
    for row in assay:
        sample_id = row["Sample Name"]
        if row["Parameter Value[Scan polarity]"] != "positive":
            raise ValueError("Mixed or unsupported polarity")
        spacing = tuple(float(row[f"Parameter Value[Pixel size {axis}]"]) for axis in "xy")
        if not all(np.isfinite(v) and v > 0 for v in spacing):
            raise ValueError("Invalid physical pixel pitch")
        ibd, xml = (
            cache / Path(row[k]).name
            for k in ("Raw Spectral Data File", "Derived Spectral Data File")
        )
        for path in (xml, ibd):
            if path.name not in inputs:
                raise ValueError(f"Input absent from official asset catalog: {path.name}")
        coords, x, c = extract(xml, ibd, targets, (spacing[0], spacing[1]))
        path = out / (ibd.stem + ".npz")
        with tempfile.NamedTemporaryFile(dir=out, suffix=".npz", delete=False) as handle:
            temp = Path(handle.name)
        try:
            np.savez_compressed(temp, coords=coords, X=x, C=c, features=np.asarray(features))
            os.replace(temp, path)
        finally:
            temp.unlink(missing_ok=True)
        group = {"Normal": "PDAC without diabetes", "diabetes mellitus": "PDAC with diabetes"}[
            groups[sample_id]
        ]
        sections.append(
            {
                "path": path.relative_to(manifest.parent).as_posix(),
                "patient_id": sample_id,
                "section_id": sample_id,
                "group": group,
                "sha256": sha256(path),
            }
        )
        measurements.append(
            {
                "section_id": sample_id,
                "group": group,
                "pixels": len(coords),
                "pixel_pitch_um": list(spacing),
                "source_xml": xml.name,
                "source_ibd": ibd.name,
                "section_bytes": path.stat().st_size,
                "section_sha256": sections[-1]["sha256"],
                "signal_nonzero_pixels": (x > 0).sum(axis=0),
            }
        )
    policy = {
        "interpretation": "Subset of archived positive-mode m/z features. Intensities are as exported in SCiLS imzML.",
        "selection": "Equally spaced rank quantiles of unique positive MAF mass_to_charge values; no labels, intensities or fold changes used in rank selection. The archive MAF itself is study-specific.",
        "target_mz": targets,
        "window_ppm_each_side": PPM,
        "abundance": "Sum of all matching centroid intensities, all acquired pixels, zero for an empty window. Chemical identities are unverified.",
        "covariates": [
            "full positive-spectrum total ion current",
            "within-section centered raster x (um)",
            "within-section centered raster y (um)",
        ],
        "tic_interpretation": "Each extracted window contributes to the all-peak TIC. The explained fraction includes this algebraic dependence and spatial gradients correlated with tissue structure.",
        "spatial_frame": "Zero-based local pixel centers from imzML raster coordinates and ISA pitch. Between-run absolute stage registration and batch/slide assignments are unavailable.",
        "patient_mapping": "One positive-mode imaging sample per study participant in the spatial subset, per the paper's cohort description and the ISA sample-to-assay crosswalk.",
        "noise": "Continuous centroid intensities. Nugget estimates quantify sensitivity under working independent pixel noise; an ion peak does not imply an identified molecule.",
    }
    spec = {
        "id": ACCESSION,
        "accession": ACCESSION,
        "kind": "public",
        "source": SOURCE,
        "platform": "MALDI-MSI",
        "measurement": "intensity",
        "units": "um",
        "preprocessing_evidence": evidence.name,
        "features": features,
        "files": sections,
        "layer": "positive-mode centroid intensities summed in closed target windows",
        "covariates": policy["covariates"],
        "selection": policy["selection"],
    }
    write_json(manifest, {"cohorts": [spec]})
    write_json(
        evidence,
        {
            "kind": "public",
            "accession": ACCESSION,
            "source": SOURCE,
            "paper": cohort["paper"],
            "files": list(inputs.values()),
            "preprocessing": policy,
            "sections": measurements,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache",
        type=Path,
        required=True,
        help="Flat directory of archives downloaded for both catalog roles",
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())


if __name__ == "__main__":
    main()
