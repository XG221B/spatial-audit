"""Input contract for independent two-group spatial cohorts and strict JSON I/O."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

import numpy as np

from .h5_backed import value_blocks


def strict_read(path: Path) -> Any:
    def invalid(value):
        raise ValueError(f"Invalid JSON number {value}")

    def finite(text):
        value = float(text)
        if not math.isfinite(value):
            invalid(text)
        return value

    return json.loads(
        Path(path).read_text(encoding="utf-8-sig"), parse_constant=invalid, parse_float=finite
    )


def clean(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(item) for item in value]
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        name = Path(tmp.name)
    try:
        name.write_text(
            json.dumps(clean(value), ensure_ascii=False, allow_nan=False, indent=2) + "\n",
            encoding="utf-8",
        )
        os.replace(name, path)
    finally:
        name.unlink(missing_ok=True)


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def archive_accession(declared, source):
    """Read a GEO Series ID only from its official accession URL."""
    url = urlsplit(source or "")
    if (
        url.hostname in {"www.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov"}
        and url.path == "/geo/query/acc.cgi"
    ):
        ids = parse_qs(url.query).get("acc", [])
        if len(ids) == 1 and re.fullmatch(r"GSE[0-9]+", ids[0]):
            if declared and declared != ids[0]:
                raise ValueError("Declared accession disagrees with the official GEO record URL")
            return ids[0]
    return declared


@dataclass
class Section:
    patient_id: str
    group: str
    coords: np.ndarray
    X: Any
    C: Any = None
    slide_id: str | None = None
    batch: str | None = None
    stage: Any = None
    section_id: str = ""
    domain: Any = None
    library_size: Any = None


@dataclass
class Cohort:
    id: str
    kind: str
    features: list[str]
    sections: list[Section]
    platform: str = "unknown"
    source: str = ""
    provenance: dict = field(default_factory=dict)

    design: str = "independent"
    stratum: str = "batch"

    def patient_strata(self):
        if self.stratum not in {"batch", "slide_id"}:
            raise ValueError("stratum must be batch or slide_id")
        strata: dict = {}
        for section in self.sections:
            level = getattr(section, self.stratum)
            if level is None or strata.setdefault(section.patient_id, level) != level:
                raise ValueError("Stratification requires one observed stratum per patient")
        return strata

    def validate(self) -> Cohort:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,119}", self.id):
            raise ValueError("Use a portable cohort ID")
        if self.kind not in {"public", "synthetic", "local"} or (
            self.kind == "public" and not self.source.startswith(("https://", "http://"))
        ):
            raise ValueError("Public cohorts require a source URL; declare provenance kind")
        if not self.features or len(set(self.features)) != len(self.features):
            raise ValueError("Feature names must be nonempty and unique")
        if self.design not in {"independent", "paired", "stratified"}:
            raise ValueError("Unknown permutation design")
        assignments: dict[str, set] = {}
        for s in self.sections:
            if (
                not isinstance(s.patient_id, str)
                or not s.patient_id.strip()
                or not isinstance(s.group, str)
                or not s.group.strip()
            ):
                raise ValueError("Missing patient or group")
            membership = assignments.setdefault(s.patient_id, set())
            membership.add(s.group)
            if len(membership) > 1 and self.design != "paired":
                raise ValueError("Within-patient contrasts require a paired design")
            n = s.X.shape[0]
            if not n or s.X.shape != (n, len(self.features)) or np.shape(s.coords) != (n, 2):
                raise ValueError("Observation or feature axes disagree")
            if not np.isfinite(s.coords).all():
                raise ValueError("Nonfinite coordinates or measurements")
            for array in value_blocks(s.X):
                if not np.isfinite(array).all():
                    raise ValueError("Nonfinite coordinates or measurements")
            if s.C is not None and (
                np.ndim(s.C) != 2 or s.C.shape[0] != n or not np.isfinite(s.C).all()
            ):
                raise ValueError("Invalid covariate matrix")
            if s.stage is not None and (
                np.shape(s.stage) != (2,) or not np.isfinite(s.stage).all()
            ):
                raise ValueError("Stage needs two finite coordinates")
            if s.domain is not None and np.shape(s.domain) != (n,):
                raise ValueError("Domain must label every observation")
            if s.library_size is not None:
                depth = np.asarray(s.library_size, dtype=float)
                if depth.shape != (n,) or not np.isfinite(depth).all() or np.any(depth < 0):
                    raise ValueError(
                        "Library size must be finite, nonnegative and aligned to observations"
                    )
        if self.design == "paired" and any(len(g) != 2 for g in assignments.values()):
            raise ValueError("Paired design requires both conditions for every patient")
        if self.design == "stratified":
            self.patient_strata()
        counts = [g for membership in assignments.values() for g in membership]
        if len(set(counts)) != 2 or min(counts.count(g) for g in set(counts)) < 2:
            raise ValueError("Need two independent groups with at least two patients each")
        return self


def load_manifest(path: Path, *, design=None, stratum=None) -> list[Cohort]:
    path = Path(path).resolve()
    manifest = strict_read(path)
    cohorts = []
    for spec in manifest["cohorts"]:
        if spec.get("units") != "um":
            raise ValueError("Manifest must explicitly declare calibrated units 'um'")
        evidence = spec.get("preprocessing_evidence")
        if evidence is not None and (
            not isinstance(evidence, str)
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.json", evidence)
        ):
            raise ValueError("preprocessing_evidence must be a portable JSON filename")
        features = spec.get("features")
        sections: list[Section] = []
        inputs = []
        seen_paths: set[Path] = set()
        for entry in spec["files"]:
            file = (path.parent / entry["path"]).resolve()
            if file in seen_paths:
                raise ValueError(f"Repeated data file within the same cohort: {file.name}")
            seen_paths.add(file)
            if file.suffix.lower() not in {".npz", ".h5ad"}:
                raise ValueError(f"Unsupported input extension {file.suffix}")
            digest = sha256(file)
            if "sha256" in entry and digest != entry["sha256"]:
                raise ValueError(f"Checksum mismatch: {file.name}")
            inputs.append({"file": file.name, "sha256": digest, "bytes": file.stat().st_size})
            if file.suffix.lower() == ".npz":
                with np.load(file, allow_pickle=False) as data:
                    library_size = data.get("library_size")
                    if entry.get("library_size_path"):
                        depth_file = path.parent / entry["library_size_path"]
                        if sha256(depth_file) != entry.get("library_size_sha256"):
                            raise ValueError("External library size checksum mismatch")
                        external = np.load(depth_file, allow_pickle=False)
                        if library_size is not None and not np.array_equal(library_size, external):
                            raise ValueError("Conflicting library sizes")
                        library_size = external
                    names = data["features"].astype(str).tolist()
                    features = features or names
                    index = {f: i for i, f in enumerate(names)}
                    if len(index) != len(names):
                        raise ValueError("Duplicate feature IDs")
                    matrix = data["X"][:, [index[f] for f in features]]
                    sections.append(
                        Section(
                            entry["patient_id"],
                            entry["group"],
                            data["coords"],
                            matrix,
                            data.get("C"),
                            entry.get("slide_id"),
                            entry.get("batch"),
                            entry.get("stage"),
                            entry.get("section_id", file.stem),
                            data.get("domain"),
                            library_size,
                        )
                    )
            elif file.suffix.lower() == ".h5ad":
                from .h5_backed import read_sections

                features, loaded = read_sections(file, entry, spec, features)
                sections.extend(loaded)
            else:
                raise ValueError(f"Unsupported input extension {file.suffix}")
        cohorts.append(
            Cohort(
                spec["id"],
                spec["kind"],
                features or [],
                sections,
                spec.get("platform", "unknown"),
                spec.get("source", ""),
                {
                    "inputs": inputs,
                    **{
                        key: spec[key]
                        for key in ("fixture_seed", "construction")
                        if spec["kind"] == "synthetic" and key in spec
                    },
                    **({"preprocessing_evidence": evidence} if evidence is not None else {}),
                    "accession": archive_accession(spec.get("accession"), spec.get("source")),
                    "measurement": spec.get("measurement"),
                    "batch_label": spec.get("batch_label", "acquisition batch"),
                    "coordinates": "um",
                    "layer": spec.get("layer", "X"),
                    "covariates": spec.get(
                        "covariates", spec.get("covariate_key", "per-file C if available")
                    ),
                    "selection": spec.get(
                        "selection", "manifest-defined feature and section selection"
                    ),
                },
                design or spec.get("design", "independent"),
                stratum or spec.get("stratum", "batch"),
            ).validate()
        )
    if not cohorts or len({c.id for c in cohorts}) != len(cohorts):
        raise ValueError("Manifest requires unique cohorts")
    return cohorts
