"""Single-H5AD audit using the same backed section reader as manifest cohorts."""

import re
from pathlib import Path
from tempfile import TemporaryDirectory

from .h5_backed import read_sections
from .report import audit
from .schema import Cohort, sha256


def portable_id(stem):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "-", stem).strip("-_.")
    return name[:120] or "local-cohort"


def audit_anndata(
    adata,
    patient_col="patient_id",
    group_col="group",
    *,
    coord_key="spatial",
    section_col="section_id",
    batch_col="batch",
    slide_col="slide_id",
    layer="X",
    covariate_key=None,
    library_size_col=None,
    cohort_id=None,
    **audit_options,
):
    """Audit AnnData or an H5AD path. Coordinates in obsm[coord_key] must be in µm.

    Scientific provenance is labeled local. In-memory objects are serialized in a
    temporary directory so both entry points use identical H5AD decoding rules.
    """

    def run(path):
        columns = {
            "patient_id": patient_col,
            "group": group_col,
            "section_id": section_col,
            "batch": batch_col,
            "slide_id": slide_col,
        }
        spec = {
            "columns": columns,
            "infer_section_ids": True,
            "coord_key": coord_key,
            "layer": layer,
            "covariate_key": covariate_key,
            "library_size_column": library_size_col,
        }
        features, sections = read_sections(path, {}, spec, None)
        cohort = Cohort(
            cohort_id or portable_id(path.stem),
            "local",
            features,
            sections,
            provenance={
                "coordinates": "um",
                "layer": layer,
                "columns": columns,
                "inputs": [
                    {"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
                ],
            },
        )
        return audit(cohort, **audit_options)

    if isinstance(adata, (str, Path)):
        return run(Path(adata))
    with TemporaryDirectory(prefix="spatial-anndata-") as tmp:
        path = Path(tmp) / "input.h5ad"
        adata.write_h5ad(path)
        return run(path)
