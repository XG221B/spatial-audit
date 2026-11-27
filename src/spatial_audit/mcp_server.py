"""Local stdio tools for spatial audit. Start with spatial-audit-mcp."""

from pathlib import Path
from typing import Any, Literal

from mcp.server.fastmcp import FastMCP

from .api import audit_anndata
from .claim_check import compare_methods
from .design_audit import design_audit, exact_floor
from .fixtures import fixtures
from .report import audit
from .schema import clean, load_manifest, strict_read

server = FastMCP("spatial-audit")


@server.tool()
def audit_cohort(manifest: str, cohort_id: str) -> dict[str, Any]:
    """Audit one named cohort from a trusted local input manifest; return five dimensions."""
    for cohort in load_manifest(Path(manifest)):
        if cohort.id == cohort_id:
            return audit(cohort)
    raise ValueError("Cohort ID absent from manifest")


@server.tool()
def check_design_collinearity(manifest: str, cohort_id: str) -> dict[str, Any]:
    """Evaluate acquisition aliasing and physical coordinate overlap for a local cohort."""
    for cohort in load_manifest(Path(manifest)):
        if cohort.id == cohort_id:
            return clean(design_audit(cohort))
    raise ValueError("Cohort ID absent from manifest")


@server.tool()
def compute_exact_p_floor(
    n0: int, n1: int, features: int = 1, alpha: float = 0.05
) -> dict[str, Any]:
    """Resolution and feasible BH rank for the doubled-tail randomization test."""
    return exact_floor(n0, n1, features, alpha)


@server.tool()
def simulate_trap_cohort(
    trap: Literal["null", "pseudorep", "composition", "stage"],
) -> dict[str, Any]:
    """Generate and evaluate a named, deterministic Synthetic Benchmark Fixture in memory."""
    cohort = next(c for c in fixtures() if c.id == "synthetic-" + trap)
    return audit(cohort)


@server.tool()
def quick_audit_h5ad(
    path: str,
    patient_col: str,
    group_col: str,
    alpha: float = 0.05,
    design: Literal["independent", "paired", "stratified"] = "independent",
    stratum: Literal["batch", "slide_id"] = "batch",
) -> dict[str, Any]:
    """Audit a local H5AD with calibrated micrometre coordinates stored in obsm/spatial."""
    return audit_anndata(
        Path(path), patient_col, group_col, alpha=alpha, design=design, stratum=stratum
    )


@server.tool()
def claim_check_methods(methods: str, report: str) -> dict[str, Any]:
    """Compare a local methods text or searchable PDF with a saved audit report."""
    return compare_methods(Path(methods), strict_read(Path(report)))


@server.tool()
def list_manifest_cohorts(manifest: str) -> dict[str, Any]:
    """List validated cohort IDs, study designs and biological-unit counts."""
    return {
        "cohorts": [
            {
                "id": c.id,
                "design": c.design,
                "sections": len(c.sections),
                "patients": len({s.patient_id for s in c.sections}),
            }
            for c in load_manifest(Path(manifest))
        ]
    }


def main():
    server.run(transport="stdio")


if __name__ == "__main__":
    main()
