"""Local stdio tools for spatial audit. Start with spatial-audit-mcp."""

from pathlib import Path
from typing import Any, Literal

from mcp.server.fastmcp import FastMCP

from .design_audit import design_audit, exact_floor
from .fixtures import fixtures
from .report import audit
from .schema import clean, load_manifest

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


def main():
    server.run(transport="stdio")


if __name__ == "__main__":
    main()
