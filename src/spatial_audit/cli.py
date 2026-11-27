"""spatial-audit command line."""

import argparse
from pathlib import Path

from .report import audit, save_report
from .schema import load_manifest


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spatial-audit")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate-manifest")
    validate.add_argument("--manifest", type=Path, required=True)
    for command in ("demo", "run", "quick-run"):
        sub = commands.add_parser(command)
        if command == "run":
            sub.add_argument("--manifest", type=Path, required=True)
        elif command == "quick-run":
            sub.add_argument("h5ad", type=Path, help="H5AD with calibrated micrometre coordinates")
            sub.add_argument("--patient-col", required=True)
            sub.add_argument("--group-col", required=True)
            sub.add_argument("--section-col", default="section_id")
            sub.add_argument("--batch-col", default="batch")
            sub.add_argument("--slide-col", default="slide_id")
            sub.add_argument("--coord-key", default="spatial", help="obsm coordinate key (µm)")
            sub.add_argument("--layer", default="X")
            sub.add_argument("--covariate-key")
            sub.add_argument("--library-size-col")
            sub.add_argument("--cohort-id")
        else:
            sub.add_argument(
                "--write-fixtures", type=Path, help="Export reusable synthetic NPZs and manifest"
            )
        sub.add_argument("--out", type=Path, default=Path("results"))
        sub.add_argument(
            "--wavelengths",
            nargs="+",
            default=["250", "600", "1400"],
            help="Increasing wavelengths in µm, or auto",
        )
        sub.add_argument("--alpha", type=float, default=0.05)
        sub.add_argument("--design", choices=("independent", "paired", "stratified"))
        sub.add_argument("--stratum", choices=("batch", "slide_id"))
        sub.add_argument("--nugget-mode", choices=("raw", "vst"), default="raw")
        sub.add_argument("--graph-mode", choices=("grid", "voronoi"), default="grid")
        sub.add_argument(
            "--leave-one-out-tic",
            action="store_true",
            help="Covariate column zero must be full TIC",
        )
        sub.add_argument("--max-exact", type=int, default=50000)
        sub.add_argument("--draws", type=int, default=4095)
    claim = commands.add_parser("claim-check")
    claim.add_argument("--methods", type=Path, required=True)
    claim.add_argument("--report", type=Path, required=True)
    claim.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "claim-check":
        from .claim_check import compare_methods
        from .schema import strict_read, write_json

        write_json(args.out, compare_methods(args.methods, strict_read(args.report)))
        return
    if args.command == "validate-manifest":
        for cohort in load_manifest(args.manifest):
            print(cohort.id, len(cohort.sections), "sections validated")
        return
    waves = "auto" if args.wavelengths == ["auto"] else tuple(map(float, args.wavelengths))
    options = {
        "wavelengths": waves,
        "max_exact": args.max_exact,
        "draws": args.draws,
        "alpha": args.alpha,
        "design": args.design,
        "stratum": args.stratum,
        "graph_mode": args.graph_mode,
        "nugget_mode": args.nugget_mode,
        "leave_one_out_tic": args.leave_one_out_tic,
    }
    if args.command == "quick-run":
        from .api import audit_anndata

        report = audit_anndata(
            args.h5ad,
            args.patient_col,
            args.group_col,
            section_col=args.section_col,
            slide_col=args.slide_col,
            batch_col=args.batch_col,
            coord_key=args.coord_key,
            layer=args.layer,
            cohort_id=args.cohort_id,
            covariate_key=args.covariate_key,
            library_size_col=args.library_size_col,
            **options,
        )
        save_report(report, args.out)
        print(report["id"], "audit saved to", args.out)
        return
    if args.command == "demo":
        from .fixtures import fixtures

        cohorts = list(fixtures())
        if args.write_fixtures:
            export_fixtures(cohorts, args.write_fixtures)
    else:
        cohorts = load_manifest(args.manifest, design=args.design, stratum=args.stratum)
    for cohort in cohorts:
        report = audit(cohort, **options)
        save_report(report, args.out)
        if (
            args.command == "demo"
            and waves == (250.0, 600.0, 1400.0)
            and args.alpha == 0.05
            and args.nugget_mode == "raw"
            and args.graph_mode == "grid"
            and args.design in (None, "independent")
            and not args.leave_one_out_tic
            and args.max_exact >= 924
        ):
            expected = {
                "synthetic-null": set(),
                "synthetic-pseudorep": {"pseudorep"},
                "synthetic-composition": {"composition", "nugget"},
                "synthetic-stage": {"design"},
            }
            observed = {k for k, v in report["dimensions"].items() if v["status"] != "green"}
            if observed != expected[cohort.id]:
                raise RuntimeError(f"Unexpected fixture alarm vector for {cohort.id}: {observed}")
        print(cohort.id, " ".join(f"{k}={v['status']}" for k, v in report["dimensions"].items()))


def export_fixtures(cohorts, root):
    import numpy as np

    from .schema import sha256, write_json

    specs = []
    for c in cohorts:
        folder = root / c.id
        folder.mkdir(parents=True, exist_ok=True)
        files = []
        for i, s in enumerate(c.sections):
            dest = folder / f"section-{i}.npz"
            data = {"coords": s.coords, "X": s.X, "features": np.array(c.features)}
            if s.C is not None:
                data["C"] = s.C
            if s.library_size is not None:
                data["library_size"] = s.library_size
            if s.domain is not None:
                data["domain"] = np.asarray(s.domain, str)
            np.savez_compressed(dest, **data)
            files.append(
                {
                    "patient_id": s.patient_id,
                    "group": s.group,
                    "slide_id": s.slide_id,
                    "batch": s.batch,
                    "stage": s.stage,
                    "section_id": s.section_id,
                    "path": dest.relative_to(root).as_posix(),
                    "sha256": sha256(dest),
                }
            )
        specs.append(
            {
                "id": c.id,
                "kind": "synthetic",
                "platform": c.platform,
                "source": c.source,
                "units": "um",
                "design": c.design,
                "stratum": c.stratum,
                **{
                    key: c.provenance[key]
                    for key in ("measurement", "layer", "selection", "fixture_seed", "construction")
                    if key in c.provenance
                },
                "files": files,
            }
        )
    write_json(root / "manifest.json", {"cohorts": specs})


if __name__ == "__main__":
    main()
