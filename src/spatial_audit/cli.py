"""spatial-audit command line."""

import argparse
from pathlib import Path

from .report import audit, save_report
from .schema import load_manifest


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spatial-audit")
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("demo", "run"):
        sub = commands.add_parser(command)
        if command == "run":
            sub.add_argument("--manifest", type=Path, required=True)
        else:
            sub.add_argument(
                "--write-fixtures", type=Path, help="Export reusable synthetic NPZs and manifest"
            )
        sub.add_argument("--out", type=Path, default=Path("results"))
        sub.add_argument("--wavelengths", type=float, nargs=3, default=[250.0, 600.0, 1400.0])
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
    if args.command == "demo":
        from .fixtures import fixtures

        cohorts = list(fixtures())
        if args.write_fixtures:
            export_fixtures(cohorts, args.write_fixtures)
    else:
        cohorts = load_manifest(args.manifest)
    for cohort in cohorts:
        report = audit(cohort, args.wavelengths, args.max_exact, args.draws)
        save_report(report, args.out)
        if (
            args.command == "demo"
            and args.wavelengths == [250.0, 600.0, 1400.0]
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
                "files": files,
            }
        )
    write_json(root / "manifest.json", {"cohorts": specs})


if __name__ == "__main__":
    main()
