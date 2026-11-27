"""Measure the operating-system process high-water mark of a full synthetic H5AD audit."""

from __future__ import annotations

import argparse
import os
import platform
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np

from spatial_audit.schema import strict_read, write_json

ROOT = Path(__file__).resolve().parent.parent


def make_files(folder, cells, genes, patients, density, seed):
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    rng = np.random.default_rng(seed)
    entries = []
    observed_nnz = 0
    for i in range(patients):
        n = cells // patients + int(i < cells % patients)
        draws = round(n * genes * density)
        raw = sparse.coo_matrix(
            (
                rng.integers(1, 5, size=draws, dtype=np.int16),
                (rng.integers(n, size=draws), rng.integers(genes, size=draws)),
            ),
            shape=(n, genes),
        ).tocsr()
        observed_nnz += raw.nnz
        coords = rng.uniform(0, 2000, size=(n, 2)).astype(np.float32)
        features = pd.DataFrame(index=[f"G{j}" for j in range(genes)])
        obj = ad.AnnData(X=raw, var=features)
        obj.obsm["spatial"] = coords
        obj.obsm["composition"] = np.sin(coords[:, [0]] / 500)
        name = f"p{i}.h5ad"
        obj.write_h5ad(folder / name)
        entries.append(
            {
                "path": name,
                "patient_id": f"P{i}",
                "group": "A" if i % 2 == 0 else "B",
                "batch": f"b{i // 2 % 2}",
                "slide_id": f"s{i // 2}",
                "stage": [float(i // 2), 0.0],
            }
        )
    spec = {
        "id": "synthetic-memory",
        "kind": "synthetic",
        "units": "um",
        "covariate_key": "composition",
        "platform": "Simulated sparse cell panel",
        "files": entries,
    }
    write_json(folder / "manifest.json", {"cohorts": [spec]})
    return observed_nnz


def worker(folder):
    from spatial_audit.report import audit, save_report
    from spatial_audit.schema import load_manifest

    started = perf_counter()
    (cohort,) = load_manifest(folder / "manifest.json")
    report = audit(cohort)
    save_report(report, folder / "output")
    elapsed = perf_counter() - started
    if sys.platform == "win32":
        import psutil

        peak = psutil.Process().memory_info().peak_wset
        counter = "Windows PeakWorkingSetSize"
    else:
        import resource

        hwm = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak = int(hwm if sys.platform == "darwin" else hwm * 1024)
        counter = "getrusage RUSAGE_SELF ru_maxrss"
    write_json(
        folder / "performance.json",
        {
            "wall_seconds": elapsed,
            "peak_rss_bytes": int(peak),
            "os_counter": counter,
            "dimensions": list(report["dimensions"]),
            "n_observations": report["n_observations"],
            "n_patients": report["n_patients"],
            "n_features": report["n_features"],
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cells", type=int, default=1_000_000)
    p.add_argument("--features", type=int, default=64)
    p.add_argument("--patients", type=int, default=10)
    p.add_argument("--density", type=float, default=0.03)
    p.add_argument("--seed", type=int, default=9401)
    p.add_argument("--out", type=Path, help="Default output depends on benchmark mode")
    p.add_argument(
        "--scaling-curve", action="store_true", help="Profile four sizes in fresh workers"
    )
    p.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.worker:
        worker(args.worker)
        return
    if (
        args.patients < 4
        or args.patients > args.cells
        or args.features < 1
        or not 0 < args.density <= 1
        or args.seed < 0
    ):
        p.error("invalid sparse benchmark dimensions")
    sizes = [10_000, 50_000, 250_000, 1_000_000] if args.scaling_curve else [args.cells]
    if args.patients > min(sizes):
        p.error("each size must include all patients")
    rows = [measure(args, cells) for cells in sizes]
    if args.scaling_curve:
        result = {"kind": "synthetic", "record_type": "scaling_curve", "rows": rows}
        name = "scaling_curve"
    else:
        result = rows[0]
        name = "memory"
    write_json(args.out or ROOT / f"results/synthetic/{name}.json", result)
    if not all(row["pass_rss"] for row in rows):
        raise SystemExit(1)


def measure(args, cells):
    with TemporaryDirectory(prefix="spatial-audit-memory-") as work:
        folder = Path(work)
        nnz = make_files(folder, cells, args.features, args.patients, args.density, args.seed)
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker", str(folder)],
            check=True,
            env=env,
        )
        result = strict_read(folder / "performance.json")
    limit = int(1.5 * 2**30)
    result.update(
        kind="synthetic",
        record_type="memory_benchmark",
        seed=args.seed,
        nnz=nnz,
        density=nnz / (cells * args.features),
        os=platform.system(),
        machine=platform.machine(),
        python=platform.python_version(),
        budget_bytes=limit,
        pass_rss=result["peak_rss_bytes"] < limit,
        scope="Fresh process, backed H5AD ingestion through serialization of all five "
        "dimensions. Fixture generation in parent process excluded. "
        "Sparse count panel; peak depends on geometry, features and sparsity.",
    )
    print(
        f"{result['n_observations']:,} observations, peak {result['peak_rss_bytes'] / 2**30:.3f} GiB"
    )
    return result


if __name__ == "__main__":
    main()
