"""Stream the corrected GSE250346 dense RNA assay out of its pinned XDR archive.

This accession-specific reader uses the byte offsets of three serialized REAL
matrices in the exact public release. Header, symbol, shape, count range and
release digest are validated. R still decodes the original metadata and
dimnames. A private reduced serialization retains the original attribute
pairlists and symbol references. Only the first matrix supplies measurements.
"""

from __future__ import annotations

import gzip
import shutil
import struct
from pathlib import Path

import numpy as np

from spatial_audit.schema import sha256

DIGEST = "f97b6e20365715fed5e96e831304f86e32e5dbefaeada6d58668a805fd4ebbdd"
# Byte offsets refer to the uncompressed, XDR-serialized RDS payload.
SLOTS = [(83, "counts"), (4519719135, "data"), (9039438166, "scale.data")]
N_TARGETS = 343
N_CELLS = 1630319
PANEL_SIZE = 32


def copy_bytes(source, destination, size):
    while size:
        data = source.read(min(4 * 1024**2, size))
        if not data:
            raise ValueError("Truncated R serialization")
        destination.write(data)
        size -= len(data)


def project(source_path: Path, directory: Path):
    if sha256(source_path) != DIGEST:
        raise ValueError("Streaming XDR projection requires the registered corrected GEO release")
    directory.mkdir(exist_ok=True, parents=True)
    reduced = directory / "metadata.rds"
    with (
        gzip.open(source_path, "rb") as source,
        reduced.open("wb") as dest,
        (directory / "counts.i32").open("wb") as counts,
        (directory / "depth.f64").open("wb") as depth,
    ):
        for offset, slot in SLOTS:
            symbol = struct.pack(">II", 0x40009, len(slot)) + slot.encode("ascii")
            copy_bytes(source, dest, offset - len(symbol) - source.tell())
            if source.read(len(symbol)) != symbol:
                raise ValueError(f"Unexpected R attribute symbol at {slot}")
            dest.write(symbol)
            header = source.read(8)
            if header != struct.pack(">II", 0x20E, N_TARGETS * N_CELLS):
                raise ValueError(f"Unexpected real matrix header at {slot}")
            # Emit a one-element list whose child is the exact original
            # attribute pairlist following this payload. It carries the
            # true dim, dimnames and all R serialization symbol references.
            dest.write(struct.pack(">II", 0x13, 1))
            if slot == "counts":
                for i in range(0, N_CELLS, 2048):
                    n = min(2048, N_CELLS - i)
                    block = np.frombuffer(source.read(n * N_TARGETS * 8), dtype=">f8").reshape(
                        n, N_TARGETS
                    )
                    if (
                        not np.isfinite(block).all()
                        or np.any(block < 0)
                        or np.any(block != np.floor(block))
                        or np.max(block) > np.iinfo(np.int32).max
                    ):
                        raise ValueError("Invalid archived RNA count value")
                    counts.write(block[:, :PANEL_SIZE].astype("<i4").tobytes())
                    depth.write(block.sum(axis=1, dtype=np.float64).astype("<f8").tobytes())
            else:
                source.seek(N_TARGETS * N_CELLS * 8, 1)
        shutil.copyfileobj(source, dest, 4 * 1024**2)
