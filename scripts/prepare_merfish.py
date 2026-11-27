"""Build GSE284005 from the 17 official count/centroid/annotation triples."""

import argparse
import csv
import gzip
import os
import re
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from fetch_and_verify import fetch

from spatial_audit.schema import sha256, write_json

ACCESSION = "GSE284005"
SOURCE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + ACCESSION
CODE = "https://zenodo.org/api/records/17100316/files/HumanMS-main.zip/content"
MANUAL = "https://vizgen.com/wp-content/uploads/2024/05/91600001_MERSCOPE-Instrument-User-Guide_Rev-H.pdf"


def samples(soft):
    """Use the deposited sample titles without assigning patients from title prefixes."""
    rows = []
    for block in re.split(r"(?m)^\^SAMPLE = ", soft)[1:]:
        gsm = block.splitlines()[0].strip()
        title = re.search(r"(?m)^!Sample_title = (\w+)$", block)
        if not title or not re.fullmatch(r"GSM\d+", gsm):
            raise ValueError("Unexpected GEO record")
        title_text = title[1]
        values = dict(re.findall(r"(?m)^!Sample_characteristics_ch1 = ([^:]+): ([^\n]+)$", block))
        urls = re.findall(r"(?m)^!Sample_supplementary_file = (\S+)", block)
        assets = [
            {
                "accession": ACCESSION,
                "name": u.split("/")[-1],
                "url": u.replace("ftp://", "https://"),
            }
            for u in urls
            if u.startswith("ftp://ftp.ncbi.nlm.nih.gov/geo/")
        ]
        if len(assets) != 3 or not all(
            a["name"].startswith(gsm + "_" + title_text + "_") for a in assets
        ):
            raise ValueError("Expected a linked count, cell-class and coordinate file per GSM")
        rows.append((gsm, title_text, values, assets))
    if (
        len(rows) != 17
        or len({(r[2]["gender"], r[2]["age"], r[2]["ms duration"]) for r in rows}) != 17
    ):
        raise ValueError("Recheck donor/sample assignment against the published cohort")
    return rows


def project_section(count_path, position_path, annotation_path, panel):
    annotation = pd.read_csv(annotation_path, sep="\t", index_col="cells", keep_default_na=False)
    coords = pd.read_csv(position_path, sep="\t", index_col="cells")
    if not coords.index.is_unique or not annotation.index.is_unique:
        raise ValueError("Nonunique cell identifier")
    with gzip.open(count_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        cell_ids = next(reader)[1:]
        if len(set(cell_ids)) != len(cell_ids):
            raise ValueError("Duplicate count column")
        keep = [
            j
            for j, name in enumerate(cell_ids)
            if name in annotation.index and name in coords.index
        ]
        names = [cell_ids[i] for i in keep]
        if not names:
            raise ValueError("No aligned cells")
        features: dict[str, np.ndarray] = {}
        full_depth = np.zeros(len(names), dtype=float)
        seen = set()
        for row in reader:
            gene = row[0]
            if gene.lower().startswith("blank"):
                continue
            if gene in seen:
                raise ValueError("Duplicate feature")
            seen.add(gene)
            values = np.asarray(row[1:], float)[keep]
            if (
                not np.isfinite(values).all()
                or (values < 0).any()
                or (values != np.floor(values)).any()
            ):
                raise ValueError("Invalid full MERFISH count library")
            full_depth += values
            if (panel and gene in panel) or (not panel and len(features) < 32):
                features[gene] = values
    panel = panel or list(features)
    counts = np.column_stack([features[f] for f in panel])
    if np.any(counts < 0) or not np.all(counts == np.rint(counts)):
        raise ValueError("The selected assay is not integer-valued counts")
    cls = annotation.loc[names, "majorCluster_final"]
    indicators = pd.get_dummies(cls, dtype=np.float32)
    return {
        "X": counts.astype(np.float32),
        "features": np.array(panel),
        "coords": coords.loc[names, ["X", "Y"]].to_numpy(float),
        "C": indicators.to_numpy(),
        "domain": annotation.loc[names, "Region_banksy"].to_numpy(str),
    }, {
        "library_size": full_depth,
        "count_columns": len(cell_ids),
        "aligned": len(names),
        "covariates": list(indicators.columns),
        "regions": sorted(annotation["Region_banksy"].unique()),
    }


def prepare(cache, manifest, receipt):
    cache.mkdir(parents=True, exist_ok=True)

    def acquire(assets):
        for asset in assets:
            fetch({"files": [asset]}, cache, local=(cache / asset["name"]).is_file())

    acquire(
        [
            {
                "accession": ACCESSION,
                "name": "samples.soft",
                "url": SOURCE + "&targ=gsm&form=text&view=brief",
            },
            {"accession": ACCESSION, "name": "HumanMS-main.zip", "url": CODE},
        ]
    )
    rows = samples((cache / "samples.soft").read_text(encoding="utf-8"))
    with ZipFile(cache / "HumanMS-main.zip") as z:
        annotation_code = z.read("HumanMS-main/st_banksy.R").decode()
    healthy = set(
        re.findall(r'(?m)^(\w+)\.meta\$Region_banksy\s*=\s*"healthy_?WM"', annotation_code)
    )
    if len(healthy) != 3 or not healthy.issubset({r[1] for r in rows}):
        raise ValueError("Review author healthy-white-matter assignment")
    sections = []
    evidence = []
    panel: list[str] = []
    for gsm, title, meta, assets in rows:
        acquire(assets)
        prefix = cache / (gsm + "_" + title)
        data, info = project_section(
            str(prefix) + "_count.tsv.gz",
            str(prefix) + "_coordinates.tsv.gz",
            str(prefix) + "_celltypes.tsv.gz",
            panel,
        )
        panel = data["features"].tolist()
        group = "control" if title in healthy else "MS"
        if (float(meta["ms duration"]) == 0) != (group == "control"):
            raise ValueError("GEO and published sample keys require adjudication")
        dest = cache / (gsm + ".npz")
        np.savez_compressed(dest, **data)
        try:
            relative_path = os.path.relpath(dest, manifest.parent)
        except ValueError:
            relative_path = str(dest.resolve())
        depth_path = cache / (gsm + ".depth.npy")
        np.save(depth_path, info.pop("library_size"), allow_pickle=False)
        sections.append(
            {
                "path": relative_path,
                "library_size_path": (Path(relative_path).parent / depth_path.name).as_posix(),
                "library_size_sha256": sha256(depth_path),
                "patient_id": gsm,
                "group": group,
                "section_id": title,
                "sha256": sha256(dest),
            }
        )
        evidence.append({"sample": gsm, "title": title, "group": group, **info})
    spec = {
        "id": ACCESSION,
        "kind": "public",
        "platform": "MERFISH",
        "source": SOURCE,
        "units": "um",
        "measurement": "counts",
        "features": panel,
        "files": sections,
        "selection": "First 32 nonblank deposited genes in first GSM order, reused across all samples; all cells shared by counts, centroids and annotated major classes.",
    }
    write_json(manifest, {"cohorts": [spec]})
    from spatial_audit.schema import strict_read

    write_json(
        receipt,
        {
            "kind": "public",
            "record_type": "input_provenance",
            "accession": ACCESSION,
            "paper": "https://pubmed.ncbi.nlm.nih.gov/41167189/",
            "group_source": CODE,
            "coordinate_source": MANUAL,
            "assignment": "Published design comprises 14 MS and 3 control donors. The 17 GSMs have distinct demographic tuples. The authors' st_banksy.R identifies three healthy WM sections; all three have zero recorded MS duration. Treat each GSM as one donor for this contrast.",
            "covariates": "Within-section major cell-class indicators, derived from archived annotation. Projection may share molecular information with the tested panel.",
            "coordinate_frame": "Published MERSCOPE centroid X,Y; instrument global micrometres. Shared inter-slide stage registration unavailable.",
            "sections": evidence,
            "files": strict_read(cache / "dataset_manifest.json")["files"],
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.cache.resolve(), args.manifest.resolve(), args.provenance.resolve())
