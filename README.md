# spatial-audit

**Patient-level reproducibility, statistical resolution, and confounding audits for multi-patient spatial omics**  
**多患者空间组学与质谱成像队列的 5 维统计体检与混杂敏感性审计工具**

Open the zero-backend **[Interactive Dashboard](docs/index.html)** directly in any browser (works offline via `file://` or deployed to GitHub Pages) to inspect all measured public cohorts, cross-language R benchmarks, local JSON uploads, and the interactive statistical-trap sandbox.

[Dashboard](docs/index.html) · [Cohort Evidence](docs/COHORTS.md) · [Input Contract](docs/INPUTS.md) · [R Benchmark](docs/SOTA.md) · [Repository Tree](docs/REPOSITORY_TREE.md) · [中文一作答辩手册](docs/AUTHOR_GUIDE_AND_DEFENSE_CN.md) · [Validation & Critique](docs/SELF_CRITIQUE_AND_REFINEMENT.md)

---

## Overview / 项目简介

Spatial transcriptomics (Visium, Xenium, Stereo-seq, CosMx SMI, MERFISH) and MALDI-MSI generate tens of thousands of spatial observations from a limited number of donors. **`spatial-audit`** evaluates multi-patient group contrasts across five complementary statistical and experimental-design dimensions, outputting strict RFC-8259 JSON, Markdown audit summaries, and interactive visualizations:

| Dimension | What it measures | Key diagnostic output |
|---|---|---|
| **1. Power & Exact Floor** | True biological replicate count (patients/donors) and exact combinatorial resolution | Legal label allocations, doubled-tail exact $p$-value floor, executed Monte Carlo floor, first feasible BH rank |
| **2. Acquisition Design** | Alignment of clinical groups with `slide_id`, `batch`, and physical stage coordinates | Patient-weighted normalized mutual information, binary Fisher exact tests, linear stage separation, `UNIDENTIFIABLE_DESIGN` warnings |
| **3. Pseudoreplication** | Sensitivity of discovery rates to the unit of analysis | Pooled spot/pixel Mann–Whitney BH call rate vs. equal-patient moderated permutation BH call rate, rate gap, and ratio |
| **4. Multiscale Composition** | Fraction of spatial signal energy explained by tissue composition, covariates, or TIC at physical wavelengths ($\mu\mathrm{m}$) | Finite-volume graph Laplacian, shared Chebyshev bandpass filters, within-band orthogonal projection shares |
| **5. HF-Nugget Sensitivity** | Susceptibility of fine-scale spatial contrasts to mean-linked short-lag white noise | Matched-path two-lag semivariogram intercept, projected white-noise gain, signed raw vs. noise-adjusted fine-band energy |

> **Academic Framing / 学术措辞声明**  
> All outputs provide an objective methodological assessment of reproducibility, statistical resolution, and confounding sensitivity under a patient-level deconfounding framework on publicly archived data. Flags characterize the declared contrast, panel, coordinate calibration, and assumptions; they do not adjudicate the scientific validity of any published study.  
> 统一学术表述：“基于公开存档数据，在患者/生物学重复级去混杂框架下的可复现性、统计功效与混杂敏感性评估”。

---

## Quick Start (3 Commands) / 三条命令快速复现

Run from the repository root in a **Python 3.13** environment (Windows 11 / Linux compatible):

```sh
python -m pip install -e ".[dev,agent,archives]"
spatial-audit demo --out work/demo --write-fixtures work/fixtures
spatial-audit run --manifest work/fixtures/manifest.json --out work/replayed
```

1. Installs `spatial-audit` along with development, MCP agent, and archive-reading dependencies.
2. Runs the four built-in deterministic **Synthetic Benchmark Fixtures** (`null`, `pseudorep`, `composition`, `stage`) and exports them with a portable manifest.
3. Replays the exported manifest through the standard `spatial-audit run` CLI pipeline, generating strict RFC-8259 JSON and Markdown reports under `work/replayed/synthetic/`.

---

## Public Cohorts & Reproducibility / 公开队列与实测范围

The curated catalog covers **20 primary human spatial cohorts** across **6 platform families** (Visium, Xenium, Stereo-seq, CosMx SMI, MERFISH, and MALDI-MSI), yielding **22 public audit reports** (including explicitly labeled panel and donor subsets for cross-method comparison). Four supplementary records preserve provenance for cohorts outside the primary independent-donor benchmark: `STDS0000056` (axolotl injury model), `GSE346380` (unresolved donor independence), `GSE225475` (unverified biopsy-to-patient crosswalk), and `GSE291246` (longitudinal paired exposure contrast). Independent COPD (`GSE318039`) and aortic valve (`GSE306288`) contrasts are included in the primary suite.

Full provenance, donor crosswalks, coordinate calibrations, and SHA-256 receipts are documented in [docs/COHORTS.md](docs/COHORTS.md) and [results/public/coverage.json](results/public/coverage.json).

### Archive Preparation Commands

```sh
# Visium AD cortex (GSE220442) & Xenium TNBC chemotherapy response (GSE348499)
python scripts/prepare_gse220442.py --download --raw-dir work/ad-raw --out work/ad
spatial-audit run --manifest work/ad/manifest.json --out work/measured
python scripts/prepare_gse348499.py --download --verify-geo --raw-dir work/xenium-raw --out work/xenium

# MALDI-MSI PDAC (MTBLS12782)
python scripts/fetch_and_verify.py --accession MTBLS12782 --role metadata --out work/msi
python scripts/fetch_and_verify.py --accession MTBLS12782 --role data --out work/msi
python scripts/prepare_msi.py --cache work/msi --manifest work/msi.json --provenance work/msi-inputs.json

# Additional platforms (MERFISH, Stereo-seq, CosMx, Visium, Xenium)
python scripts/prepare_merfish.py --cache work/merfish --manifest work/merfish.json --provenance work/merfish-inputs.json
python scripts/prepare_stereo.py --fetch --cache work/stereo --manifest work/stereo.json --counts-manifest work/stereo-counts.json --provenance work/stereo-inputs.json
python scripts/prepare_cosmx.py --cache work/cosmx --manifest work/cosmx.json --provenance work/cosmx-inputs.json
python scripts/prepare_visium.py --cache work/psoriasis --manifest work/psoriasis.json --provenance work/psoriasis-inputs.json
python scripts/prepare_kuppe.py --fetch --cache work/heart --manifest work/heart.json --provenance work/heart-inputs.json
python scripts/prepare_colitis.py --cache work/colitis --manifest work/colitis.json --provenance work/colitis-inputs.json
python scripts/prepare_renal.py --cache work/kidney --manifest work/kidney.json --provenance work/kidney-inputs.json
```

Additional archive adapters under `scripts/` follow the same `--cache DIRECTORY --manifest FILE --provenance FILE` convention (pass `--help` for exact flags):

| Public Contrast | Adapter Script | Notes |
|---|---|---|
| `S-BSST1410`, `GSE267226`, `Zenodo 10012934` | `prepare_biostudies_lung.py`, `prepare_pasc.py`, `prepare_ipf.py` | Supports `--fetch` and `--counts-manifest FILE` for R count comparators |
| `GSE336890`, `GSE190732`, `GSE206391` | `prepare_cortical_biopsies.py`, `prepare_covid_visium.py`, `prepare_skin.py` | Uses official cached assets via `fetch_and_verify.py` |
| `GSE250346`, `GSE282123` | `prepare_fibrosis.py`, `prepare_uc_xenium.py` | `prepare_fibrosis.py` R export requires `Rscript` |
| `GSE318039`, `GSE306288` | `prepare_lung_valve.py` | Pass `--accession GSE318039` or `--accession GSE306288` |

For pre-prepared local H5AD or NPZ cohorts, use `scripts/run_benchmark_suite.py` (`--visium-dir`, `--xenium-dir`, or `--manifest`). Running `python scripts/run_benchmark_suite.py` with no arguments rebuilds the offline dashboard bundles from `results/`.

---

## Methods Claim Check / 论文方法学文本核对

Extract explicit unit counts and statistical strategy mentions from UTF-8 text or text-layer PDFs and compare them against a measured audit report:

```sh
spatial-audit claim-check --methods examples/GSE348499_excerpt.txt --report results/public/GSE348499-subset.json --out work/claim-matrix.json
```

Each match preserves its local quotation and page index, flagging scope differences for human review. See [subset scope review](results/public/claim_check.json) and [full-cohort agreement](results/public/claim_check_full.json).

---

## AI Agent Integration (MCP)

With the `agent` extra installed, `spatial-audit-mcp` exposes the audit engine over local **MCP stdio** for any compatible AI client:

```json
{
  "mcpServers": {
    "spatial-audit": { "command": "spatial-audit-mcp", "args": [] }
  }
}
```

| MCP Tool | Arguments | Purpose |
|---|---|---|
| `audit_cohort` | `manifest`, `cohort_id` | Run the full 5-dimension audit on a trusted local manifest |
| `check_design_collinearity` | `manifest`, `cohort_id` | Inspect acquisition-stage separation and categorical aliasing |
| `compute_exact_p_floor` | `n0`, `n1`, `features`, `alpha` | Compute exact permutation resolution and minimum feasible BH rank |
| `simulate_trap_cohort` | `trap` (`null`, `pseudorep`, `composition`, or `stage`) | Generate and audit a deterministic synthetic control |

---

## Benchmarks, Validation & Manuscript / 基准验证与论文构建

```sh
python scripts/run_phase_diagram.py     # 10,080 design states / 1,935,360 synthetic cohorts
python scripts/measure_memory.py        # 1M-observation sparse CPU memory profile (0.278 GiB peak RSS)
python scripts/verify_dois.py           # Bibliographic registry verification (33 references)
python scripts/run_benchmark_suite.py   # Rebuild docs/data/ JSON and offline JS bundles
python scripts/verify_release_gate.py   # Full release gate: pytest, ruff, mypy, JSON, LaTeX, wheel/sdist
```

- **Cross-language R benchmark**: Runs original `edgeR`, `DESeq2`, `DESpace2`, `nnSVG`, `SPARK-X`, and `Cardinal` implementations in isolated R subprocesses with independent BH reconstruction from raw $p$-values ([docs/SOTA.md](docs/SOTA.md)).
- **Browser QA**: Install `.[browser]` and run `python -m playwright install chromium` before running `pytest` to execute offline Chromium checks on `docs/index.html`.
- **Manuscript**: Requires `pdflatex` and `bibtex` on `PATH` to compile [paper/main.pdf](paper/main.pdf) and [paper/supplement.pdf](paper/supplement.pdf) with auto-generated metrics and vector figures.

---

## Deployment & License / 部署与许可

- **GitHub Pages**: Push to `main` and enable **Settings → Pages → GitHub Actions**. The workflow in `.github/workflows/pages.yml` publishes `docs/` directly.
- **Docker**: Build with `docker build -t spatial-audit .` and run `docker run --rm spatial-audit demo` (mount a volume with `--out` to persist reports).
- **PyPI**: Configure a Trusted Publisher for the `pypi` GitHub environment; publishing a GitHub Release triggers `.github/workflows/publish.yml`.
- **License**: Software is released under the [MIT License](LICENSE). Public datasets remain subject to their original archive terms and citations.
