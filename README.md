# spatial-audit

**多患者空间组学的可复现性与混杂敏感性审计**  
**Patient-level reproducibility and confounding audits for spatial omics**

打开 [离线看板](docs/index.html) 即可查看实测结果。网页支持本地 JSON、物理尺度曲线、采集位置、跨方法比较和统计陷阱沙盒。所有数据随静态页面分发，浏览器上传的文件在本机读取。

Open the dashboard directly from a clone, or deploy the same `docs/` directory to Pages. Public archive records, measured reports and synthetic controls have explicit provenance. A gray cell means the input has yet to be audited. Select **Synthetic fixtures** for the built-in demonstrations. [Dashboard](docs/index.html) · [Repository tree](docs/REPOSITORY_TREE.md) · [Study evidence](docs/COHORTS.md) · [Input contract](docs/INPUTS.md) · [R benchmark](docs/SOTA.md) · [中文答辩手册](docs/AUTHOR_GUIDE_AND_DEFENSE_CN.md) · [Refinement record](docs/SELF_CRITIQUE_AND_REFINEMENT.md)

## 三条命令复现 / Quick start

在 Python 3.13 环境中进入仓库根目录运行以下三条命令。第二条生成四个有标签的夹具及其数据清单，第三条通过标准输入路径重新审计。JSON 与 Markdown 保存在对应的 public 或 synthetic 子目录。Run from the repository root in a Python 3.13 environment.

```sh
python -m pip install -e ".[dev,agent,archives]"
spatial-audit demo --out work/demo --write-fixtures work/fixtures
spatial-audit run --manifest work/fixtures/manifest.json --out work/replayed
```

单独运行 `spatial-audit demo` 会执行固定种子的四个设计检查。严格零假设、伪重复、成分和均值混杂、载物台分离各有预先声明的灯号。随机队列的误报率与检出率在独立相图实验中估计。

## Public data / 公开队列与实测范围

The catalog covers twenty human spatial studies across Visium, Xenium, Stereo-seq, CosMx, MERFISH and MALDI-MSI. Four supplementary records retain provenance for the axolotl injury series STDS0000056, GSE346380 with unresolved donor independence, GSE225475 with an unverified biopsy-to-patient crosswalk, and GSE291246 with a longitudinal exposure contrast. GSE318039 and GSE306288 supply independent lung and valve contrasts in the primary benchmark. The evidence table links the accession, original spatial observation, statistical method and study-specific limitations. The dashboard shows which analyses have actually run. Panels and subsets from one study share patients; count accessions when describing independent archive coverage. [Catalog](docs/COHORTS.md) · [Machine-readable coverage](results/public/coverage.json)

GSE220442 的重建从官方计数、条码、空间比例因子和病人分组出发。32 特征面板用于跨方法比较，完整归档面板另行出报告。皮层分层指示量反映解剖组成。GSE348499 的化疗应答从 GEO 样本级元数据读取，跨方法比较子集使用每组按 GSM 排序的前三位供者。MTBLS12782 使用配对 imzML/ibd、实际像素间距、共同质量窗、全谱 TIC 和坐标漂移项。真实检查依赖每个存档的具体设计与可取得的数据层。

```sh
python scripts/prepare_gse220442.py --download --raw-dir work/ad-raw --out work/ad
spatial-audit run --manifest work/ad/manifest.json --out work/measured
python scripts/prepare_gse348499.py --download --verify-geo --raw-dir work/xenium-raw --out work/xenium
python scripts/fetch_and_verify.py --accession MTBLS12782 --role metadata --out work/msi
python scripts/fetch_and_verify.py --accession MTBLS12782 --role data --out work/msi
python scripts/prepare_msi.py --cache work/msi --manifest work/msi.json --provenance work/msi-inputs.json
python scripts/prepare_merfish.py --cache work/merfish --manifest work/merfish.json --provenance work/merfish-inputs.json
python scripts/prepare_stereo.py --fetch --cache work/stereo --manifest work/stereo.json --counts-manifest work/stereo-counts.json --provenance work/stereo-inputs.json
python scripts/prepare_cosmx.py --cache work/cosmx --manifest work/cosmx.json --provenance work/cosmx-inputs.json
python scripts/prepare_visium.py --cache work/psoriasis --manifest work/psoriasis.json --provenance work/psoriasis-inputs.json
python scripts/prepare_kuppe.py --fetch --cache work/heart --manifest work/heart.json --provenance work/heart-inputs.json
python scripts/prepare_colitis.py --cache work/colitis --manifest work/colitis.json --provenance work/colitis-inputs.json
python scripts/prepare_renal.py --cache work/kidney --manifest work/kidney.json --provenance work/kidney-inputs.json
```

Additional archive adapters use `--cache DIRECTORY --manifest FILE --provenance FILE`. Run the script with `--help` for its exact interface. For staged caches, pass the asset role listed for that study in the catalog to `fetch_and_verify.py --accession ACCESSION --role ROLE --out DIRECTORY`. The catalog includes study-specific acquisition and filtering instructions.

| Public contrast | Adapter under `scripts/` | Additional arguments |
|---|---|---|
| S-BSST1410, GSE267226, Zenodo 10012934 | `prepare_biostudies_lung.py`, `prepare_pasc.py`, `prepare_ipf.py`, respectively | `--fetch` for acquisition, `--counts-manifest FILE` for R |
| GSE336890, GSE190732, GSE206391 | `prepare_cortical_biopsies.py`, `prepare_covid_visium.py`, `prepare_skin.py`, respectively | Official assets in cache |
| GSE250346, GSE282123 | `prepare_fibrosis.py`, `prepare_uc_xenium.py`, respectively | The fibrosis R export requires Rscript |
| GSE318039, GSE306288 | `prepare_lung_valve.py` | `--accession GSE318039` or `--accession GSE306288` |

For existing prepared H5ADs, use a portable manifest or generate the compatible Visium/Xenium manifest with `scripts/run_benchmark_suite.py --visium-dir DIRECTORY --manifest FILE` or `--xenium-dir DIRECTORY --manifest FILE`. The latter command can take `--all-patients`. For an existing NPZ manifest the same runner accepts `--manifest FILE` directly. Input files are checksummed before use. To rebuild the checked-out dashboard from the results tree, run `python scripts/run_benchmark_suite.py`. Raw archive files remain in the chosen cache. See the provenance JSON in `results/public/` for actual byte counts, feature selection and biological mapping.

## Methods / 方法与解释

| Dimension | Output |
|---|---|
| Exact floor | Number of legal patient allocations, doubled-tail resolution and first feasible BH rank |
| Acquisition design | Patient-weighted categorical association, independent binary Fisher tests, stage support, rank and local-coordinate diagnostics |
| Pseudoreplication | Pooled rank-test BH calls, patient-level moderated permutation calls, rate gap and ratio |
| Multiscale composition | Finite-volume graph diffusion in µm, shared Chebyshev bands, within-band orthogonal projection |
| HF nugget | Two-lag short-range noise estimate, projected filter gain, raw and corrected patient-level fine energy |

Statistical flags describe the declared contrast, panel, coordinate scale, endpoint and assumptions. Discovery-rate differences can reflect changes in estimand and sampling weight. Biological signal, tissue boundaries and correlated instrument noise can contribute to short-lag variation. Exact enumeration assumes exchangeable independent patients. A reproducible Monte Carlo test is explicitly labeled when the assignment space exceeds the budget. See the exact implementation thresholds and output contract. [Methods](docs/INPUTS.md)

统一学术表述为“基于公开存档数据，在患者/生物学重复级去混杂框架下的可复现性、统计功效与混杂敏感性审计”。输入缺少可核验的患者映射、空间单位或对照定义时，保留证据缺口。固定夹具的命中情况仅适用于列明的构造。

## Methods claim check / 原文核对

```sh
spatial-audit claim-check --methods examples/GSE348499_excerpt.txt --report results/public/GSE348499-subset.json --out work/claim-matrix.json
```

The rule-based English extractor accepts UTF-8 text and PDFs with a text layer. It preserves a local quotation and page index for explicit numerical units and mentions of statistical strategies. Comparison with the audited subset includes a scope-review state. Review the original figure, cohort, negation and modality before interpreting a difference. [Subset requiring scope review](results/public/claim_check.json) · [Full-cohort numerical agreement](results/public/claim_check_full.json)

## AI agents

安装 agent extra 后，可让支持 MCP stdio 的客户端以本地用户身份启动 `spatial-audit-mcp`。配置示例适用于使用 `mcpServers` 键的客户端。Specify the installed executable from the desired environment if the program is outside the client's PATH. Use trusted local manifests.

```json
{
  "mcpServers": {
    "spatial-audit": { "command": "spatial-audit-mcp", "args": [] }
  }
}
```

| Tool | Arguments |
|---|---|
| audit_cohort | manifest, cohort_id |
| check_design_collinearity | manifest, cohort_id |
| compute_exact_p_floor | n0, n1, features, alpha |
| simulate_trap_cohort | trap = null, pseudorep, composition or stage (a string) |

## Benchmark, validation and manuscript

```sh
python scripts/run_phase_diagram.py
python scripts/measure_memory.py
python scripts/verify_dois.py
python scripts/run_benchmark_suite.py
python scripts/verify_release_gate.py
```

The phase diagram uses independent synthetic realizations at each specified design cell. The lookup returns the smallest *tested* donor count satisfying the reported uncertainty rule, or a blank when none qualifies. The memory measurement covers a fresh process from backed file loading through all five dimensions and serialization. Its empirical bound applies to the documented shape, sparsity and host. [Lookup CSV](docs/data/synthetic/sample_size_lookup.csv) · [Memory evidence](results/synthetic/memory.json)

The R benchmark needs R and its version-compatible Bioconductor libraries. The driver records successful runs, applicability, tested families, wall time, process-tree RSS and method versions. Its estimands include marginal group effects, within-section spatial variability and differential spatial patterns. [Install and replicate the six methods](docs/SOTA.md)

For browser QA install `.[browser]` and run `python -m playwright install chromium`. The full gate additionally needs pdfLaTeX and BibTeX on PATH. The repository contains a Resource-format editable manuscript, supplement, verified bibliography, generated macros and six vector PDF figures. Manuscript authorship and journal submission details are left to the submitting team. [Manuscript PDF](paper/main.pdf) · [Supplement PDF](paper/supplement.pdf) · [Registry checks](paper/generated/doi_verification.json)

## Deployment / 发布

- **GitHub Pages** — push the repository to the main branch and choose **Settings → Pages → GitHub Actions**. The supplied workflow publishes the static directory. Double-clicking `docs/index.html` also works offline.
- **Docker** — `docker build -t spatial-audit .` then `docker run --rm spatial-audit demo`. Mount an output directory and pass `--out` to keep generated reports. The Python image runs the audit and MCP dependencies. Install R packages using the separate comparator recipe.
- **PyPI** — configure a trusted publisher for this repository and the `pypi` GitHub environment. A published GitHub Release triggers testing, wheel/sdist builds and OIDC publication. Verify package-name availability and set the intended release version before publishing.
- **License** — software is distributed under MIT. Public data retain their archive terms, and original studies are cited in the evidence catalog. [MIT License](LICENSE)

