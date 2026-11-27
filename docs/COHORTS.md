# Public cohort evidence / 公开队列证据

目录纳入标准为可公开定位的空间数据产品和来自存档或论文的生物学分组。检索日期为 2026-11-27。研究论文中的叙述、方法与本项目实际计算的结果各有出处。每个登记号单独计数，跨研究仍可能存在患者重叠。主清单基于供者映射、物理坐标和独立组对比的可取得性选择，无法用来估计整个领域的设计问题发生率。

The catalog indexes the public spatial assay from each study. Statistical methods are linked to their original endpoint. Use the specified donor map and coordinate calibration before running a contrast. GEO RAW.tar may contain submitted processed outputs. Official raw reads can have additional access conditions. Asset links were transcribed from archive metadata. Local downloads receive their own checksum manifest.

| 目录 | 平台 | 存档与分组 | 原论文的空间观察 | 论文分析及定位 | 生物学单位与使用条件 |
|---|---|---|---|---|---|
| human | Visium | [GSE220442](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442) · Alzheimer disease vs control, middle temporal gyrus | The authors reported cortical-layer-specific expression differences in middle temporal gyrus and investigated selected markers by in situ hybridization. | Visium DE used Seurat FindMarkers for AD vs control. Validation compared RNA puncta in annotated cells using Mann–Whitney tests. Methods, Differential expression analysis and RNAscope quantification. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9773466/) | Preserve cortical layer and donor ID. Multiple spots within one cortical layer share a donor. |
| human | Xenium | [GSE348499](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348499) · Pathologic complete response vs residual disease after neoadjuvant chemotherapy for TNBC | Pretreatment cell-state communities and cancer-cell programs were associated with response. | Group-wise gene comparisons used Wilcoxon rank-sum tests with FDR adjustment. Cell-state analyses and follow-up survival associations used other endpoints. Results, Tumour archetypes and chemotherapy response; Methods. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13293894/) | Use pretreatment spatial donors. The single-cell and spatial cohorts have different sizes. |
| human | MALDI-MSI | [MTBLS12782](https://www.ebi.ac.uk/metabolights/MTBLS12782) · Pancreatic ductal adenocarcinoma with vs without pre-existing diabetes | Spatial metabolism and transcriptional measurements associated diabetes with cholesterol-related programs in pancreatic tumors. | The paper describes OPLS-DA, variable-importance and fold-change selection with permutation assessment of the model; Fig. 2 compares tumor-region metabolites. These are distinct analytical endpoints. Methods, Metabolomic data analysis; Fig. 2. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13176823/) | Select the imaging assay. Keep ion polarities, tissue regions and clinical cohorts distinct. An assay row is not necessarily a new patient. |
| human | Stereo-seq | [STT0000036](https://db.cngb.org/stomics/project/STT0000036/spatialGeneExpression) · Mismatch repair-proficient vs deficient colorectal cancer; response vs stable disease in the prespecified anti-PD-1 subset | In this study, spatial organization of the tumor–stroma boundary differed across mismatch-repair and treatment-response groups. | Leiden clustering identified tissue structures, and Moran’s I evaluated spatial organization. The paper compares specific patient and tissue-region strata. Fig. 1 and Supplementary Fig. 2; Methods, Stereo-seq data analysis. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11599708/) | Use the treated dMMR subset for a response comparison; treatment-naive pMMR vs dMMR is a separate contrast. The processed spatial archive includes a second tumor section for CRCP59. Resolve repeated sections before between-patient comparisons. The linked FASTQ collection has separate access conditions. |
| human | Visium | [GSE202011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202011) · Psoriasis severity and healthy controls; paired lesional and nonlesional skin | Spatial cell ecosystems stratified skin by disease severity, including at some nonlesional sites. | PseudobulkExpression produced sample- and cluster-level counts. DESeq2 collapseReplicates combined technical replicates before differential analysis. Methods, Pseudobulk differential gene expression analysis. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10502701/) | Use one predefined independent-donor contrast. Paired lesion status requires a paired design. |
| human | Visium | [GSE189184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189184) · Active ulcerative colitis vs checkpoint-inhibitor colitis vs healthy controls | Checkpoint-inhibitor colitis contained spatial immune microdomains with distinct chemotactic signatures and damaged epithelial regions. | Generalized linear models related cell-type decomposition scores and pathway AUCs to spatial clusters or disease condition. STAR Methods, Visium spatial transcriptomics analysis. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979251/) | Distinguish activity status and therapy, with sections nested within patients. The study includes several molecular assays. |
| human | Xenium | [GSE282123](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282123) · Colitis types and healthy controls; vedolizumab response in a selected cohort | Distinct inflammatory neighborhoods were associated with vedolizumab nonresponse; an epithelial neighborhood was associated with response. | Spatial pseudobulk DE used DESeq2. Xenium cell-frequency comparisons used Mann–Whitney or Kruskal–Wallis tests with FDR control according to comparison. Methods, Statistics and spatial transcriptomics pseudobulk analysis. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13367960/) | Resolve panel, donor, segmentation and dataset before forming a contrast. Check possible participant overlap with other colitis cohorts. |
| human | Visium | [GSE206391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206391) · Psoriasis vs atopic dermatitis vs lichen planus; lesional and matched nonlesional biopsies | Cytokine transcripts and nearby responder signatures formed localized patterns in inflamed skin. | Layer-wise spatial enrichment of cytokines was tested with Wilcoxon signed-rank tests; density-based clustering identified nearby responder transcripts. Methods, Cytokine quantification and spatial clustering. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9747967/) | Some capture areas contain several donors; duplicates and anatomical matches require the original mapping. |
| human | Xenium | [GSE294965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294965) · ANCA-associated GN, lupus nephritis, anti-GBM GN and reference kidney | The authors mapped progression from cellular crescents to sclerosis and localized PDGF- and TGFβ-related programs. | Cell-type marker comparisons used Wilcoxon tests and multiple-testing adjustment. Cell-signaling inference and mechanistic experiments address separate questions. Extended Data (cell-type markers); Methods, Cell–cell interactions. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12571875/) | Use the original spatial paper’s data-availability statement. Each slide can contain several biopsies. The additional multimodal slide must be checked for repeat material. |
| human | Xenium | [GSE336890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336890) · Acute interstitial nephritis vs acute tubular injury vs reference kidney | The investigators described inflammatory neighborhoods around injured proximal tubules with complement and NAMPT signals. | Proximal-tubule DE in AIN vs ATI used Wilcoxon with Bonferroni adjustment. Reference ST tissue values were presented descriptively for the specified cortex analysis. Proximal tubule results; Methods, Statistical analysis of cell densities and spatial neighbors. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13616907/) | Use the per-cell biopsy assignment and distinguish cortex from medulla. Reference tissue is a small descriptive group. |
| human | MERFISH | [GSE284005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284005) · Chronic active multiple sclerosis lesions vs healthy white matter; paired domains within lesions | CD8 T cells near inflamed foamy microglia were localized to chronic active lesion rims. | BANKSY defined spatial domains and SPATA2 modelled spatial gradients. NMF characterized T-cell neighborhoods; marker analyses used Seurat FindAllMarkers. Fig. 1 and Fig. 2; STAR Methods, Spatial gradient analysis. [Original study](https://pubmed.ncbi.nlm.nih.gov/41167189/) | The original MERFISH results explicitly include healthy control tissues. Match to supplementary donor identifiers: run-region names alone do not establish patient identity. |
| human | Xenium | [GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346) · Fibrotic lung vs unaffected donors; multiple regions per donor | Airspace remodeling was associated with a progression of epithelial and macrophage changes. | Logit-transformed proportions were regressed on pathology; limma-voom models of pathology-aggregated counts included TMA origin. Methods, Differential expression across pathologies and pseudotime. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906353/) | Retain donor, pathology, TMA and anatomical region identifiers. Subsamples of one lung share a donor. |
| human | Visium | [GSE190732](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190732) · Acute, chronic and prolonged COVID-19 vs non-COVID pneumonia controls | The authors identified CCL18- and CCL21-associated tissue niches in prolonged lung pathology. | Gene-set enrichment used fgsea on Wilcoxon-ranked genes. Figures comparing disease stages report the tests used in each panel. Methods, Gene set enrichment analysis; Fig. 1. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9922044/) | The control diagnosis is pneumonia. Interpret disease-duration categories with the timing of autopsy sampling. |
| human | Visium | [S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410) · Idiopathic pulmonary fibrosis vs unaffected human lung | Aberrant alveolar epithelial cells localized within fibrotic niches enriched in TGFβ-related signals. | The investigators pooled counts by donor for pseudobulk DE. DESeq2 also contrasted control alveolar regions with IPF alveolar or fibrotic regions. Fig. 1; Methods, Differential expression analysis; Data availability. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11319205/) | This accession is the public processed human subset. The associated mouse GEO study uses a different accession. |
| human | CosMx SMI | [GSE327879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327879) · Severe COVID-19 vs non-pulmonary control autopsies | The authors described different macrophage-centered niches in fibrotic and alveolar tissue. | BuildNicheAssay clustered local cell-type neighborhoods across both slides. FindMarkers DE for gene ontology used Bonferroni-adjusted p values. Methods §§2.6–2.8. [Original study](https://doi.org/10.70401/mc.2026.0014) | GEO states that cases and controls occupy separate slides. Preserve slide labels and specimen identity. |
| human | Visium | [GSE267226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267226) · Post-acute COVID lung fibrosis vs controls, human spatial subset | Human fibrotic regions showed spatial association between CD8 T-cell signatures and dysplastic epithelial repair. | The human spatial figures compare mapped cell-type gene signatures and enriched regions. Wilcoxon rank-sum tests are reported for expression comparisons in the broader study. Extended Data, human lung ST maps and mouse validation; Methods, Spatial transcriptomics. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12456998/) | Only human Visium belongs in this contrast; the experimental perturbations were performed in animals. |
| human | Visium | [10.5281/zenodo.6578047](https://zenodo.org/records/6578047) · Myocardial infarction vs control donor hearts; distinct regions of an infarcted heart | Spatial neighborhoods distinguished injury, repair and remodeling after ischemic damage. | One-sided Wilcoxon rank-sum tests compared estimated cell-type composition of a niche against other niches, with adjusted p values. Fig. 2 and Methods, Spatial transcriptomic profiling. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9364862/) | Select Visium H5AD assets and carry heart identity through infarct, border and remote regions. |
| human | Visium | [10.5281/zenodo.10012934](https://zenodo.org/records/10012934) · Idiopathic pulmonary fibrosis vs control, processed Visium subset | The investigators described fibrotic, macrophage and immune niches in lung tissue. | The study estimated cell-type abundance with cell2location, compared abundances using Pearson correlations and used Squidpy neighborhood enrichment. Fig. 2 and Methods, Spatial neighborhood enrichment. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/) | Combine technical sections only within donor. Visium and validation Xenium cover different sampling units. |
| human | Visium | [GSE318039](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318039) · Advanced COPD vs healthy donor lung | The study localized cytotoxic T-cell interactions with alveolar type 2 and goblet cells in different lung compartments. | The ligand and receptor expression comparisons in the respective cell populations used Wilcoxon tests with FDR correction. Figures on AT2 and goblet cell crosstalk, ligand and receptor expression panels. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13536628/) | GEO identifies each donor and acquisition batch. Disease status also captures the difference between donor tissue sources. |
| human | Visium | [GSE306288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306288) · Female vs male diseased human aortic valves in the Visium assay | The investigators reported localized COMP expression near calcification in male valves and increased SERPINE1 in fibrotic female valve tissue. | The spatial analysis used Moran’s I for spatial clustering and Wilcoxon rank-sum tests for sex-associated expression. Methods, spatial transcriptomics cell scoring and spatial distribution; Fig. 5. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13220916/) | Match by numbered patient, preserve A/B repeat sections, include only spatial GSMs with deposited tissue positions and restrict to diseased valves. |
| supplementary_donor_crosswalk | Visium | [GSE225475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225475) · Psoriatic lesions vs healthy skin | The investigators localized keratinocyte inflammatory programs and fibroblasts associated with immune-cell recruitment. | The spatial analyses included deconvolution, expression maps and module scores. Cluster marker genes within the six Visium datasets were tested with two-sided Wilcoxon rank-sum tests and BH correction. Main text and Supplementary Data 4 legend (MOESM3, Description of Additional Supplementary Files). [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10261041/) | Visium QC labels identify biopsies. Supplementary Data 1 and 2 do not link each Visium biopsy to an independent donor. Cluster-marker tests address a different estimand from a case-control donor contrast. |
| supplementary_longitudinal | Xenium | [GSE291246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291246) · Basal cell carcinoma before vs after PD-1 blockade; response and tissue-site strata | The paper observed spatial co-clustering of T- and B-cell pseudoclones and associations with treatment exposure. | Neighborhood enrichment assessed spatial co-occurrence. The authors correlated paired pseudoclone sizes across spatial and treatment strata. Results, Spatial transcriptomics of the BCR and TCR repertoires; Fig. 4. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12067265/) | A longitudinal biopsy, lymph node and skin section can belong to one donor. A pre/post analysis requires matched inference. |
| animal | Stereo-seq | [STDS0000056](https://db.cngb.org/stomics/datasets/STDS0000056/data) · Axolotl post-injury time points and uninjured tissue | The authors described injury-induced progenitor states at the wound site and a possible regenerative lineage. | RNA velocity and Monocle pseudotime were used to infer cell-state trajectories. These analyses estimate direction and ordering. Results, ReaEGCs and intermediate progenitors; Fig. 4. [Original study](https://pubmed.ncbi.nlm.nih.gov/36048929/) | Animal benchmark. Multiple rostrocaudal sections can come from one brain. Map each section to an animal before between-animal inference. |
| specimen | CosMx SMI | [GSE346380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346380) · Lung adenocarcinoma vs normal lung specimens | The archive provides expression, FOV and segmentation data for characterizing spatial heterogeneity; the series summary does not specify a tested spatial association. | No linked publication or named inferential test could be substantiated from this record. Series design describes imaging and segmentation. GEO Series summary and overall design. [GEO design record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346380) | The archive establishes specimen counts. Independent donor identities and pairing require additional documentation. |

## File access

The machine-readable catalog lists archive metadata and selected spatial product URLs. For GEO, the short SOFT endpoint carries accession, release status and the archive-supplied file list. MetaboLights indexes separate sample, assay and binary-imaging files. The two Zenodo entries provide release-specific object links and archive MD5 checksums. The STOmicsDB table for STT0000036 exposes its public-file flags and directories for individual tissue sections. A product without a listed direct download remains a catalog record. [Catalog](data/catalog.json)

CNSA labels CNP0002199 as Controlled. STT0000036 supplies the public human Stereo-seq entry. [CNSA access record](https://db.cngb.org/cnsa/ajax/project/public_view/?q=CNP0002199)

GSE346380 has an explicit specimen design. At the time of review, no matching primary paper with a specified inferential test could be linked. STDS0000056 uses an animal model. Both records appear in the supplement category.

## Reproduce the deposited MSI analysis

From the repository root with the Python package installed, choose a cache with sufficient space for the positive-mode binary spectra, XML indices and extracted sections. These commands use paths relative to the current directory.

```sh
python scripts/fetch_and_verify.py --accession MTBLS12782 --role metadata --out cache/MTBLS12782
python scripts/fetch_and_verify.py --accession MTBLS12782 --role data --out cache/MTBLS12782
python scripts/prepare_msi.py --cache cache/MTBLS12782 --manifest cache/msi.json --provenance results/public/MTBLS12782_inputs.json
spatial-audit run --manifest cache/msi.json --out results
```

The adapter selects equally spaced ranks of the unique m/z targets in the deposited metabolite assignment file. Within each pixel it sums centroid intensities inside each closed mass-tolerance window. Full-spectrum TIC and the two centered raster coordinates enter the nuisance matrix. The window panel, tolerance, calibration, group crosswalk, per-file hashes and exported section sizes are written by the adapter. The published assignment file reflects the original study’s feature selection. These windows measure signals around reference masses and carry no independent molecular identity confirmation.

报告中的分组标签来自官方 ISA 样本表。原表的 Normal 指同为胰腺癌且无糖尿病的组。读取器核对每组 imzML 与 IBD 的 UUID、嵌入式 MD5、像素间距和数组边界，另行计算 SHA-256。各切片内栅格坐标可作漂移代理，跨次运行的仪器坐标系仍需独立标定。TIC 包含目标离子窗自身的贡献，能量解释占比也会反映这部分数学相关性。

## Reproduce the treated colorectal Stereo-seq contrast

The adapter reads every deposited section's bin metadata. It retains the dMMR sections with archived anti-PD1 treatment and maps the clinical response categories using the deposited `group2` field. Section identifiers map to individual patients in this filtered contrast. The evidence includes all screening decisions, source URLs and SHA-256 digests. [STOmics archive](https://db.cngb.org/stomics/project/STT0000036/spatialGeneExpression) · [Archive file index](https://db.cngb.org/stsp/ajax/gene_expression_tableview/?q=STT0000036&page=1&per_page=30&expression_type=Gene%20Expression) · [Original article and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11599708/).

```sh
python scripts/prepare_stereo.py --fetch --cache cache/STT0000036 --manifest cache/stereo.json --counts-manifest cache/stereo-counts.json --provenance results/public/STT0000036_inputs.json
spatial-audit run --manifest cache/stereo.json --out results
```

With a complete existing cache, omit `--fetch`. The primary audit uses the deposited log-normalized expression of a predefined marker panel. The adapter reconstructs an optional integer layer by inverting the normalization and comparing its forward transformation against every matrix entry. Every bin must reproduce both its archived total and its detected-feature count before a cohort count manifest can be issued. The source provides local bin locations, domain labels and estimated cell-type composition. Treat the composition fractions as model-derived covariates that share measurements with the expression outcomes. [Machine-readable evidence](../results/public/STT0000036_inputs.json).

中文复现说明：样本筛选使用归档的错配修复、治疗和疗效字段。论文 Methods 提供 bin 的物理尺寸。主审计保留归档表达尺度；计数模型另用经过逐条正向重建校验的整数层。切片内部坐标只支持局部空间分析，跨切片仪器位置需要独立标定。患者级对比的精确检验分辨率可从审计输出中读取。


## Reproduce the annotated CosMx FOVs

The adapter matches each annotated cell to its GEO count row by the slide, FOV, and within-field cell ID. It checks the authors' composite identifier and independently fits the horizontal and vertical pixel pitches to the recorded field origins on each slide. Each FOV supplies a separate spatial graph. Cases and controls occupy separate slides, so differences between arms also include slide effects. Expression-derived cell classes provide nuisance indicators. Local coordinates support within-field distances. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327879) · [Input receipts and calibration](../results/public/GSE327879_inputs.json)

```sh
python scripts/fetch_and_verify.py --accession GSE327879 --role audit_input --out cache/GSE327879
python scripts/prepare_cosmx.py --cache cache/GSE327879 --manifest manifests/cosmx.json --provenance results/public/GSE327879_inputs.json
spatial-audit run --manifest manifests/cosmx.json --out results
```

## Reproduce the lesional psoriasis Visium contrast

Use the archived GEO `individual`, disease state, disease status and sequencing batch. The analysis retains the first deposited lesional section per PsO donor and one section per healthy donor. PsA and non-lesional skin are excluded from this comparison. All retained matrix barcodes must match in-tissue entries in their own deposited full spot map. The v1/v2 capture area has alternating columns in successive rows; the adapter converts those indices to within-array physical positions using the nominal pitch. Sequencing batch completely distinguishes the two arms. The archive supplies no verified shared instrument-stage registration for this audit. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202011) · [Selection and file receipts](../results/public/GSE202011_inputs.json)

```sh
python scripts/fetch_and_verify.py --accession GSE202011 --role audit_input --out cache/GSE202011
python scripts/prepare_visium.py --cache cache/GSE202011 --manifest manifests/psoriasis.json --provenance results/public/GSE202011_inputs.json
spatial-audit run --manifest manifests/psoriasis.json --out results
```

For a different original v1/v2 capture area with independently verified donor labels, supply a tab-separated mapping of `section_id`, `patient_id`, `group`, `matrix`, `positions` and optional `batch` and `slide_id`. Matrix and position file paths are relative to the TSV. This adapter uses the skin RNA feature panel specified in its source. Select a biologically appropriate panel before using it for a new study.

## Reproduce the MERFISH contrast

The adapter fetches the original MERSCOPE counts, centroids, region labels and published analysis-source archive, then verifies sample identities, cell intersections and section-level group assignments. Its evidence output records the tissue identities and domain and cell-type annotation inputs. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284005)

```sh
python scripts/prepare_merfish.py --cache cache/GSE284005 --manifest manifests/merfish.json --provenance results/public/GSE284005_inputs.json
spatial-audit run --manifest manifests/merfish.json --out results
```

## Reproduce the Kuppe myocardial-infarction subset

The contrast includes all unambiguously named pure fibrotic-zone sections and all control sections in the Zenodo record. Mixed-zone sections are outside this contrast. The adapter cross-checks each selection against the H5AD study title, disease ontology and the official file MD5. It retains prespecified marker genes present in every section. The evidence file lists candidate genes absent from the archived feature axes. [Zenodo study files](https://zenodo.org/records/6578047) · [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9364862/).

```sh
python scripts/prepare_kuppe.py --fetch --cache cache/Kuppe --manifest cache/kuppe.json --provenance results/public/zenodo-6578047-FZ_inputs.json
spatial-audit run --manifest cache/kuppe.json --out results
```

归档 H5AD 的空间嵌入按轴压缩到单位区间；一张切片还交换了坐标轴并翻转方向。适配器对每个条形码逐一匹配厂商标准版阵列布局，交叉验证显示变换，再用标称六角网格间距换算切片内物理坐标。厂商参考样本只贡献条形码地址。细胞类型权重来自原研究的去卷积估计，其解释占比依赖该模型。[Geometry evidence](../results/public/zenodo-6578047-FZ_inputs.json).

## Reproduce the GSE267226 post-acute lung fibrosis contrast

The deposited H5 matrices identify their capture chemistry as Visium V4. The spatial supplementary archive contains count matrices and tissue images. The adapter retrieves the manufacturer layout from two independent CytAssist datasets and verifies that both provide identical addresses for the entire barcode set. All study barcodes must map uniquely. Nominal array coordinates support local spatial statistics. Absolute acquisition positions and histological registration are unavailable. [GEO archive](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267226) · [Original article, Extended Data Fig. 1c–g](https://pmc.ncbi.nlm.nih.gov/articles/PMC12456998/) · [Manufacturer V4 glossary](https://www.10xgenomics.com/support/software/space-ranger/latest/getting-started/space-ranger-glossary).

```sh
python scripts/prepare_pasc.py --fetch --cache cache/GSE267226 --manifest cache/pasc.json --counts-manifest cache/pasc-counts.json --provenance results/public/GSE267226_inputs.json
spatial-audit run --manifest cache/pasc.json --out results
```

所有病例和对照均使用同一预设基因面板，只保留每份探针矩阵中都有且基因符号唯一的标志物。零文库深度的位置保留在条形码校验证据中，进入统计计算的 spot 均具有正的文库总量。主报告用全库按深度标准化后的对数值，计数模型可读取单独清单指向的整数层。此归档未提供逐 spot 的组成估计。[实测输入与验证](../results/public/GSE267226_inputs.json).


## Reproduce the prolonged COVID-19 Visium contrast

For each retained autopsy case, GEO supplies a filtered count matrix and disease group. The spatial transcriptomics Methods identify three donors per disease stage and three non-COVID pneumonia controls. This contrast uses prolonged disease and the pneumonia controls. The original HDF5 chemistry attribute and the deposited library protocol identify Visium v1. Two manufacturer example tissues have identical complete maps from spot barcode to the row and staggered column of a capture area. The adapter verifies that identity, checks every retained barcode and uses only array addresses from those references. Tissue masks, images and stage positions of the reference tissues have no role in the reconstruction. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190732) · [Study and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9922044/) · [Input hashes and crosswalk](../results/public/GSE190732_inputs.json)

```sh
python scripts/fetch_and_verify.py --accession GSE190732 --role audit_input --out cache/GSE190732
python scripts/prepare_covid_visium.py --cache cache/GSE190732 --manifest manifests/covid.json --provenance results/public/GSE190732_inputs.json
spatial-audit run --manifest manifests/covid.json --out results
```

## Reproduce the pulmonary fibrosis Visium contrasts

Both projects use the same prespecified lung marker list. Each adapter intersects this list with the study's unambiguous feature identifiers, keeps the archive's donor–section hierarchy and writes a separate manifest for integer-count tools. The BioStudies audit normalizes to the complete Space Ranger library. The Mayr audit retains archived log-normalized expression and compares each panel entry with a reconstruction from the complete deposited raw count matrix. Its evidence file records per-section numerical discrepancies. [Shared candidate panel](../data/lung_marker_panel.txt).

**Mayr et al.** The combined H5AD contains spots on direct-placement Visium and CytAssist slides, including adjacent sections from the same donor. The adapter checks donor, technology and section assignments against Supplementary Figure S1A and Table S2, using the deposited `patient`, `sampleID`, `technology`, `treatment` and `VisiumSlides` fields. Technical sections contribute within the same donor. Estimated cell-type composition uses normalized cell2location lower posterior abundance estimates. [Zenodo H5AD](https://zenodo.org/records/10012934) · [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/) · [Supplementary figures and tables](https://pmc-oa-opendata.s3.amazonaws.com/PMC11313858.1/sciadv.adl5473_sm.pdf).

```sh
python scripts/prepare_ipf.py --fetch --cache cache/Mayr --manifest cache/mayr.json --counts-manifest cache/mayr-counts.json --provenance results/public/zenodo-10012934_inputs.json
spatial-audit run --manifest cache/mayr.json --out results
```

**BioStudies S-BSST1410.** The public sample sheet assigns human sections to `subject_alias`, with specimen, slide, capture area and repeat fields. The Space Ranger ZIP provides filtered H5 counts, tissue positions and scale factors for each capture area. The deconvolution archive supplies spot-resolved abundance estimates using the published Habermann reference. Only uniquely matched, in-tissue, positive-library spots with deconvolution estimates enter this contrast. The input receipts record exclusions by section. [BioStudies record](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410) · [Official file index](https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BSST1410) · [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11319205/).

```sh
python scripts/prepare_biostudies_lung.py --fetch --cache cache/S-BSST1410 --manifest cache/biostudies.json --counts-manifest cache/biostudies-counts.json --provenance results/public/S-BSST1410_inputs.json
spatial-audit run --manifest cache/biostudies.json --out results
```

这些研究的真实样本数按存档供者编号确定。对同一供者内的多次取材，报告同时记录切片和组织块，使读者能检查组间取样差异。阵列坐标来自原始行列地址，并与图像像素距离交叉校验。Space Ranger 的 `spot_diameter_fullres` 使用理论图像直径进行换算。未经跨片登记的坐标只描述片内几何。[Space Ranger 缩放参数说明](https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs).


## Reproduce the pulmonary fibrosis Xenium contrast

The corrected GEO RDS contains the native Xenium nucleus centroids and raw RNA counts. The adapter matches its patient and disease-status fields to Supplementary Table 1, checks each core's TMA/run against the donor-level record, and chooses the lexicographically first core for each donor. The inventory preserves agreement flags and the selection for every core. Native spatial coordinates are in micrometres. The study's adjusted coordinates position cores for visualization. Fibrotic lungs span multiple diagnoses, so the clinical group contrast covers pulmonary fibrosis collectively. Cell-lineage indicators reflect the authors' molecular annotations. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346) · [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906353/) · [Acquisition and donor crosswalk](../results/public/GSE250346_inputs.json)

```sh
python scripts/fetch_and_verify.py --accession GSE250346 --role audit_input --out cache/GSE250346
python scripts/prepare_fibrosis.py --cache cache/GSE250346 --manifest manifests/fibrosis.json --provenance results/public/GSE250346_inputs.json
spatial-audit run --manifest manifests/fibrosis.json --out results
```

Install R, Matrix and SeuratObject for this adapter. Its default reader checks the corrected release SHA-256 and streams the original XDR counts by contiguous cell blocks. It projects only the registered targets, separately sums all RNA targets for each nucleus, and preserves the serialized dimnames and metadata for R to decode. Both dimension labels and cell order must agree. Every library total is cross-checked against the archived `nCount_RNA`. A private temporary directory on the cache volume holds the projected counts and metadata during reconstruction. `--native-rds` uses conventional full-object R deserialization when sufficient RAM is available. `--rscript` selects the local R executable.

## Reproduce the baseline colitis and lung/valve contrasts

The GSE282123 adapter joins the study's public annotated Xenium Dataset 1 to its raw count counterpart by acquisition and cell barcode. It checks the clinical supplement, keeps healthy controls and baseline pre-vedolizumab ulcerative colitis, and retains the patient membership of each tissue core and technical run. A shared treatment-response interpretation would require a different endpoint. Its input receipt records sources for the datasets and clinical table. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282123) · [Clinical joins and input files](../results/public/GSE282123_inputs.json)

```sh
python scripts/prepare_uc_xenium.py --cache cache/GSE282123 --manifest manifests/uc.json --provenance results/public/GSE282123_inputs.json
spatial-audit run --manifest manifests/uc.json --out results
python scripts/prepare_lung_valve.py --accession GSE318039 --cache cache/GSE318039 --manifest manifests/GSE318039.json --provenance results/public/GSE318039_inputs.json
spatial-audit run --manifest manifests/GSE318039.json --out results
python scripts/prepare_lung_valve.py --accession GSE306288 --cache cache/GSE306288 --manifest manifests/GSE306288.json --provenance results/public/GSE306288_inputs.json
spatial-audit run --manifest manifests/GSE306288.json --out results
```

The lung and valve adapter reads the official GEO sample descriptions and downloads the corresponding Visium files. The lung contrast compares advanced COPD with donor controls. The valve contrast compares female with male donors with calcific aortic valve disease and aggregates repeated sections within each patient. Each input receipt records the marker panel, donor mapping, matrix files and array geometry. [Lung receipt](../results/public/GSE318039_inputs.json) · [Valve receipt](../results/public/GSE306288_inputs.json)
