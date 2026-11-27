---
title: 'spatial-audit: A five-dimensional statistical and experimental-design audit engine for multi-patient spatial omics'
tags:
  - Python
  - spatial transcriptomics
  - mass spectrometry imaging
  - reproducibility
  - pseudoreplication
  - confounding
  - bioinformatics
authors:
  - name: Chaohao Liu
    orcid: 0009-0003-4341-2494
    corresponding: true
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 27 November 2026
bibliography: paper.bib
---

# Summary

Modern spatial biology assays measure thousands of genes or metabolites at tens of thousands of locations across tissue sections, yet clinical and experimental cohorts typically include only a small number of patients. Without careful checks, analyses can mistake correlated cells for independent patients, confuse clinical group effects with slide or instrument batch effects, or attribute tissue-composition and short-range measurement noise to genuine spatial biology.

`spatial-audit` is an open-source Python package, command-line toolkit, local Model Context Protocol (MCP) server, and offline interactive web observatory that performs an automated five-dimensional statistical and experimental-design audit for multi-patient spatial transcriptomics (Visium, Xenium, Stereo-seq, CosMx SMI, MERFISH) and MALDI mass spectrometry imaging (MALDI-MSI) cohorts:

1. **Exact permutation resolution (`power`)**: Validates recorded donor identifiers and group labels (while true biological independence remains an experimental-design assumption), computes the exact two-sided permutation $p$-value floor $\min(1, 2/\binom{n_a+n_b}{n_a})$ and executed Monte Carlo floor, and reports the smallest feasible Benjamini--Hochberg (BH) rejection rank [@bh].
2. **Acquisition design & stage collinearity (`design`)**: Quantifies alignment between group labels and categorical acquisition factors (`slide_id`, `batch`) via patient-weighted normalized mutual information and binary Fisher exact tests, and tests for complete linear separation on physical stage coordinates.
3. **Pseudoreplication inflation (`pseudorep`)**: Contrasts pooled observation-level Mann--Whitney BH discovery rates against equal-patient moderated permutation testing with empirical Bayes variance shrinkage refit symmetrically under every label permutation [@limma; @pseudorep].
4. **Multiscale composition & covariate sharing (`composition`)**: Builds a no-flux finite-volume graph Laplacian on occupied micrometre control volumes, applies shared Chebyshev polynomial bandpass filters [@wavelets] across physical wavelengths, and computes the feature-wise mean of within-band orthogonally projected energy shares attributable to cell-type composition, anatomical regions, or total ion current (TIC).
5. **High-frequency short-lag noise sensitivity (`nugget`)**: Estimates short-range white-noise variance via a nonnegative matched-path two-lag semivariogram intercept scaled by the projected filter's white-noise gain, and evaluates how fine-scale patient contrasts shift after signed noise subtraction.

Each audit emits strict RFC-8259 JSON (representing undefined values as `null`) alongside a Markdown summary; strict JSON reports can be uploaded directly into the zero-backend browser observatory (`docs/index.html`) for local inspection without network transmission.

# Statement of need

Multi-patient spatial omics studies routinely profile $10^4$--$10^6$ spots, cells, or pixels from fewer than a dozen biological individuals. Treating correlated within-donor observations as independent replicates severely inflates nominal discovery rates [@pseudorep], while small donor counts impose hard combinatorial limits on exact permutation $p$-values. Simultaneously, slide or stage aliasing, macro-anatomical tissue composition, TIC self-inclusion, and abundance-linked short-lag noise can confound spatial group comparisons. Computational biologists, biostatisticians, and manuscript reviewers need a standardized, reproducible, modality-agnostic diagnostic report that surfaces these five failure modes before downstream biological interpretation.

# State of the field

Established packages such as `edgeR` [@edger], `DESeq2` [@deseq2], `DESpace2` [@despace2], `nnSVG` [@nnsvg], `SPARK-X` [@sparkx], and `Cardinal` [@cardinal] address distinct downstream inferential estimands: patient pseudobulk differential expression, domain-level spatial interactions, within-section spatial variability, or imaging mass spectrometry comparisons. None of these tools unify donor-level combinatorial resolution, acquisition-stage identifiability, pooled-versus-patient pseudoreplication auditing, micrometre-scale covariate orthogonalization, and short-lag nugget sensitivity across both transcriptomic and MSI modalities. Building `spatial-audit` as a dedicated, lightweight audit layer—rather than extending a single transcriptomic or R/Bioconductor ecosystem—enables consistent cross-platform diagnostics, backed HDF5/NPZ streaming, strict JSON provenance, and a zero-backend static observatory while still integrating original-package R comparators via an isolated subprocess benchmarking harness.

# Software design

`spatial-audit` is implemented in Python 3.13 (`numpy`, `scipy`, `pandas`, `anndata`, `statsmodels`) with portable relative paths and full Windows 11 and Linux compatibility. Core logic is separated into cohesive modules (`schema.py`, `design_audit.py`, `pseudorep.py`, `spectral_audit.py`, `claim_check.py`, `report.py`, `cli.py`, and `mcp_server.py`). Key architectural choices include: (i) backed HDF5/H5AD streaming so million-observation cohorts run within sub-gigabyte memory budgets (0.278 GiB peak RSS on 1,000,000 sparse observations across 10 donors and 64 features at 2.96% density); (ii) strict separation of public archive outputs and deterministic synthetic fixtures (`null`, `pseudorep`, `composition`, `stage`) with SHA-256 input receipts; (iii) a rule-based `claim-check` command for comparing audit outputs against methods text; and (iv) a local MCP stdio server (`spatial-audit-mcp`) exposing the same deterministic audit functions to agentic workflows.

# Research impact statement

`spatial-audit` ships checksummed archive adapters and measured audits spanning 20 distinct primary human cohorts across six platforms—including Visium Alzheimer's cortex [@chenAD], Xenium triple-negative breast cancer [@yanTNBC], CosMx severe COVID-19 lung [@geo327879], and MALDI-MSI pancreatic ductal adenocarcinoma [@msiPDAC]—together with four supplementary archive records and 77 native method–cohort runs with independently reconstructed BH calls. An independent factorial simulation of 10,080 design states (1,935,360 synthetic cohorts) maps conditional operating characteristics and transparently documents that under strong abundance-dependent heteroscedasticity with zero batch alignment, signed nugget adjustment still reaches up to 60.4% null rejection (28.1% on the prespecified display slice); accordingly, nugget flags are explicitly scoped as noise-sensitivity indicators rather than calibrated tests of latent biological structure. Software reliability is verified by an automated release gate encompassing `pytest` unit, integration, and Playwright Chromium browser tests, `ruff`, `mypy`, strict JSON checks, offline DOI verification, and wheel/sdist installation tests.

# Empirical validation

<!-- BEGIN measured validation -->

In GSE220442, expanding the panel for the same 6 cortical donors from 32 to the frozen 2,512-feature archived panel preserved all five screening categories. Pooled BH calls rose from 87.5% to 93.3%. With three donors per arm, the two-sided permutation floor of 0.10 guarantees zero donor-level BH calls at 0.05 on both panels. Cortical-layer covariates accounted for 35.9% and 33.9% of energy at 1,400 micrometres. Neither panel supported the 250-micrometre nugget endpoint.

For cortical SDF4, the pooled adjusted p was 3e-56 and the donor-adjusted p was 1. In the breast cancer panel, ABAT yielded pooled and donor-adjusted p values of 1.6e-30 and 0.992. For ABCB4, the pCR-minus-RD energy difference changed from -3.68e-05 to -7.98e-06 after short-lag subtraction. The ABCB4 contrast compares pretreatment spatial variation between eventual response groups. Without measured cell-type covariates in these sections, the contribution of cellular composition to this difference remains unresolved.

In pancreatic cancer, the without-minus-with-diabetes energy difference for m/z 178.00810 changed from 0.743 to 0.010. The archive identifies this feature by mass window. Leave-one-ion-out TIC modestly changed the joint TIC and raster projection from 45.09% to 44.88% at 1,400 micrometres. Both projections retain the raster coordinates and may capture gradients shared by tissue and acquisition.

The 35-donor lung fibrosis Xenium contrast GSE250346 yielded 9/32 donor-level BH calls, 13 edgeR calls and 16 DESeq2 calls. All three procedures called 9 shared features (ACKR1, AGER, APLNR, ATF6, AXIN2, BCL2L11, BMP4, BMPR2, CA4). AGER and CA4 link the panel to alveolar and vascular compartments, respectively. Whole-tissue summaries combine changes in compartment abundance and expression within cells. The shared calls represent cross-method concordance on the same donors.

In GSE348499, conditional randomization reused the observed donors and the same 32 features. Within RD (19 donors), donor-mean testing called 43/1,600 feature-splits (4/50 splits with calls) and pooled cells yielded 1,493/1,600 feature-splits (50/50 splits with calls). Within pCR (25 donors), donor-mean testing called 27/1,600 feature-splits (1/50 splits with calls) and pooled cells yielded 1,521/1,600 feature-splits (50/50 splits with calls).

With seven donors per arm and a twofold baseline increase under Poisson sampling, the raw test rejected in 10/192 cohorts (5.2%) and the Anscombe test in 74/192 (38.5%). The pointwise 95% Wilson interval for the latter was 31.9% to 45.6%. Across all 27 settings without batch alignment (10,368 simulated cohorts shared between modes), the raw maximum was 21/192 (10.9%). At dispersion 0.35 with the same donor count and baseline shift, alternative detection was 56.2% on the raw scale and 24.5% after transformation. A square-root transformation changes the amplitude of a fixed additive count pattern as baseline abundance changes. These single-feature simulations tested an absolute donor mean contrast. The production independent-donor test uses a moderated statistic across the registered feature panel.

Across 4 independently launched CPU audits with 10 donors and 64 sparse features, wall time increased from 2.7 to 34.5 seconds and peak RSS from 179.2 to 283.8 MiB between 10,000 and 1,000,000 observations. Each size has one host-specific measurement.

<!-- END measured validation -->

Paired contrasts support within-donor sign flips and stratified contrasts preserve group counts within each stratum. Optional clipped Voronoi geometry uses irregular sampling sites. The ion-specific TIC projection retains raster covariates while excluding each target ion from its own TIC. The targeted CosMx search found specimen or assay labels in GSE348495, GSE312757 and GSE311636. The inspected records did not establish the independent donor groups required for a new benchmark entry.

# AI usage disclosure

Large language models and agentic coding assistants were used during software architecture drafting, test generation, documentation editing, and manuscript polishing. All statistical implementations, cohort adapters, numerical outputs, and bibliographic records were independently verified through automated deterministic test suites, cross-language R comparator reconstruction, checksummed archive receipts, and release-gate validation, with final responsibility retained by the author.

# References
