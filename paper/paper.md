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

`spatial-audit` ships checksummed archive adapters and measured audits spanning 20 distinct primary human cohorts across six platforms—including Visium Alzheimer's cortex [@chenAD], Xenium triple-negative breast cancer [@yanTNBC], CosMx severe COVID-19 lung [@geo327879], and MALDI-MSI pancreatic ductal adenocarcinoma [@msiPDAC]—together with four supplementary archive records and 77 independently verified R comparator fits. An independent factorial simulation of 10,080 design states (1,935,360 synthetic cohorts) maps conditional operating characteristics and transparently documents that under strong abundance-dependent heteroscedasticity with zero batch alignment, signed nugget adjustment still reaches up to 60.4% null rejection (28.1% on the prespecified display slice); accordingly, nugget flags are explicitly scoped as noise-sensitivity indicators rather than calibrated tests of latent biological structure. Software reliability is verified by an automated release gate encompassing 304 `pytest` unit, integration, and Playwright Chromium browser tests, `ruff`, `mypy`, strict JSON checks, offline DOI verification, and wheel/sdist installation tests.

# AI usage disclosure

Large language models and agentic coding assistants were used during software architecture drafting, test generation, documentation editing, and manuscript polishing. All statistical implementations, cohort adapters, numerical outputs, and bibliographic records were independently verified through automated deterministic test suites, cross-language R comparator reconstruction, checksummed archive receipts, and release-gate validation, with final responsibility retained by the author.

# References
