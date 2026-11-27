# Cover Letter — Journal of Biomedical Informatics (JBI)

Dear Editor-in-Chief and Editorial Board of the *Journal of Biomedical Informatics*,

I am pleased to submit the manuscript titled **"Spatial Audit: a generalizable five-dimensional statistical and experimental-design audit framework for multi-patient spatial omics"** for consideration as an Original Research / Methodological Review article in the *Journal of Biomedical Informatics* under the traditional subscription publishing model.

Clinical and translational spatial omics studies—spanning Visium, Xenium, Stereo-seq, CosMx SMI, MERFISH, and MALDI mass spectrometry imaging (MALDI-MSI)—increasingly profile tens of thousands to millions of spatial measurements across small human patient cohorts. While mature bioinformatics packages exist for downstream differential expression or within-section spatial variability testing, the biomedical informatics community currently lacks a unified, modality-agnostic pre-inference diagnostic framework to evaluate whether a multi-patient spatial dataset supports valid donor-level inference.

In this work, we introduce **`spatial-audit`**, a generalizable five-dimensional statistical and experimental-design audit framework and open-source software resource. It evaluates:
1. **Exact combinatorial and Monte Carlo donor-level permutation resolution** (supporting independent, paired, and batch/slide-stratified clinical designs);
2. **Acquisition-design collinearity** across slides, batches, and physical instrument stage coordinates (with patient-centroid separation and leave-one-patient-out cross-validated AUC);
3. **Pseudoreplication sensitivity** contrasting pooled observation-level tests against equal-patient moderated permutation tests and library-size normalized contrasts;
4. **Physical-scale covariate confounding** using finite-volume graph diffusion (square-grid and rigid-motion-invariant Voronoi), shared Chebyshev bandpass filtering, and within-band orthogonal projection (including leave-one-ion-out total ion current for MALDI-MSI); and
5. **High-frequency short-lag noise sensitivity**, including an empirical comparison of raw versus Anscombe variance-stabilized scales.

We demonstrate the framework's methodological generalizability across **20 distinct human clinical archives from six spatial platforms**, **77 independently verified R comparator fits across six established packages**, **real-data empirical-null donor splits** (showing 93.3%–95.1% false-positive BH calls under pooled cell testing versus 1.7%–2.7% under donor-level permutation), and **1,935,360 simulated cohorts** that transparently characterize calibration limits under abundance-linked heteroscedasticity. The resource runs end-to-end on 1,000,000 sparse observations in 34.5 s within 283.8 MiB peak memory and ships with a zero-backend offline interactive web observatory, strict RFC-8259 JSON provenance, a methods claim-checker, and local MCP agent tools.

I confirm that this manuscript is original, has not been published previously, and is not under consideration elsewhere. All code, checksummed provenance records, automated tests (333 passing tests across Windows, Linux, Docker, and Chromium), and the interactive observatory are publicly available under the MIT license at `https://github.com/XG221B/spatial-audit` and `https://xg221b.github.io/spatial-audit/`. The author declares no competing financial or non-financial interests and received no external funding.

Thank you for your time and consideration.

Sincerely,

**Chaohao Liu**  
Independent Researcher  
ORCID: `0009-0003-4341-2494`  
GitHub: `https://github.com/XG221B/spatial-audit`
