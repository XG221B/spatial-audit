# Cover Letter — Journal of Biomedical Informatics (JBI)

Dear Editor-in-Chief and Editorial Board of the *Journal of Biomedical Informatics*,

I am pleased to submit the manuscript titled **"Spatial Audit: a generalizable five-dimensional statistical and experimental-design audit framework for multi-patient spatial omics"** for consideration as an **Original Research** article in the *Journal of Biomedical Informatics* under the traditional subscription publishing model.

Clinical and translational spatial omics studies—spanning Visium, Xenium, Stereo-seq, CosMx SMI, MERFISH, and MALDI mass spectrometry imaging (MALDI-MSI)—increasingly profile tens of thousands to millions of spatial measurements across small human patient cohorts. While mature bioinformatics packages exist for downstream differential expression or within-section spatial variability testing, the biomedical informatics community currently lacks a unified, modality-agnostic pre-inference diagnostic framework to evaluate whether a multi-patient spatial dataset supports valid donor-level inference.

In this work, we introduce **`spatial-audit`**, a generalizable five-dimensional statistical and experimental-design audit framework and open-source software resource. It evaluates:
1. **Exact combinatorial and Monte Carlo donor-level permutation resolution** (supporting independent, paired, and batch/slide-stratified clinical designs);
2. **Acquisition-design collinearity** across slides, batches, and physical instrument stage coordinates (with patient-centroid separation and leave-one-patient-out cross-validated AUC);
3. **Pseudoreplication sensitivity** contrasting pooled observation-level tests against equal-patient moderated permutation tests and library-size normalized contrasts;
4. **Physical-scale covariate confounding** using finite-volume graph diffusion (square-grid and rigid-motion-invariant Voronoi), shared Chebyshev bandpass filtering, and within-band orthogonal projection (including leave-one-ion-out total ion current for MALDI-MSI); and
5. **High-frequency short-lag noise sensitivity**, including an empirical comparison of raw versus Anscombe variance-stabilized scales.

We evaluated 20 primary human archives across six platforms. The R comparison comprised 77 native method–cohort runs with independently reconstructed BH calls. Eleven archives supported at least one covariate projection and nine supported evaluation of short-range noise for at least one feature.

In the breast cancer empirical-null experiment, labels were randomized within response groups, conditional on the observed donors and 32 features. In residual disease, donor testing made 43/1,600 feature-split calls (4/50 splits with calls), compared with 1,493/1,600 for pooled cells (50/50 splits). In pathological complete response, the corresponding numbers were 27/1,600 (1/50 splits) and 1,521/1,600 (50/50 splits). These two denominators show how a small feature-level call fraction can still affect some donor splits.

I confirm that this manuscript is original, has not been published previously, and is not under consideration elsewhere. The Windows release gate runs 333+ tests, including Chromium browser checks. CI is configured for Linux and Windows and includes a separate Docker demo smoke test. The MIT-licensed code, provenance records and interactive observatory are available at `https://github.com/XG221B/spatial-audit` and `https://xg221b.github.io/spatial-audit/`. The author declares no competing financial or non-financial interests and received no external funding.

Thank you for your time and consideration.

Sincerely,

**Chaohao Liu**  
Independent Researcher  
ORCID: `0009-0003-4341-2494`  
GitHub: `https://github.com/XG221B/spatial-audit`
