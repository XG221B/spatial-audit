# Statement of Significance — Journal of Biomedical Informatics (JBI)

## Problem or Issue
Multi-patient spatial transcriptomics and mass spectrometry imaging studies compare clinical groups using $10^4$–$10^6$ spatial measurements from small donor cohorts, making downstream biomarker discovery highly vulnerable to pseudoreplication, slide or stage aliasing, tissue-composition confounding, and abundance-linked short-range noise.

## What is Already Known
Existing bioinformatics packages target specific downstream estimands—such as pseudobulk differential expression, domain interactions, or within-section spatial variability—but do not provide a unified, modality-agnostic pre-inference audit of biological replication, acquisition identifiability, and physical-scale confounding.

## What This Paper Adds
We introduce *Spatial Audit*, an open-source five-dimensional audit engine, benchmarked across 20 human clinical cohorts from six platforms, 1,935,360 simulated cohorts, real-data empirical-null splits, and six original R comparators, paired with a zero-backend interactive web observatory and strict machine-readable provenance.
