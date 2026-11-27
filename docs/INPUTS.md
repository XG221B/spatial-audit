# Input contract and measured endpoints

A manifest is a UTF-8 JSON object with a nonempty `cohorts` array. Every cohort declares a portable ID, a provenance `kind` (`public`, `synthetic`, or `local`), `units: "um"`, platform and files. Public records include the authoritative source URL. Paths resolve relative to the manifest. SHA-256 is optional in a new input and is checked when supplied. Biological identity and calibration must come from the acquisition records. The default independent design assigns each patient to one group. Set `design` to `paired` for two conditions measured in every patient.

## Complete NPZ example

Generate reusable files with `spatial-audit demo --write-fixtures work/fixtures`. Each section NPZ contains a finite two-dimensional `X`, an aligned `coords` array with two columns, and a one-dimensional `features` array of Unicode strings. Optional `C` has one row per observation; optional `domain` is a per-row categorical factor for DESpace2. Pickled arrays are disabled. The accompanying generated manifest maps each file to patient, group and acquisition metadata.

## H5AD

Use one or multiple H5ADs. Names are aligned across files. A missing requested feature is an input error. Supply explicit feature names to reuse a targeted panel. Otherwise the first file sets the panel, optionally truncated by `feature_limit`. The following is a structural example; replace the file and accession URL with those of the intended public study, or use the archive preparation scripts for a runnable public manifest.

| Key | Meaning / default |
|---|---|
| files[].path | Required input path, one row per file |
| files[].patient_id, group, section_id | Optional constant labels for an entire file |
| columns | Mapping from canonical names patient_id, group, section_id, slide_id, batch to obs column names |
| layer | Input measurement layer; defaults to X |
| coord_key | obsm spatial coordinates, default spatial |
| covariates | Ordered list of numeric obs columns |
| covariate_key | Alternative obsm matrix, such as estimated cell-type fractions |
| domain | Categorical obs label for each spatial observation |
| stage_columns | Pair of obs columns containing the stage position, constant per acquisition section |
| files[].stage | Alternative pair of physical stage coordinates |
| preprocessing_evidence | Optional portable JSON filename for a companion archive reconstruction or input receipt |
| measurement | Declared scale such as counts or log_normalized, carried into the noise-sensitivity report |
| library_size_column | Full per-observation library depth when available, used by compatible comparators |
| selection | Short description of the selected patients, features and normalization |

Use one section per physical field with a coherent coordinate system. Convert pixels with documented instrument or image scale factors before constructing the manifest. Explicit stage coordinates must share the acquisition reference frame. A section with several donors can use per-observation patient labels. The reader checks within-section acquisition metadata and can stream selected feature blocks from dense, CSR and CSC HDF5 storage. The audit preserves the supplied measurement scale. Raw count, normalized abundance and intensity have different interpretations for both patient means and the noise model.

The GSE220442 `--visium-dir` importer generates `coordinate_calibration=classic_visium_100um`. Each legacy section uses its integer Space Ranger array row and column on the classic 6.5-mm capture area to recover 100 µm interspot spacing. For calibrated H5ADs without those addresses, the importer checks the observed pitch and reconstructs the hexagonal lattice within 0.5 µm before using the stored positions. This mode applies specifically to the classic Visium array. Image pixels, HD bins and Xenium centroids have separate coordinate conventions.

## Signals and thresholds

| Dimension | Measured quantity | Traffic-light rule |
|---|---|---|
| Power | Doubled smaller permutation tail, exact allocation count and first possible strict BH rank | Red if even the last rank cannot pass, yellow if the first possible rank exceeds one |
| Design | Patient-weighted categorical normalized mutual information, descriptive stage separation and AUC, and a sparse least-squares estimability check with all observed acquisition main effects | Red for exact fixed-effect aliasing or complete observed stage separation with more than eight independent patients. Geometric stage separation in smaller cohorts receives yellow in the absence of exact aliasing. Yellow for high mutual information (0.5), high in-sample AUC (0.85), missing acquisition fields or separated local boxes. The JSON distinguishes exact model aliasing from disjoint observed stage support. |
| Pseudoreplication | Pooled Mann–Whitney and BH versus equal-patient mean contrast and moderated permutation BH | Red when the pooled call fraction exceeds the patient fraction by 0.25 or more, yellow from 0.10; a zero denominator produces a null ratio |
| Composition | Within-band energy fraction explained by measured covariates | Red at 0.50, yellow at 0.25 or when missing data or spatial resolution limit interpretation |
| Nugget | Before/after residual fine-band contrast, sign of mean and estimated noise changes | A feature needs a material contrast, concordant noise/mean directions and at least half attenuation. Red when at least 20% of evaluable features flag; sparse evidence or individual flags give yellow |

The thresholds are screening conventions. The shipped audits use a strict nominal and BH threshold of 0.05. Use `--alpha` to select a different threshold. A red discovery-rate gap has no causal interpretation by itself. For all modalities, a measured linear separation can identify absent sampled support while the group column still has full rank in a linear adjustment. The output records both diagnostics. Missing external stage coordinates leave only local geometry for inspection. Power reports both the exhaustive floor and the attainable floor under the executed Monte Carlo budget, together with their first feasible strict BH ranks. The status uses the executed budget. A band requires connected two-dimensional spatial coverage, wavelength at least π times grid pitch and an estimable residual after excluding projections that exhaust effective band dimension. Per-band patient and section counts give its denominator. The nugget keeps signed energy differences; the measurement layer determines the physical interpretation. A log-normalized library yields sensitivity on its transformed scale. With the classic 100 µm Visium lattice the default 250 µm band is below the conservative resolution limit. An explicitly redeclared analysis such as `--wavelengths 400 800 1600` requests a different estimand; the per-section field-of-view and effective-dimension checks still apply.

## Outputs

Every cohort yields a versioned JSON result plus a Markdown report under its provenance directory. The five metric dictionaries include feature-level p/q values, nulls for undefined numbers, the exact or sampled permutation mode, per-patient residual band energies and input provenance. Basis ranks and noise gains are section specific. Repeated cohorts from one accession can illustrate feature-panel or donor-subset sensitivity; the dashboard and paper identify those variants. Public metadata, synthetic calibration and software release evidence have distinct paths.

## Single-file audits and optional designs

Run `spatial-audit validate-manifest --manifest manifest.json` to check input mappings, checksums and biological units. For a single calibrated H5AD, run `spatial-audit quick-run cohort.h5ad --patient-col donor --group-col condition --out reports`. Supply `--section-col` when a patient-condition contains several physical sections. Without a section column, the reader treats each patient-condition as one spatial field. Direct reports use the local provenance category. The [quickstart notebook](../examples/quickstart_tutorial.ipynb) exercises both entry points.

- Paired audits average sections within each patient-condition and jointly flip each patient's vector of feature differences. The exact doubled-tail floor with `n` complete pairs is `min(1, 2 / 2**n)`. Unmatched patients require an explicit cohort selection before this analysis.
- Stratified audits accept `design: "stratified"` and `stratum: "batch"` or `"slide_id"` per cohort, or the corresponding CLI options. Permutations preserve the arm counts within each stratum. Each patient must have one observed stratum. Single-arm strata contribute one allocation.
- The depth diagnostic divides each section's feature means by its mean library size. It uses total registered-panel abundance when a full depth vector is unavailable across the cohort. Patient summaries weight sections equally. The group depth ratio divides the second alphabetical group's mean patient depth by the first group's.
- The acquisition AUC uses one section-averaged centroid per patient-condition. Each held-out patient's scores come from a linear least-squares fit on the other patients. Categorical association uses the Pearson statistic with seeded patient-profile permutations for tables with more than two levels. The missing fraction counts sections. The missingness test compares patients with any missing section across groups, using paired discordance for paired designs.
- `--wavelengths` accepts one or more increasing physical wavelengths, or `auto` for log-spaced bands within the measured pitch and footprint limits. `--graph-mode voronoi` uses the convex hull as the reflecting boundary. Nugget estimates require matched collinear two-step paths. Irregular meshes without such paths have no corrected energy estimate.
- `--nugget-mode vst` applies the Anscombe transform to nonnegative volume means before computing both fine residual energy and its noise correction. Raw and stabilized contrasts have different units. Composition shares retain their input scale. `--leave-one-out-tic` treats the first supplied covariate as full total ion current in the same units as the ions. For each ion, the projection uses TIC minus that ion, retaining the other covariates.
