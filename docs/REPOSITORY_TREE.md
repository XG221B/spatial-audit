# Repository tree / 仓库目录

Scientific outputs are included below. The ignored local distribution directory holds the wheel and source tarball built by `python -m build`. The Python package configuration is in [pyproject.toml](../pyproject.toml).

```text
spatial-audit/
├── .github/
│   └── workflows/
│       ├── ci.yml
│       ├── pages.yml
│       └── publish.yml
├── data/
│   ├── GSE220442/
│   │   ├── archive_fingerprints.json
│   │   ├── feature_panel.tsv
│   │   └── README.md
│   ├── GSE348499/
│   │   ├── feature_panel.txt
│   │   ├── geo_samples.json
│   │   └── README.md
│   └── lung_marker_panel.txt
├── docs/
│   ├── data/
│   │   ├── public/
│   │   │   ├── bundle.js
│   │   │   ├── GSE189184.json
│   │   │   ├── GSE189184_inputs.json
│   │   │   ├── GSE190732.json
│   │   │   ├── GSE190732_inputs.json
│   │   │   ├── GSE202011.json
│   │   │   ├── GSE202011_inputs.json
│   │   │   ├── GSE206391.json
│   │   │   ├── GSE206391_inputs.json
│   │   │   ├── GSE220442-full-panel.json
│   │   │   ├── GSE220442.json
│   │   │   ├── GSE220442_reconstruction.json
│   │   │   ├── GSE250346.json
│   │   │   ├── GSE250346_inputs.json
│   │   │   ├── GSE267226.json
│   │   │   ├── GSE267226_inputs.json
│   │   │   ├── GSE282123.json
│   │   │   ├── GSE282123_inputs.json
│   │   │   ├── GSE284005.json
│   │   │   ├── GSE284005_inputs.json
│   │   │   ├── GSE294965.json
│   │   │   ├── GSE294965_inputs.json
│   │   │   ├── GSE306288.json
│   │   │   ├── GSE306288_inputs.json
│   │   │   ├── GSE318039.json
│   │   │   ├── GSE318039_inputs.json
│   │   │   ├── GSE327879.json
│   │   │   ├── GSE327879_inputs.json
│   │   │   ├── GSE336890.json
│   │   │   ├── GSE336890_inputs.json
│   │   │   ├── GSE348499-subset.json
│   │   │   ├── GSE348499.json
│   │   │   ├── GSE348499_reconstruction.json
│   │   │   ├── MTBLS12782.json
│   │   │   ├── MTBLS12782_inputs.json
│   │   │   ├── S-BSST1410.json
│   │   │   ├── S-BSST1410_inputs.json
│   │   │   ├── STT0000036.json
│   │   │   ├── STT0000036_inputs.json
│   │   │   ├── zenodo-10012934.json
│   │   │   ├── zenodo-10012934_inputs.json
│   │   │   ├── zenodo-6578047-FZ.json
│   │   │   └── zenodo-6578047-FZ_inputs.json
│   │   ├── synthetic/
│   │   │   ├── bundle.js
│   │   │   ├── phase_diagram.svg
│   │   │   ├── sample_size_lookup.csv
│   │   │   ├── synthetic-composition.json
│   │   │   ├── synthetic-null.json
│   │   │   ├── synthetic-pseudorep.json
│   │   │   └── synthetic-stage.json
│   │   ├── audit_summary.json
│   │   ├── catalog.json
│   │   └── sota_registry.json
│   ├── app.js
│   ├── AUTHOR_GUIDE_AND_DEFENSE_CN.md
│   ├── COHORTS.md
│   ├── index.html
│   ├── INPUTS.md
│   ├── release_gate.json
│   ├── REPOSITORY_TREE.md
│   ├── SELF_CRITIQUE_AND_REFINEMENT.md
│   ├── SOTA.md
│   └── styles.css
├── examples/
│   └── GSE348499_excerpt.txt
├── paper/
│   ├── figures/
│   │   ├── benchmark.pdf
│   │   ├── composition.pdf
│   │   ├── exact-floor.pdf
│   │   ├── nugget.pdf
│   │   ├── phase-diagram.pdf
│   │   └── replication.pdf
│   ├── generated/
│   │   ├── cohorts.tex
│   │   ├── doi_verification.json
│   │   └── metrics.tex
│   ├── main.pdf
│   ├── main.tex
│   ├── references.bib
│   ├── supplement.pdf
│   └── supplement.tex
├── results/
│   ├── public/
│   │   ├── claim_check.json
│   │   ├── claim_check_full.json
│   │   ├── coverage.json
│   │   ├── dataset_manifest.json
│   │   ├── GSE189184.json
│   │   ├── GSE189184.md
│   │   ├── GSE189184_inputs.json
│   │   ├── GSE190732.json
│   │   ├── GSE190732.md
│   │   ├── GSE190732_inputs.json
│   │   ├── GSE190732_readiness.json
│   │   ├── GSE202011.json
│   │   ├── GSE202011.md
│   │   ├── GSE202011_inputs.json
│   │   ├── GSE206391.json
│   │   ├── GSE206391.md
│   │   ├── GSE206391_inputs.json
│   │   ├── GSE220442-full-panel.json
│   │   ├── GSE220442-full-panel.md
│   │   ├── GSE220442.json
│   │   ├── GSE220442.md
│   │   ├── GSE220442_reconstruction.json
│   │   ├── GSE225475_readiness.json
│   │   ├── GSE250346.json
│   │   ├── GSE250346.md
│   │   ├── GSE250346_inputs.json
│   │   ├── GSE267226.json
│   │   ├── GSE267226.md
│   │   ├── GSE267226_inputs.json
│   │   ├── GSE267226_readiness.json
│   │   ├── GSE282123.json
│   │   ├── GSE282123.md
│   │   ├── GSE282123_inputs.json
│   │   ├── GSE284005.json
│   │   ├── GSE284005.md
│   │   ├── GSE284005_inputs.json
│   │   ├── GSE291246_readiness.json
│   │   ├── GSE294965.json
│   │   ├── GSE294965.md
│   │   ├── GSE294965_inputs.json
│   │   ├── GSE306288.json
│   │   ├── GSE306288.md
│   │   ├── GSE306288_inputs.json
│   │   ├── GSE318039.json
│   │   ├── GSE318039.md
│   │   ├── GSE318039_inputs.json
│   │   ├── GSE327879.json
│   │   ├── GSE327879.md
│   │   ├── GSE327879_inputs.json
│   │   ├── GSE336890.json
│   │   ├── GSE336890.md
│   │   ├── GSE336890_inputs.json
│   │   ├── GSE348499-subset.json
│   │   ├── GSE348499-subset.md
│   │   ├── GSE348499.json
│   │   ├── GSE348499.md
│   │   ├── GSE348499_reconstruction.json
│   │   ├── MTBLS12782.json
│   │   ├── MTBLS12782.md
│   │   ├── MTBLS12782_inputs.json
│   │   ├── S-BSST1410.json
│   │   ├── S-BSST1410.md
│   │   ├── S-BSST1410_inputs.json
│   │   ├── sota.json
│   │   ├── sota.md
│   │   ├── STT0000036.json
│   │   ├── STT0000036.md
│   │   ├── STT0000036_inputs.json
│   │   ├── STT0000036_refinement.json
│   │   ├── zenodo-10012934.json
│   │   ├── zenodo-10012934.md
│   │   ├── zenodo-10012934_inputs.json
│   │   ├── zenodo-6578047-FZ.json
│   │   ├── zenodo-6578047-FZ.md
│   │   └── zenodo-6578047-FZ_inputs.json
│   └── synthetic/
│       ├── memory.json
│       ├── nugget_refinement.json
│       ├── phase_diagram.json
│       ├── synthetic-composition.json
│       ├── synthetic-composition.md
│       ├── synthetic-null.json
│       ├── synthetic-null.md
│       ├── synthetic-pseudorep.json
│       ├── synthetic-pseudorep.md
│       ├── synthetic-stage.json
│       └── synthetic-stage.md
├── scripts/
│   ├── build_paper.py
│   ├── export_fibrosis_panel.R
│   ├── fetch_and_verify.py
│   ├── measure_memory.py
│   ├── prepare_biostudies_lung.py
│   ├── prepare_colitis.py
│   ├── prepare_cortical_biopsies.py
│   ├── prepare_cosmx.py
│   ├── prepare_covid_visium.py
│   ├── prepare_fibrosis.py
│   ├── prepare_gse220442.py
│   ├── prepare_gse348499.py
│   ├── prepare_ipf.py
│   ├── prepare_kuppe.py
│   ├── prepare_lung_valve.py
│   ├── prepare_merfish.py
│   ├── prepare_msi.py
│   ├── prepare_pasc.py
│   ├── prepare_renal.py
│   ├── prepare_skin.py
│   ├── prepare_stereo.py
│   ├── prepare_uc_xenium.py
│   ├── prepare_visium.py
│   ├── project_fibrosis_rds.py
│   ├── run_benchmark_suite.py
│   ├── run_phase_diagram.py
│   ├── run_sota_benchmark.py
│   ├── run_sota_benchmark.R
│   ├── verify_dois.py
│   └── verify_release_gate.py
├── src/
│   └── spatial_audit/
│       ├── __init__.py
│       ├── archive.py
│       ├── claim_check.py
│       ├── cli.py
│       ├── design_audit.py
│       ├── fixtures.py
│       ├── h5_backed.py
│       ├── mcp_server.py
│       ├── pseudorep.py
│       ├── report.py
│       ├── schema.py
│       ├── spectral_audit.py
│       └── visium_geometry.py
├── tests/
│   ├── test_biopsy_archives.py
│   ├── test_cosmx.py
│   ├── test_dashboard.py
│   ├── test_engine.py
│   ├── test_fibrosis.py
│   ├── test_h5_backed.py
│   ├── test_io.py
│   ├── test_mcp.py
│   ├── test_phase_diagram.py
│   ├── test_public_evidence.py
│   ├── test_release.py
│   ├── test_sota.py
│   ├── test_spectral.py
│   ├── test_uc_xenium.py
│   ├── test_visium.py
│   ├── test_visium_archive.py
│   └── test_xenium_archive.py
├── .dockerignore
├── .gitattributes
├── .gitignore
├── Dockerfile
├── LICENSE
├── pyproject.toml
└── README.md
```
