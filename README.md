# Supplementary code and records

Title: Proliferation-adjusted loss of squamous differentiation-associated gene expression in invasive cervical squamous cell carcinoma compared with high-grade precancer: a multi-cohort transcriptomic analysis

Contents
- scripts/: analysis scripts (Python 3.12; dependencies in requirements.txt).
- plans/: dated analysis plans and addenda, byte-identical to the versions whose hashes appear in the lock files, except for the three files listed under "Editorial changes to plan files". They record the analysis chronology, including failed earlier candidates.
- records/: frozen-score lock files, source-file download manifests with checksums, and the record of the post hoc checks.
- tables/: Supplementary Tables S1-S22 as tab-separated files (SHA256 hashes are in the Index sheet of the workbook).
- full_result_outputs/: complete outputs with the working file names used by the scripts, including all-target rows of each testing family and technical replicates; table_number_map.tsv maps working numbers to the final tables.

Raw data are public and are not duplicated. Accessions are listed in Table 1 and Supplementary Table S01 of the article; the manifests identify the exact source files. Scripts expect the project layout data/raw, data/processed, and results relative to the parent of scripts/. Reproduction keeps the saved locks and does not create a new holdout experiment.

The saved post hoc checks (scripts/revision_posthoc_checks.py) produced the outputs of Tables S4-S8 from processed matrices. The processed GSE14404 gene-level matrices must be regenerated from the raw files with scripts/analyze_GSE14404_independent_tissue.py; those checks use the saved GSE14404 score file.

Editorial changes to plan files
One working note was reworded in each of the three files below before this repository was made public. Analysis rules, cohorts, candidates, and results are unchanged. The original files are retained by the authors and are available on request. The recorded hashes refer to the originals, so the plan-hash checks in scripts/verify_results.py and scripts/verify_reference_transfer.py will not pass on these copies.
- plans/ANALYSIS_PLAN_v1.md: original SHA256 d250661502c12d9b3bcd569cf2a7c27b0e4be0bff76d5eab7c75b7d7bd4a144c (records/candidate_lock.json)
- plans/ADDITIONAL_COHORT_RULES_v2.md: original SHA256 ff05af6f01b42ba0c8306aac8ad2eeb43c4603dc32b0a9f833ae93f0aa92ab37 (records/additional_cohort_candidate_lock.json)
- plans/REFERENCE_METHOD_TRANSFER_v1.md: original SHA256 f2c63a0972dea4ae1ced9fac6fd309f24f91e9ff3474002a2fab9dd71217492a

Command order
```bash
.venv/bin/python scripts/download_data.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_analysis.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_exploratory_extension.py
.venv/bin/python scripts/lock_new_cohort.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/validate_RNAseq.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/check_RNAseq_scale.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/summarize_singlecell.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/review_publication_robustness.py
.venv/bin/python scripts/verify_results.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_genomewide_extension.py
.venv/bin/python scripts/download_extension_data.py
.venv/bin/python scripts/download_extension_data.py --singlecell
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/lock_meta_discovery.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_additional_bulk.py
.venv/bin/python scripts/request_archs4_RNA.py GSE147009
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/analyze_archs4_RNA.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/summarize_original_gene_followups.py
.venv/bin/python scripts/request_archs4_RNA.py GSE223804
.venv/bin/python scripts/request_archs4_RNA.py GSE149763
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/compare_RNA_processing.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/analyze_singlecell_extension.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/analyze_stage_and_composition.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/clarify_prospective_endpoint.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/program_robustness_checks.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/benchmark_discovery_strategy.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/functional_context.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/epithelial_context_checks.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/extension_summary_figures.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_extension_cohorts.py
.venv/bin/python scripts/verify_results.py
.venv/bin/python scripts/verify_extension.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/appraise_program_distinctiveness.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_program_appraisal.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/nested_discovery_rule_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/decompose_captured_RNA_mixture.py
# Extract cervical_EMTAB12305_sources_part1/2/3.zip first; the provided source
# bytes include complete metadata and 27 sparse matrix/feature/barcode files.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/analyze_EMTAB12305.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/EMTAB12305_epithelial_states.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/render_singlecell_source_comparison.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_reference_transfer.py
OPENBLAS_NUM_THREADS=1 python scripts/revision_posthoc_checks.py
```
