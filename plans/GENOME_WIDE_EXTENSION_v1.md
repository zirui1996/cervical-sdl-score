# Genome-wide adjustment and reproducibility extension

Written 2026-10-29 before running the extension, after all prior candidate analyses and validation outcomes were known. This is a new exploratory extension, not a retrospectively prespecified primary endpoint. It does not replace either candidate lock or redefine validation success.

## Question

How does conditioning on a transcriptomic proliferation score change cervical high-grade precursor versus invasive-carcinoma expression associations, and does it improve or reduce cross-study agreement? The aim is a genome-wide reproducibility assessment, not a new diagnostic classifier or a claim that adjustment necessarily helps.

## Fixed design

Use the existing primary samples in GSE63514, GSE138080 and GSE7803, the existing processed log-expression matrices, frozen Composite score and original discovery universe of 6,685 genes. Exclude the predefined score genes as previously specified. Do not select new markers from this extension.

1. Refit unadjusted expression ~ group and adjusted expression ~ group + Composite score in each cohort. In cohort-level tables use the complete measurable subset of the fixed discovery universe, excluding genes with missing or zero-variance primary expression. HC3 inference and within-cohort/model BH correction cover that subset; these new external FDR values are separate from the earlier 257-gene extension's FDR.
2. In discovery report full-universe adjusted versus unadjusted FDR counts, overlaps, direction reversals and coefficient changes. Do not call loss of significance proof that proliferation explains a gene, or gain of significance proof of a new mechanism.
3. Define a common measurable universe as genes from the frozen discovery universe with complete and variable expression in all three cohorts. Compare standardized group coefficients (native coefficients divided by within-study gene SD) across each of the three study pairs. Report Spearman correlations and sign concordance for unadjusted and adjusted models over all common genes without selecting by external P values. Near-zero effects can have unstable signs, so sign concordance is descriptive.
4. Quantify uncertainty in the paired difference in adjusted versus unadjusted Spearman correlations using 1,000 stratified participant bootstraps independently within each study and pathology group, fixed seed 20261030. Refit both models and recompute expression SD in every resample. Use 2.5/97.5 percentiles. These are exploratory, unadjusted pointwise intervals conditional on the fixed gene universe, not simultaneous confirmatory tests. Resample participants, not genes; no arbitrary significance threshold determines what is reported.
5. Report continuous coefficient-change distributions and discovery-only FDR retention. Keep the three earlier exploratory genes in their existing evidence framework; no new three-gene score is introduced.
6. Assess sampling sensitivity in discovery across the full 6,685 genes: primary 40 CIN3 versus 20 laser-captured cancers; include all 28 cancers without dissection adjustment; include all 28 with the whole-section indicator. Also compare the 20 laser-captured versus 8 whole-section cancers directly with and without Composite adjustment. Label this an observational sampling association; dissection is not randomized and is confounded with specimen characteristics. The within-cancer contrast is not proof of technical contamination or a causal dissection effect.

## Reporting

Report all study pairs and both models whether agreement improves, worsens or remains uncertain. Do not claim novel statistical methodology, causal separation of proliferation from invasion, temporal malignant progression, or clinical prediction. Preserve every original failed candidate and the unconfirmed new RNA cohort. Package scripts, full coefficient tables, bootstrap draws, sources and figures with the manuscript.
