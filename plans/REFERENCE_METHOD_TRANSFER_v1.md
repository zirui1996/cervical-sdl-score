# Method transfer and new single-cell source appraisal

Date: 2026-10-29. The reference is the authors' endometrial carcinoma
null-controlled manuscript. Its findings are read as manuscript-reported results,
not independently rerun or assumed published. Existing cervical outcomes are known.
This is a dated post hoc extension; no old candidate or primary result is changed.

## Transfer of principles

Use explicit cohort roles, selection-aware controls, full reselection when a
cohort is held out, patient/specimen-level inference, and compartment localization.
Do not turn the cervical project into another high-AUC/random-panel critique.
Existing reference panels remain supportive, with their documented matching flaw.
Do not gate validation cohorts on having significant differential expression:
HSIL-versus-invasion is a closer contrast than tumor-versus-normal and exclusion
after weak results would select validation on the desired outcome.

## Leave-one-discovery-study-out reselection

For each of the three old microarrays, recompute the original signed sqrt(n)
Stouffer ranking using only the other two studies, over the unchanged 6,388-gene
universe. BH the full universe; require same direction in both training studies,
FDR < .05 and mean absolute standardized coefficient >= .5; select first ten by
combined P and gene name, without a fallback if fewer than ten qualify. Freeze
fold-specific genes/signs and test the left-out study with exactly the original
equal within-sample rank program and score ~ group + Composite HC3 model. Retain
every fold, multiplicity-adjust over three fold-level tests, and report overlap
with the existing program. This evaluates the discovery rule after the fact,
not the original ten-gene set as an untouched validation or a clinical classifier.
No new LASSO, SVM-RFE, weights, threshold or AUC optimization is introduced.

## Additional data source audit

OMIX001005/PRJCA008573 contains 13 reported cervical single-cell specimens with
normal, HSIL and early/advanced SCC. Public listing is not proof of successful
download: sampled official count-file links returned 404/timeout. Save the audit.

The official EBI API for E-MTAB-12305 reports ten assays: H1/H2 HSIL, T1-T4
cervical tumors, N1-N3 normal-adjacent tissue and L1 metastatic lymph node. Verify
every mapping against SDRF. Patient mapping and histology must be retained;
adjacent tissue and lymph node cannot become HSIL comparators. This study's
description emphasizes immune and mesenchymal cells: first inspect whether its
released count matrices include adequate epithelial cells before claiming it
can test epithelial localization. Multiple assays may come from the same patient.

If available, analyze native sparse count matrices plus features/barcodes for
H1/H2, T1-T4 and N1-N3, sequentially in bounded memory. L1 is excluded from
cervical-tissue comparisons. Use the same QC as the existing GSE208653 extension
(200-7500 detected genes, <=20% mitochondrial fraction, positive counts), the same
canonical broad marker lists and confidence cutoffs, and retain ambiguous cells.
Check nonnegative integer units before summing. No integrated or cell-cycle-
regressed values enter expression inference. No target genes define compartments.

Summarize count-sum pseudobulk by specimen/compartment. Apply the frozen ten-gene
score, original ranking universe and >=8 target coverage. Require >=30 cells for
a compartment score. Show each specimen, source-group direction and target
localization. Primary descriptive contrast: H1/H2 versus T1-T4 in the epithelial
compartment; All cells is explicitly a capture-dependent pseudo-bulk, not actual
bulk RNA. Normal-adjacent specimens supply context only. If the epithelial
compartment is absent or poorly sampled, declare the source unsuitable for that
question instead of treating immune expression as epithelial validation.

The two HSIL specimens do not support a large adjusted donor model. No cell-wise
P values, temporal trajectories, inferred causal interactions, CNV-defined
malignancy or treatment-response claims. Cross-source cell totals are never
patient totals; the old discordant source remains displayed even if a new source
agrees. This plan concerns fixed-feature localization, not a new feature search.
