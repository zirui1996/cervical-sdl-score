# Proliferation-adjusted expression differences between cervical precancer and invasive carcinoma

Fixed on 2026-10-29 after outcome-blind metadata review, before any expression-stage statistical testing in these cohorts. This is a dated local analysis plan, not a public preregistration. Purpose: feasibility of a pure computational, independently validated study. No acceptance or novelty is presumed.

## Primary question and cohorts

Which epithelial-enriched expression differences between CIN3 and invasive squamous carcinoma persist after adjustment for a prespecified proliferation score and reproduce across independent microarray cohorts?

Discovery GSE63514: primary set includes 40 CIN3 and 20 laser-captured cancer samples; exclude the 8 whole-section cancers from the primary test. All 128 samples are retained for descriptive stage summaries. Sensitivities: all CIN3/cancers with a dissection-method covariate; CIN2+3 versus laser-captured cancer. Source metadata labels cancers as cervical squamous epithelial cancer. A cross-sectional lesion comparison is not longitudinal progression, conization upgrading, or prospective diagnosis.

Primary external GSE138080: 15 CIN2/3 and 10 SCC, laser-capture enriched. This is a broader precursor-grade generalization, not an exactly matched CIN3 endpoint. Normal samples are descriptive only. Secondary external GSE7803: 7 HSIL and 21 SCC. Exclude 3 explicitly marked Test Sample arrays, which are technical/test samples, before expression testing. Keep the 10 normals for descriptive analysis. The cohorts derive from distinct source studies (SUCCEED; Amsterdam/VUmc; Michigan study). Preserve sample identifiers and source notes; do not claim a completed patient-identity overlap audit beyond available public metadata.

## Processing and proliferation score

Use author-processed matrices, preserve downloaded bytes and metadata, check values against source processing notes and transform positive linear intensity to log2 only if its scale requires it (99th percentile >100 or range >50). Do not normalize or batch-correct across cohorts. Exclude ambiguous gene/probe mappings and summarize valid probes by gene-level median. Keep a complete mapping audit.

Freeze the discovery gene universe to symbols measurable by all three platform annotations (GPL570/GPL4133/GPL96), without examining external outcome-linked expression. Require complete discovery values, positive variance, pooled discovery median log2 intensity >=4 and IQR>=0.25. Exclude all prespecified proliferation-score genes from the tested target universe to prevent mathematical self-correlation.

Proliferation marker definition: MSigDB Hallmark E2F_TARGETS union G2M_CHECKPOINT, release 2025.1.Hs. Score: mean within-sample fractional expression rank of available marker genes over the annotation-based common platform gene universe; standardize the score within each analyzed cohort. Record gene coverage and score-group overlap. Sensitivity scores: E2F only, G2M only, and a fixed marker panel (MKI67, PCNA, MCM2, MCM3, MCM4, MCM5, MCM6, MCM7, TOP2A, CDK1, CCNB1, CCNB2). These scores are transcriptional proxies, not measured proliferation.

## Primary testing and candidate lock

Per gene fit expression ~ intercept + invasive_indicator + proliferation_score, using HC3 robust covariance and a two-sided t reference with residual degrees of freedom. Report invasion coefficient, 95% CI, raw P and BH FDR over the entire tested gene universe. Also report unadjusted group difference and FDR; changes in significance alone do not establish attenuation or mediation.

Qualifying candidates: adjusted FDR<0.05 and absolute adjusted log2 expression coefficient>=0.5. Freeze at most 5 candidates ranked by adjusted P (gene symbol breaks ties). Write gene identifiers, discovery direction, full statistics, gene universe, plan hash, code hash and UTC timestamp to a candidate lock BEFORE loading external expression values. If no gene qualifies, do not promote nonsignificant genes as validated discoveries; report the null outcome.

## External confirmation and robustness

Apply the same primary expression-on-group-plus-proliferation model separately in each external cohort. Do not select genes, direction, score definitions or cutpoints from external results. Correct candidate P values by Holm within each external cohort, including unmeasured/non-estimable locked candidates as P=1. Replication in a cohort requires corrected P<0.05 and direction concordance. Stronger support requires the same candidate to meet this rule in both external cohorts. Do not pool studies to rescue failed replication.

Report the design condition number, group/proliferation collinearity and score ranges. Common-support sensitivity restricts samples to the intersection of both groups' score ranges, only when >=5 samples remain in each group. Other sensitivities: alternate proliferation scores; a quadratic score term; leave-one-out candidate coefficients; whole-section adjustment in discovery; and CIN2+3 discovery. All are supportive and cannot replace a failed primary external test. Per-cohort unadjusted candidate AUC is descriptive only, with direction fixed in discovery and stratified bootstrap CI; no optimized thresholds or clinical diagnostic claims.

## Cellular localization and novelty

Reuse GSE171894's existing analyzed object for descriptive donor-by-lineage summaries of all locked candidates. Preserve its provisional lineage annotations and author-normalized expression scale. Do not call epithelial cells confirmed malignant cells, or treat individual cells as patients. No metastasis or temporal-transition inference from this localization cohort. If a target is too sparsely measured, report that limitation.

Search all locked candidates and the closest designs after selection without replacing genes to chase novelty or significance. A lack of retrieved identical papers is not proof that no prior study exists. Distinguish a known gene/new tested question from a genuinely unstudied molecule. Public-data originality and journal suitability must be reassessed from actual results.

## Deliverables and stopping rule

Keep sources with timestamps and SHA256; exact phenotype records; processing logs; every tested gene's results; candidate lock; all external, failed and sensitivity results; code and package versions; figures; and an honest feasibility report. Do not draft a positive biomarker manuscript unless both methodological quality and validation support it. No causal, treatment-response, or future-progression claims are licensed by this design.
