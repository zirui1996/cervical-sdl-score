# Joint discovery before opening additional tissue expression

Written 2026-10-29 after source/phenotype review and before opening GSE166466 or GSE149763 expression matrices. The matrices may already be downloading; no new expression values or gene-wise outcomes have been inspected. This is a post hoc redesign of discovery, disclosed alongside the failed original single-discovery screen and subsequent three-gene evidence. The older external cohorts are discovery inputs here, not untouched validation.

## Rationale

Ranking genes only in one discovery study can prioritize cohort-specific effects. Use all three previously analyzed microarray studies as a joint discovery stage, then test a frozen small set in the newly acquired tissue cohorts. This adds an independently testable design rather than replacing genes after their new-cohort results are seen.

## Fixed joint-discovery rule

Use all 6,388 complete, variable genes from the original frozen 6,685-gene discovery universe that are measured in all three older microarray cohorts, as established by the genome-wide extension. No new-cohort data determine the universe or ranking. Keep cohort-specific expression ~ pathology + Composite models and HC3 two-sided P values.

For each gene, transform each study's two-sided P into a signed normal score using the coefficient sign. Combine the three scores with fixed sqrt(n) weights (n=60,25,28), normalized by sqrt(sum(n)); use a two-sided standard-normal tail for the combined P. Correct all 6,388 combined P values by BH. This is a weighted signed Stouffer evidence summary, not a pooled common effect or a heterogeneity-free causal estimate.

Qualify genes only if all three coefficients have the same nonzero sign, combined FDR<0.05 and the mean absolute within-study standardized coefficient is at least 0.5. Order by combined P then symbol, freeze at most 10. If fewer qualify, keep that number; do not relax thresholds. Record every tested gene, selected directions, sample weights, plan hash, source-table hash and timestamp. Do not claim new statistics or novelty of these genes.

## New-cohort testing

Test the union of these at most ten genes and the three earlier frozen genes. Preserve the earlier three-gene family as a historical primary follow-up, and report (a) Holm over the meta-ranked set within each new study, (b) Holm over the full union within each new study, and (c) Holm over the full union across both new studies. Unmeasured/non-estimable targets get P=1. Candidate replication requires same direction and corrected P<0.05 under the explicitly identified family; no threshold or set is changed based on new results. The combined two-study family is the conservative headline for any new positive candidate claim. Author RPKM and unadjusted fits remain sensitivities, never substitute for a failed primary adjusted model.

Other cohorts have fixed context roles: stage specificity (GSE51993), provisional distinct prospective endpoint (GSE75132), specimen-level single-cell localization (GSE208653). They are not replacement invasive-cancer confirmation datasets. Do not train a new diagnostic panel on these small validation studies.
