# Single-study versus joint-discovery ranking benchmark

Fixed 2026-10-29 after the thirteen frozen target/program outcomes were known in new cohorts, before testing all other new-cohort genes. This is a post hoc benchmark; it is not presented as the original primary endpoint. No result changes the ten-gene set or prior validation decisions.

Question: does joint discovery yield more reproducible gene-wise association rankings in new studies than a single discovery study? Compare the frozen original GSE63514 signed HC3 normal-score ranking with the already fixed three-study weighted-Stouffer ranking. Use all measurable genes in the 6,388-gene joint-discovery universe, not genes selected by a new-cohort P value. External outcomes are signed normal scores from expression ~ pathology + Composite, same HC3/t inference, in GSE166466, GSE149763 and GSE147009. Do not fit classifiers or choose a new threshold.

Report Spearman correlation of each frozen ranking with each new-cohort ranking and their paired difference (joint minus single). Also report all three individual older-study correlations as context. Normalize signs identically. Rank-based measures tolerate platform effect scales but do not establish clinical prediction or biology.

For uncertainty, perform 1,000 fixed stratified participant-bootstrap attempts within each new cohort with seed 20261031, holding the development rankings and test-gene universe fixed. Recompute the adjusted HC3 statistics and rank correlations in each resample. Record rank-deficient or non-finite fits rather than silently adding more attempts; report valid counts and pointwise percentile intervals. These are exploratory intervals, conditional on frozen development data/universe and unadjusted for the three cohort comparisons. If improvement is unclear or inconsistent, report that explicitly. Small RNA cohorts may have broad intervals and unstable resamples.

The benchmark adds a directly testable justification for the joint-discovery strategy; it must not be used to retroactively relabel the new program as prospectively designed or to hide any failed cohort.
