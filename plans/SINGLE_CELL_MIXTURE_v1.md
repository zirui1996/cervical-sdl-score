# Captured-RNA mixture sensitivity

Written 2026-10-29 before this decomposition was computed; existing bulk and single-cell results are known. This is post hoc.


Another useful extension goes beyond the reference manuscript's three-marker
composition covariate. In GSE208653, restrict to the broad compartments having
at least 30 cells in every HSIL/SCC specimen. Keep ambiguous cells as their own
category if eligible. For each specimen, normalize each compartment's count-sum
pseudobulk to CPM; its weight is the compartment's share of total captured UMI
counts over eligible compartments, NOT a true tissue cell fraction. Check that
the weighted mixture exactly reconstructs the restricted all-cell pseudobulk.
Display the retained cell and UMI coverage relative to the full specimen.

For all four cross-group specimen pairs (2 HSIL x 2 SCC), decompose the fixed
score difference using the two-order Shapley average. For score f(w,E), composition
is [(f(wS,EH)-f(wH,EH))+(f(wS,ES)-f(wH,ES))]/2; profile is
[(f(wH,ES)-f(wH,EH))+(f(wS,ES)-f(wS,EH))]/2. They must sum to the observed
restricted-mixture contrast despite nonlinear within-sample ranking. Also show
every specimen at the common mean weight vector. Do not interpret a small net
contrast with opposing components as a meaningful percentage explained.

These are descriptive synthetic mixtures of captured RNA. They do not estimate
the fraction of real bulk change caused by tissue composition, true in-vivo
fractions, causal mediation or malignant epithelial regulation. Broad compartments
can hide within-compartment state changes. Do not bootstrap cells or count the
four crossing pairs as four independent patient pairs.
