# Fixed epithelial marker-state context in the additional source

Written 2026-10-29 after the E-MTAB-12305 broad-compartment program results were
known: all four tumor epithelial scores exceed both HSIL scores, whereas the
older GSE208653 source disagrees. This follow-up is post hoc and descriptive.

Apply exactly the pre-existing three epithelial marker lists and decision rules
in epithelial_context_checks.py, without selecting markers using either source's
program association. Restrict to the six primary E-MTAB-12305 HSIL/SCC specimens
and the already assigned epithelial cells. Require marker mean >= .3, >=2 detected
markers, and top-minus-second >= .15; otherwise Mixed_unspecified. Summarize state
fractions and count-sum pseudobulk fixed-program scores per specimen/state; below
30 cells do not report a score. Compare with existing GSE208653 descriptions.

These are marker-high expression states, not pathology-confirmed glandular tissue,
malignant identities or measured tissue admixture. Some broad markers, such as
KRT7 and AGR2, can be expressed in squamous tumors. No cell-wise tests, calibrated
cell proportions, causal explanation or feature/threshold changes. Retain the
older discordant source. If no state is adequately sampled in all specimens,
report that limitation rather than pooling cells to claim validation.
