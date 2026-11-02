# Epithelial subtype and differentiation context

Written 2026-10-29 after the primary validation/GO results and review of prior cervical FCGBP literature. This is a new exploratory sensitivity prompted by a concrete biological alternative: FCGBP has been localized to endocervical columnar/glandular epithelium rather than ectocervical squamous epithelium (see https://pmc.ncbi.nlm.nih.gov/articles/PMC10888570/). Neither FCGBP reduction nor a bulk score automatically identifies an invasive squamous-cell mechanism.

Freeze independent marker sets, excluding all thirteen tested genes and proliferation-score genes:
- Glandular/secretory context: PAX8, TFF3, AGR2, MUC5B, MUC5AC, KRT7.
- Squamous basal identity: KRT5, KRT14, TP63, DSG3, DSC3, SFN.
- Terminal squamous differentiation: KRT4, KRT13, IVL, SPRR1A, SPRR3, KRT1, DSG1.

In each of the seven tissue cohorts, calculate marker mean fractional-rank proxies over the same fixed gene-ranking universe, record coverage, and fit separate additional-covariate models of the frozen program and thirteen genes: original group + Composite + one of these proxies. Apply BH across all extension target/model/cohort tests. These proxies may represent biology downstream of lesion state and are not proven confounders; loss of significance is not proof of admixture, nor persistence proof of cell-intrinsic invasion. Record variance inflation and avoid mechanistic interpretation when collinearity or sample size makes a model unstable.

For GSE208653, use the already retained epithelial cells and fixed independent markers to describe marker-high epithelial states (glandular, squamous basal, differentiated squamous, mixed/unspecified). Require >=2 detected markers, mean normalized marker score >=0.3 and margin >=0.15, otherwise retain mixed/unspecified. Keep specimen as the unit, require >=30 cells for program summaries, and report all states/donors including ADC separately. States are expression labels, not validated histology or malignant identity. Do not use program genes, known pathology labels or desired direction to assign states. This analysis cannot upgrade four squamous/HSIL specimens to powered confirmation.
