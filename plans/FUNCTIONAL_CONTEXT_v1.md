# Exploratory functional annotation

Written 2026-10-29 before downloading/testing GO sets. Use MSigDB c5.go.bp.v2025.1.Hs.symbols.gmt, the same release as the frozen Hallmark source. This is functional context, not a mechanistic experiment or classifier.

Background is the fixed 6,388 joint-discovery genes. Foregrounds are (1) all ten frozen program genes, (2) all jointly qualified positive-direction genes, and (3) all jointly qualified negative-direction genes using the unchanged meta-discovery thresholds. Test all GO biological-process sets containing 10 to 500 background genes with the hypergeometric upper tail; include P=1 tests, and BH-correct across all foreground-by-term tests jointly. Report complete results, foreground sizes, overlaps and corrected P. Do not relabel the program as epithelial differentiation or a specific pathway unless its own annotation evidence supports that wording. Annotation enrichment does not establish regulation, causal mechanism or cell specificity.
