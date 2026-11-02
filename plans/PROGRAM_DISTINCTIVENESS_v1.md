# Fixed-program distinctiveness and uncertainty extension

Written 2026-10-29 before running the analyses below. All original target outcomes,
including discordant RNA and single-cell findings, are already known. This is an
explicitly post hoc scientific appraisal, not a prospectively registered validation.
The existing ten genes, signs, scoring rule and primary results remain fixed.

## Question

Does the fixed ten-gene set carry more invasion-associated evidence than other
equally sized gene sets chosen without using any new-study outcomes? Distinguish
an external association from evidence that the exact gene membership is special.

## Fixed reference construction

- Seed 2026102907; 2,000 draws per reference family; ten distinct genes per draw.
- Exclude the union of the original three targets and ten program genes from all
  reference pools. The discovery universe is the existing 6,388 common measured
  genes. No screening on new expression, effects, P values or coexpression.
- Family A: uniform random ten-gene sets from this discovery universe.
- Family B: uniform sets from the existing jointly qualified, same-direction
  negative genes (combined BH < .05; mean absolute standardized beta >= .5).
- Family C: match each target to one of its 20 nearest available negative
  qualifying genes, using Euclidean distance in standardized discovery-only
  features: magnitude of combined Z, mean absolute standardized beta, mean
  within-sample expression rank and mean within-cohort rank SD. Scales are SDs
  among the negative qualifying pool. Choose ten unique genes using random
  target order and uniform choice among unused nearest neighbours. Save candidate
  distances and realized matching diagnostics; imperfect matching is a limitation.
- Score every reference set with the same negative equal-rank rule as the frozen
  all-negative ten-gene program. No weights or thresholds are trained. The ranking
  universe stays exactly the original common platform universe intersecting the
  measured genes of each study. Require at least eight measured genes per panel;
  retain unscorable draws and report their number rather than redrawing them.
- Also display a nonoverlapping broad negative-discovery score, using every
  qualifying negative gene except the 13 frozen targets. This descriptive baseline
  is not a new selected signature and does not replace the primary program.

## Outcomes, fixed before computation

Use the three newly analyzed tissue studies GSE166466, GSE149763 and GSE147009
(the latter two are small RNA studies). Fit score ~ group + Composite, using HC3
and residual-df t inference as before. Primary reference statistic: signed t for
the fixed invasion orientation. Also report standardized beta and the fraction
of panels with same-direction nominal P < .05. Do not select the best panel.

In GSE166466 also evaluate the direct late-minus-early contrast from
Normal/HSIL/Cancer + Composite, vector (0,-2,1,0), using signed HC3 t. This makes
four contexts times three reference families = 12 descriptive comparisons.
Report finite-panel rank and the add-one exceedance fraction
(1 + # reference t >= fixed t)/(1 + # scorable panels), plus Holm over these 12
fractions. These are conditional reference-panel comparisons, NOT phenotype
permutation P values, false-positive guarantees, clinical validation, or proof of
biological specificity. Gene-gene correlation is not matched. Show all families.

## RNA uncertainty

Retain individual 95% HC3 intervals. Show their relation to zero and to the point
estimate from GSE166466; overlap is descriptive, not a test of equal effects.
Optionally summarize two strictly new RNA cohorts using an exploratory
Paule-Mandel random-effects estimate and modified Knapp-Hartung interval with
variance factor max(1,q), t df=k-1. Also show all three RNA studies separately
labelled as including the previously examined GSE223804. Never combine discovery
microarrays with validation to manufacture confirmation. Do not infer platform
causation from between-study patterns, compute observed power, or claim practical
equivalence without a justified margin. Preserve covariance, small-k and score-SD
comparability caveats. Pooling is secondary estimation, not a new replication.

## Decision

Use findings to narrow the research claim even if the program does not outperform
discovery-matched sets. No gene replacement, alternative score selection, new
patient counts, or mechanism/clinical efficacy claims follow from this exercise.
