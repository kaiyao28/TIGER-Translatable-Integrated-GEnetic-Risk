# Getting started with TIGER

TIGER is deliberately transparent: installation changes how functions are
loaded, not how the method is reviewed. The complete implementations remain as
plain R source in the repository's `R/` directory.

Before applying TIGER, prepare a correctly centred liability-scale PRS and a
liability-scale PRS R-squared. The main interface accepts one row per person,
with optional group, RV-status, and high-impact genotype columns. High-impact
genotypes may represent a single variant coded 0/1/2 or a multi-variant system
such as the six APOE genotypes.

```r
library(TIGER)

K <- 0.01
SP <- 0.50
r2_liability <- 0.10
individuals <- data.frame(
  ID = c("P001", "P002"),
  PRS_liability = c(-0.2, 0.4)
)

results <- tiger_probabilities(
  individuals, K = K, SP = SP,
  method = "PAIR (summary)", r2_liability = r2_liability
)
```

Before calculation, `check_tiger_inputs()` can validate IDs, liability-scale
PRS values, RV matching, and high-impact genotype matching. Hypothetical CSV
schemas are available through `tiger_template_file()`.

When modelling a high-impact component separately, exclude that variant or LD
region from the PRS and its R-squared estimate. The RV probability prior is
`0.50` for a balanced case-control analysis. Use a population-level prior only
when the reference evidence and intended interpretation are population-level.

Use `?tiger_probabilities`, `?TIGER-high-impact`, and `?TIGER-plots` for the
main interfaces. Full input, liability-scale PRS, method, and plotting guides
are available from the
[TIGER repository](https://github.com/kaiyao28/TIGER-Translatable-Integrated-Genetic-Risk).
