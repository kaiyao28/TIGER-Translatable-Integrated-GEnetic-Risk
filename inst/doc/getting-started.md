# Getting started with TIGER

TIGER is deliberately transparent: installation changes how functions are
loaded, not how the method is reviewed. The complete implementations remain as
plain R source in the repository's `R/` directory.

Before applying TIGER, prepare a correctly centred liability-scale PRS and
liability-scale PRS R-squared. The main interface accepts one row per person,
with optional APOE and RV-status columns. If APOE is modelled separately,
exclude the APOE region from PRS construction and prepare the matching
reference PRS in the same way.

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

The default RV background is `prevalence = 0.50`, representing the balanced
50:50 case/control sample calculation. Supply a justified population-level
probability only when the reference evidence and intended interpretation are
population-level.

The repository provides the full guides under `docs/`: `GETTING_STARTED.md`,
`USAGE.md`, `LIABILITY_PRS_GUIDE.md`, `AD_APOE_GUIDE.md`, and `METHODS.md`.
