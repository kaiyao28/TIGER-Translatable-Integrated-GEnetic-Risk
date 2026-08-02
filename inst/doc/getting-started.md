# Getting started with TIGER

TIGER is deliberately transparent: installation changes how functions are
loaded, not how the method is reviewed. The complete implementations remain as
plain R source in the repository's `R/` directory.

Before applying TIGER, prepare a correctly centred liability-scale PRS and
liability-scale PRS R-squared. Keep individual PRS, APOE or other common
high-impact status, and RV carrier records separate. If APOE is updated
separately, exclude the APOE region from PRS construction and prepare the
matching reference PRS in the same way.

```r
library(TIGER)

K <- 0.01
SP <- 0.50
r2_liability <- 0.10
prs_liability <- seq(-1, 1, by = 0.25)

p_prs <- pair_probability_summary(
  prs_liability,
  K = K,
  SP = SP,
  r2_liability = r2_liability
)
```

The default RV background is `prevalence = 0.50`, representing the balanced
50:50 case/control sample calculation. Supply a justified population-level
probability only when the reference evidence and intended interpretation are
population-level.

The repository provides the full guides under `docs/`: `GETTING_STARTED.md`,
`USAGE.md`, `LIABILITY_PRS_GUIDE.md`, `AD_APOE_GUIDE.md`, and `METHODS.md`.
