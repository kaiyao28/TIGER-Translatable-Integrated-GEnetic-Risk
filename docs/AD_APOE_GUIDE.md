# Applying TIGER to Alzheimer’s disease and APOE

APOE is a common, high-impact genetic component and should be handled separately
from both the polygenic score and the rare-variant layer.

## Before calculating the PRS

If APOE will be added separately, normally exclude the APOE region from:

- discovery GWAS variants used to build the PRS;
- LD reference data used by the PRS method;
- target and population-reference PRS calculations; and
- estimation of PRS R².

Define and report the excluded genomic interval, genome build and variant list.
The exact interval is an analysis decision: use a defensible LD-aware definition
for the study population rather than assuming one universal boundary. Apply the
same exclusion to every dataset. If a published PRS already contains the APOE
region and cannot be rebuilt, do not add APOE as though it were independent;
that would double count some of the same information.

## Inputs

The main TIGER interface requires one APOE genotype value per person in the
individual table. If APOE starts in a separate `ID, APOE_genotype` file, join it to the
prepared PRS input by ID and stop on missing or duplicated records. TIGER also
requires either:

1. observed frequencies of the six APOE genotypes among comparable cases and
   controls (preferred); or
2. e2/e3/e4 allele frequencies among cases and controls, with Hardy-Weinberg
   equilibrium used to derive genotype frequencies.

The evidence should match the target disorder definition, ancestry, age context
and population as closely as possible. Frequencies are proportions from 0 to 1,
not odds ratios or percentages.

## Minimal application

```r
source("R/probability_methods.R")
source("R/high_impact_variants.R")
source("R/rare_variant_probability.R")

# Supply justified values from the intended AD population and evidence source.
K <- 0.01
SP <- 0.50
r2l <- 0.10
SP_RV <- 0.50
prs_liability <- c(-0.20, 0.10, 0.35)
apoe <- c("e3/e3", "e3/e4", "e4/e4")

p_prs <- pair_probability_summary(prs_liability, K, SP, r2l)

apoe_reference <- apoe_genotype_reference(
  case_allele_frequencies = c(e2 = 0.04, e3 = 0.66, e4 = 0.30),
  control_allele_frequencies = c(e2 = 0.08, e3 = 0.79, e4 = 0.13)
)
p_prs_apoe <- high_impact_method_probability(
  prs_liability, apoe, apoe_reference,
  K = K, SP = SP, method = "PAIR (summary)", r2_liability = r2l
)

# If independently supported RVs are present, add them last.
p_rv <- intrinsic_rv_probability(odds_ratio = 4, prevalence = SP_RV)
p_prs_apoe_rv <- apply_rv_probability(p_prs_apoe, p_damaging = p_rv)
```

The allele frequencies above are illustrative only and must not be treated as
recommended AD parameters.

`SP_RV = 0.50` gives TIGER's balanced-sample RV calculation. Use
a population-level background probability only when the frequency/effect
evidence and intended probability are explicitly population-level.

## Interpretation and checks

The method-specific APOE update derives genotype-specific population prevalence
and sample prevalence values and recalculates the selected PRS conversion. It assumes that APOE information
has been removed from the PRS. Adding an RV afterward additionally assumes that
the RV effect is independent of the PRS and APOE component.

Before use, check:

- whether `K` is appropriate for the population, age horizon and outcome;
- whether `SP` represents the intended target context;
- whether APOE genotype or allele frequencies match the target population;
- whether HWE is reasonable if allele frequencies are used;
- whether the PRS and APOE/RV evidence overlap; and
- calibration in an independent dataset representative of the intended use.

TIGER estimates disorder probability under stated assumptions. It does not
provide diagnosis, prognosis, treatment recommendations or clinical validity.
