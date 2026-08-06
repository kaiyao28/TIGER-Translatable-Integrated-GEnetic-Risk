# TIGER method guide

TIGER separates the initial PRS probability conversion from two distinct,
optional components: common high-impact genotypes such as APOE, and rare
variants. Users may calculate PRS only, PRS + RV, PRS + APOE/high-impact, or a
combined PRS + APOE/high-impact + RV probability. APOE is not classified as an
RV in TIGER.

## Input prerequisite

All PRS probability methods in TIGER require scores prepared on the scale
expected by the selected method. The TIGER examples use a centred
**liability-scale PRS** and matching liability-scale R². Prepare the target and
population-reference scores consistently before probability conversion; do not
use an unconverted raw or merely standardized PRS. See
[`LIABILITY_PRS_GUIDE.md`](LIABILITY_PRS_GUIDE.md) and the runnable TIGER example at
[`examples/liability_conversion_example.R`](../examples/liability_conversion_example.R).

## Default worked-example settings

```r
K <- 0.01
SP <- 0.50
r2_liability <- 0.10
SP_RV <- 0.50
```

These defaults describe an illustrative 1%-prevalence, balanced-target
scenario. They are not universal disease parameters. Replace or justify them
for the intended population, phenotype and application.

The default worked-example method is **PAIR (summary)**:

```r
p_prs <- pair_probability_summary(
  prs_liability,
  K = K,
  SP = SP,
  r2_liability = r2_liability
)
```

## PRS probability methods

### PAIR (summary)—default

`pair_probability_summary()` derives theoretical case/control PRS moments from
population prevalence `K` and liability-scale PRS R². It is phenotype
independent and is the default in the TIGER examples and plot labels.

### PAIR (sample)

`pair_probability_sample()` uses case/control PRS means and standard deviations
estimated from a labelled calibration sample. Moments estimated in the people
being scored give an in-sample descriptive result, not prospective validation.

### BPC

`bpc_probability()` applies Bayesian polygenic score Probability Conversion to
a correctly centred liability-scale PRS. Follow the detailed
[`LIABILITY_PRS_GUIDE.md`](LIABILITY_PRS_GUIDE.md).

### GenoPred

`genopred_probability()` applies the quantile conversion described by Pain et
al. The legacy name `pain_probability()` is retained for backward
compatibility.
It requires target and population-reference PRSs, observed-scale R², `K`, `SP`
and a justified quantile specification.

## Common high-impact variants

`high_impact_prevalence_parameters()` derives genotype-specific population
prevalence and sample-prevalence values from case/control genotype frequencies.
`high_impact_method_probability()` then recalculates BPC, GenoPred, PAIR
(summary), or PAIR (sample) using those values. This is the canonical TIGER
high-impact-variant update. The main `tiger_probabilities()` interface exposes
this through `include_high_impact`, `high_impact_col`, and
`high_impact_reference`. Before using this route, users must remove the
separately modelled variant and an appropriate LD-aware region from PRS
construction and R² estimation, then document the excluded interval or variant
set. TIGER cannot verify this preprocessing step from a probability input table.
`apply_high_impact_probability()` remains available as a lower-level direct
likelihood-ratio update when that alternative model is explicitly intended.
`apoe_genotype_reference()` supports APOE allele frequencies under an explicit
HWE assumption. `biallelic_genotype_reference()` provides the corresponding
0/1/2-effect-allele reference for one biallelic common high-impact variant.
Observed genotype frequencies are preferred. See
[`AD_APOE_GUIDE.md`](AD_APOE_GUIDE.md) and
[`HIGH_IMPACT_VARIANTS.md`](HIGH_IMPACT_VARIANTS.md).

## Rare variants

`intrinsic_rv_probability()` converts an odds ratio to an intrinsic damaging or
protective probability. `apply_rv_carriers()` combines zero, one or multiple
carried RVs and applies them to a PRS or PRS-plus-APOE probability. This
combination assumes independent effects.

For a damaging effect with `OR >= 1`, TIGER uses

```text
p_damaging = prevalence * (OR - 1) / {1 + prevalence * (OR - 1)}.
```

For a protective effect with `OR < 1`, TIGER applies the inverted effect to the
complementary disease-free probability:

```text
p_protective = (1 - prevalence) * (1 - OR) /
               {1 - prevalence * (1 - OR)}.
```

The individual update is

```text
P_after = (1 - p_protective) *
          {1 - (1 - P_before) * (1 - p_damaging)}.
```

Damaging RVs therefore act on the remaining disease-free probability and
protective RVs multiply the result downward. Same-direction probabilities from
multiple carried RVs are combined as `1 - product(1 - p_j)`.

The RV functions default to `prevalence = 0.50` (`SP_RV = 0.50`). This defines
the balanced 50:50 case/control sample calculation and gives
a sample-level overall RV probability. If the reference frequencies/effects
and intended output are population-level, supply a justified population-level
background probability instead. Do not mix sample-level case/control evidence
with a population-level interpretation without an explicit adaptation.

## Canonical calculation order

The audited TIGER workflow uses the following order:

1. Convert the liability-scale PRS with one selected method to obtain `P_PRS`.
2. If a separately modelled common high-impact genotype is included, calculate
   genotype-specific `K_g` and `SP_g` and rerun that same conversion to obtain
   `P_PRS_high_impact`.
3. Apply the same independently specified RV operator either to `P_PRS` or to
   `P_PRS_high_impact`.

This produces directly comparable PRS, PRS + RV, PRS + high-impact, and PRS +
high-impact + RV probabilities. Do not apply both the method-specific APOE
recalculation and `apply_high_impact_probability()` to the same person. They are
alternative high-impact updates. Use the same RV reference, carrier matrix, and
background prevalence when comparing RV effects before and after the
high-impact update.

See [`DATA.md`](DATA.md) for all input formats and [`REFERENCES.md`](REFERENCES.md)
for the underlying publications.

Run `examples/basic_example.R` for a compact comparison of the available PRS
probability conversions. Use `examples/example_data_and_plots.R` for the full
PRS, APOE and RV workflow.
