# Modelling one common high-impact variant

TIGER can model a single biallelic common variant separately from the PRS. The
variant is represented by three mutually exclusive genotypes containing zero,
one, or two copies of the effect allele. This is the same high-impact update
used for APOE, but with three rather than six genotype categories.

## Reference input

Observed genotype frequencies among comparable cases and controls are
preferred:

```r
single_variant_reference <- prepare_high_impact_reference(data.frame(
  High_impact_genotype = c("0", "1", "2"),
  Case_freq = c(0.64, 0.32, 0.04),
  Control_freq = c(0.81, 0.18, 0.01)
))
```

Each frequency column must sum to one. Genotype labels may instead be study
labels such as `AA`, `AG`, and `GG`, provided the individual data use exactly
the same labels.

The reference may also describe mutually exclusive categories jointly defined
by multiple linked or interacting variants, as in the six-genotype APOE
example. In that setting, provide one category per individual and one reference
row per possible category. The case and control frequencies must describe the
joint categories directly and each frequency column must still sum to one.

If only effect-allele frequencies in cases and controls are available, TIGER
can construct the three genotype frequencies:

```r
single_variant_reference <- biallelic_genotype_reference(
  case_effect_allele_frequency = 0.20,
  control_effect_allele_frequency = 0.10,
  genotype_labels = c("0", "1", "2")
)
```

This helper assumes Hardy-Weinberg equilibrium separately in cases and
controls. Do not use it when that assumption is unsuitable. A marginal
per-allele odds ratio alone does not remove the need to justify how genotype
frequencies are obtained.

## Separate population and sample updates

For genotype `g`, TIGER derives:

```text
K_g  = K  f_case,g / {K  f_case,g + (1 - K)  f_control,g}
SP_g = SP f_case,g / {SP f_case,g + (1 - SP) f_control,g}
```

`K_g` is the genotype-specific population prevalence. `SP_g` is the
genotype-specific case proportion in the target sample. Keeping these terms
separate allows the same BPC, GenoPred, PAIR (summary), or PAIR (sample)
conversion to be recalculated for each genotype without treating sample
ascertainment as population prevalence.

Inspect the derived parameters when preparing a new reference:

```r
high_impact_prevalence_parameters(
  single_variant_reference,
  K = 0.01,
  SP = 0.50
)
```

## Probability calculation

Before calculation, exclude the separately modelled variant and a justified
LD-aware region from PRS construction and R² estimation. Record the genome
build, excluded interval or variant set, and apply the same exclusion to target
and reference PRSs. TIGER cannot check this preprocessing step from the final
PRS table. Supply one genotype per individual:

```r
individuals <- data.frame(
  ID = c("P001", "P002", "P003"),
  PRS_liability = c(-0.4, 0.1, 0.6),
  Variant_genotype = c("0", "1", "2")
)

results <- tiger_probabilities(
  data = individuals,
  K = 0.01,
  SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = 0.10,
  include_high_impact = TRUE,
  high_impact_col = "Variant_genotype",
  high_impact_reference = single_variant_reference
)

results[c("ID", "Probability_PRS", "Probability_PRS_HIGH_IMPACT")]
```

Set `include_rv = TRUE` and provide the usual RV inputs to additionally obtain
`Probability_PRS_HIGH_IMPACT_RV`. Use the separate `include_apoe` route for the
six-genotype APOE model; it returns `Probability_PRS_APOE` and, with RVs,
`Probability_PRS_APOE_RV`.

## Assumptions

- Case and control frequencies must refer to comparable phenotype definitions,
  ancestry, age context, and variant coding.
- All variants defining the high-impact categories must be excluded from the
  PRS and its R² estimate to prevent double counting.
- HWE is an explicit assumption only when allele frequencies are converted to
  genotype frequencies.
- External calibration remains necessary before clinical interpretation.
