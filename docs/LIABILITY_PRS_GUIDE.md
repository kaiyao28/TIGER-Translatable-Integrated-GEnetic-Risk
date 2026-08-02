# Preparing liability-scale PRS inputs

Liability-scale PRS preparation is required before using any TIGER probability
conversion method. BPC, GenoPred, PAIR (summary), and PAIR (sample) all expect
the target PRS supplied through `prs_liability` to be on the liability scale.
Do not pass an unconverted PLINK score, log-odds score, or independently
standardised PRS.

## Required inputs

Prepare two sets of scores:

1. target-sample PRSs for the individuals whose probabilities will be
   estimated; and
2. PRSs from an ancestry-matched population reference sample.

Both must use the same SNP set, effect alleles, effect estimates, and
population-reference allele frequencies for centring. The reference scores are
used to obtain the matching liability-scale PRS variance explained.

You must also specify:

- `K`, the population disorder prevalence; and
- `SP_score`, the ascertainment represented during PRS construction and used
  for the liability-scale conversion.

`SP_score` is passed through the function argument named `SP`. It is distinct
from the target-context `SP` later used to calibrate predicted probabilities.
They may have the same value, such as 0.50 in a balanced simulation, but should
be chosen for their separate purposes.

## 1. Construct the PRS consistently

Construct the PRS in the usual way. Use PRS-CS, SBayesR, or another validated
shrinkage method to estimate posterior SNP effect sizes from the GWAS summary
statistics. These are the same posterior weights normally required to calculate
a PRS. TIGER does not require a separate set of SNP effects.

Where required by the shrinkage method, calculate the effective GWAS sample
size. For one case-control cohort:

For one case-control cohort:

```r
neff <- effective_sample_size(n_cases, n_controls)
```

which evaluates:

```text
Neff = 4 / (1 / N_cases + 1 / N_controls)
```

For a meta-analysis, use the appropriate sum of cohort-specific effective
sample sizes. PRS-CS can use the GWAS log-odds effects and `Neff`. For an
SBayesR workflow requiring standardised observed-scale inputs, TIGER provides:

```r
prepared_sumstats <- prepare_sbayesr_summary_statistics(
  beta = sumstats$beta,
  standard_error = sumstats$se,
  n_effective = sumstats$neff
)
```

Use the selected PRS software as usual, then apply the resulting posterior SNP
effect sizes unchanged to both the target and population-reference samples.

## 2. Centre target and reference scores identically

Use population-reference allele frequencies when scoring both datasets. A
schematic PLINK command is:

```bash
plink1.9 \
  --bfile TARGET_PREFIX \
  --read-freq REFERENCE_ALLELE_FREQUENCIES.frq \
  --score POSTERIOR_EFFECTS.txt VARIANT_COLUMN ALLELE_COLUMN BETA_COLUMN sum center \
  --out TARGET_PRS
```

Run the equivalent command for the population reference. Do not select SNPs or
centre scores independently in the two datasets. Differences in SNP inclusion,
allele direction, weights, or centring change both individual PRSs and the
reference variance.

## 3. Convert both PRSs to the liability scale

When the observed-scale score was constructed under 50% ascertainment, use
`SP_score = 0.5`. Replace `K` with the justified population prevalence:

```r
SP_score <- 0.5

converted <- prepare_liability_prs_inputs(
  target_prs_observed = target_prs_observed,
  reference_prs_observed = reference_prs_observed,
  K = 0.01,
  SP = SP_score,
  center_on_reference = FALSE
)

target_prs_liability <- converted$target_prs_liability
reference_prs_liability <- converted$reference_prs_liability
r2_liability <- converted$r2_liability
```

Use `center_on_reference = FALSE` when the scores were already centred using
the same population-reference allele frequencies. Set it to `TRUE` only when
the inputs were not already centred and subtraction of the empirical reference
mean is the intended procedure.

The scale factor can also be applied explicitly:

```r
scale_factor <- liability_prs_scale_factor(K = 0.01, SP = SP_score)
target_prs_liability <- target_prs_observed * scale_factor
reference_prs_liability <- reference_prs_observed * scale_factor
r2_liability <- var(reference_prs_liability)
```

The converted target PRS and reference-derived `r2_liability` belong together
and must not be mixed with scores or R² estimates produced using a different
model, ancestry, SNP set, or centring procedure.

## Alternative when reference-derived R² is unavailable

If liability-scale R² cannot be estimated from an appropriate population
reference, a compatible reported observed-scale R² may be used. A leave-one-out
or independently evaluated PRS estimate is preferable to an in-sample estimate
because it reduces overfitting bias.

Convert the reported estimate using the population prevalence and the case
proportion of the sample in which R² was evaluated:

```r
K <- 0.01
SP_evaluation <- 0.50
reported_r2_observed <- 0.03

r2_liability <- observed_to_liability_r2(
  r2_observed = reported_r2_observed,
  K = K,
  SP = SP_evaluation
)
```

This route is appropriate only when the reported score is sufficiently
comparable to the target PRS in ancestry, variants, weights, phenotype, and
evaluation design. Record the source of R² and the `K` and `SP_evaluation`
values used for conversion. The resulting `r2_liability` can be supplied to
BPC, PAIR (summary), and the liability-based GenoPred workflow. GenoPred obtains
its required observed-scale value with:

```r
r2_observed <- liability_to_observed_r2(
  r2_liability = r2_liability,
  K = K,
  SP = SP
)
```

Here, the final `SP` is the target-context probability prior. It is not
necessarily the same as `SP_evaluation` from the study that reported R².
This alternative supplies R² only. GenoPred still requires a consistently
constructed liability-scale population-reference PRS distribution.

## 4. Use the prepared inputs across TIGER methods

```r
# BPC
p_bpc <- bpc_probability(
  target_prs_liability, K, SP, r2_liability
)

# GenoPred
r2_observed <- liability_to_observed_r2(r2_liability, K, SP)
p_genopred <- genopred_probability(
  target_prs_liability,
  reference_prs_liability,
  r2_observed,
  K,
  SP
)

# PAIR (summary)
p_pair_summary <- pair_probability_summary(
  target_prs_liability, K, SP, r2_liability
)
```

PAIR (sample) also requires liability-scale target PRSs. Its case and control
means and standard deviations must be estimated from PRSs constructed and
scaled in the same way, preferably in an independent calibration sample.

## Quality-control checklist

- Same SNP set, effect alleles, and weights in target and reference samples.
- Same population-reference allele frequencies used for centring.
- Correct effective discovery sample size for the chosen PRS method.
- Correct population prevalence `K`, conversion ascertainment `SP_score`, and
  target probability prior `SP`.
- Liability transformation applied identically to target and reference PRSs.
- Reference-derived liability-scale R² is finite, positive, and plausible.
- If reported observed-scale R² is used, its source, evaluation design,
  `K`, and `SP_evaluation` are documented.
- Target, discovery, and population-reference ancestries are appropriately
  matched.
- Converted PRS values are not independently re-standardised afterward.

Run [`examples/liability_conversion_example.R`](../examples/liability_conversion_example.R)
for a complete synthetic example.
