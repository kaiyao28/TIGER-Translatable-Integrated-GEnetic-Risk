# Preparing liability-scale PRS inputs for BPC

This guide follows the BPC authors' application instructions. Consult the
authors' current repository and publication for authoritative software-version
and empirical-analysis details.

Official BPC v1.0.1 software archive:
https://doi.org/10.5281/zenodo.15721085

Liability-scale preparation is a prerequisite for the TIGER examples. Do not
pass an unconverted raw PLINK score, log-odds score or arbitrarily standardized
PRS to functions whose `prs_liability` argument requires the liability scale.

## Required BPC inputs

`bpc_probability()` requires:

1. `prs_liability`: one or more Bayesian PRS values on the liability scale;
2. `K`: population disorder prevalence;
3. `prior`: prior disorder probability in the target context;
4. `r2_liability`: PRS variance explained on the liability scale.

The same SNP set must be present in the GWAS/PRS model, target cohort and
population reference sample. Alleles and effect directions must be harmonised.

## 1. Effective sample size

For one case/control cohort:

```r
neff <- effective_sample_size(n_cases, n_controls)
```

which evaluates

```text
Neff = 4 / (1 / N_cases + 1 / N_controls).
```

For a GWAS meta-analysis, the BPC guidance recommends the sum of cohort-specific
effective sample sizes where those cohort values are available.

## 2. Posterior mean effects

### PRS-CS

Use GWAS log-odds betas and `Neff` as input to PRS-CS. The use of `Neff`
represents the effects on a standardised observed scale with 50% ascertainment.
The BPC publication evaluated PRS-CS and SBayesR; other Bayesian methods require
their own validation.

### SBayesR

Before running SBayesR, transform log-odds betas and standard errors:

```r
prepared <- prepare_sbayesr_summary_statistics(
  beta = sumstats$beta,
  standard_error = sumstats$se,
  n_effective = sumstats$neff
)

sumstats$beta_50_50 <- prepared$beta_50_50
sumstats$se_50_50 <- prepared$standard_error_50_50
```

The equations are:

```text
beta_50_50 = beta / (SE × sqrt(Neff))
SE_50_50   = 1 / sqrt(Neff).
```

Supply these transformed statistics and `Neff` to SBayesR.

## 3. Calculate identically centred PRSs

Apply the same posterior mean effects to the target individuals and an
ancestry-matched population reference sample such as 1000 Genomes.

The BPC authors use PLINK 1.9 and population-reference allele frequencies to
centre scores consistently. A schematic command is:

```bash
plink1.9 \
  --bfile TARGET_PREFIX \
  --read-freq REFERENCE_ALLELE_FREQUENCIES.frq \
  --score POSTERIOR_EFFECTS.txt VARIANT_COLUMN ALLELE_COLUMN BETA_COLUMN sum center \
  --out TARGET_PRS
```

Run the corresponding command for the population reference using the same:

- SNP identifiers and SNP set;
- effect alleles;
- posterior mean effects;
- allele-frequency reference used for centring.

Do not independently select SNPs in target and reference data. Missing or
mismatched SNPs can alter both individual PRS values and the reference variance.

## 4. Convert PRSs to the liability scale

When the Bayesian method used `Neff`, BPC treats the PRS as being on the
standardised observed scale with `P = 0.5`. Convert both target and reference
PRSs using:

```r
scale_factor <- bpc_liability_scale_factor(K = 0.01, P = 0.5)
target_prs_liability <- target_prs_observed * scale_factor
reference_prs_liability <- reference_prs_observed * scale_factor
```

Alternatively, prepare both inputs together:

```r
prepared <- prepare_bpc_inputs(
  target_prs_observed = target_prs_observed,
  reference_prs_observed = reference_prs_observed,
  K = 0.01,
  P = 0.5,
  center_on_reference = FALSE
)
```

Use `center_on_reference = FALSE` when PLINK already centred scores using the
population-reference allele frequencies. Set it to `TRUE` only when input scores
have not already been centred and subtracting the empirical reference mean is
the intended centring procedure.

## 5. Estimate liability-scale R²

BPC estimates liability-scale R² without reference phenotypes:

```r
r2_liability <- var(reference_prs_liability)
```

`prepare_bpc_inputs()` returns the same estimate as
`prepared$r2_liability`. The reference sample must be ancestry-matched and use
the same score construction as the target sample.

## 6. Calculate BPC probabilities

```r
probability <- bpc_probability(
  prs_liability = prepared$target_prs_liability,
  K = 0.01,
  prior = 0.5,
  r2_liability = prepared$r2_liability
)
```

For a randomly selected population member, the prior may equal `K`. In a
help-seeking or otherwise selected target setting, the appropriate prior may be
higher and should be justified independently of genotype.

## Quality-control checklist

- Same SNP set in model, target and reference data.
- Alleles harmonised and ambiguous mismatches resolved.
- Same posterior mean effects in target and reference PRSs.
- Same allele-frequency reference used for centring.
- Correct `Neff`, including cohort-specific aggregation for meta-analysis.
- Correct population prevalence `K` for the intended disorder and population.
- Liability transformation applied to both target and reference PRSs.
- Reference PRS variance finite, positive and plausible.
- Target and GWAS/reference ancestries appropriately matched.
- Prior defined for the target context and reported transparently.

The BPC authors report evaluation in European-ancestry individuals. Applications
to other ancestries require appropriate reference data and independent
calibration assessment.
