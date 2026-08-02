# TIGER framework. GNU GPL v3 or later. See LICENSE.
# Minimal observed-to-liability PRS conversion example.
# Run from the TIGER repository root:
# Rscript examples/liability_conversion_example.R

source("R/tiger.R")
load_tiger()

target <- read.csv(
  tiger_example_file("example_target_prs_observed.csv"),
  stringsAsFactors = FALSE
)
reference <- read.csv(
  tiger_example_file("example_reference_prs_observed.csv"),
  stringsAsFactors = FALSE
)

K <- 0.01
P <- 0.50

# These example scores are treated as already centred with the same
# population-reference allele frequencies. In a real analysis, target and
# reference scores must use the same SNPs, alleles, posterior effects and
# centring procedure.
converted <- prepare_bpc_inputs(
  target_prs_observed = target$PRS_observed,
  reference_prs_observed = reference$PRS_observed,
  K = K,
  P = P,
  center_on_reference = FALSE
)

target$PRS_liability <- converted$target_prs_liability
reference$PRS_liability <- converted$reference_prs_liability
r2_liability <- converted$r2_liability

cat("Liability conversion multiplier:",
    round(converted$liability_scale_factor, 6), "\n")
cat("Reference-derived liability-scale PRS R2:",
    round(r2_liability, 6), "\n\n")
print(target)

# The converted target PRS and matching reference-derived R2 can now be used in
# PAIR (summary), the default TIGER worked-example method.
target$Probability_PAIR_summary <- pair_probability_summary(
  prs_liability = target$PRS_liability,
  K = K,
  prior = P,
  r2_liability = r2_liability
)

cat("\nPAIR (summary) probabilities:\n")
print(target)
