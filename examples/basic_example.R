# TIGER framework. GNU GPL v3 or later. See LICENSE.
# Minimal worked example. Run from the repository root with:
# Rscript examples/basic_example.R

source("R/tiger.R")
load_tiger()

set.seed(2026)
K <- 0.01                  # population lifetime prevalence
SP <- 0.50              # sample prevalence in the target/test context
r2_liability_target <- 0.10

# In practice these observed-scale values come from the same posterior mean
# effects, SNP set and population-reference centring procedure. Synthetic values
# are used here only to demonstrate the application functions.
r2_observed <- liability_to_observed_r2(r2_liability_target, K, SP)
reference_prs_observed <- stats::rnorm(500, 0, sqrt(r2_observed))
target_prs_observed <- seq(-0.8, 0.8, length.out = 9) /
  bpc_liability_scale_factor(K, SP = 0.5)
prepared <- prepare_bpc_inputs(
  target_prs_observed, reference_prs_observed,
  K = K, SP = 0.5, center_on_reference = FALSE
)
target_prs <- prepared$target_prs_liability
reference_prs <- prepared$reference_prs_liability
r2_liability <- prepared$r2_liability

probabilities <- data.frame(
  PRS_liability = target_prs,
  BPC = bpc_probability(target_prs, K, SP, r2_liability),
  GenoPred = genopred_probability(
    target_prs, reference_prs, r2_observed, K, SP = SP,
    corrected = TRUE
  ),
  PAIR_summary = pair_probability_summary(
    target_prs, K, SP, r2_liability
  )
)

# Example damaging RV with OR 8 and protective RV with OR 0.4. SP = 0.50 is the
# balanced-sample RV default; change it only for a justified application level.
p_damaging <- intrinsic_rv_probability(8, prevalence = SP)
p_protective <- intrinsic_rv_probability(0.4, prevalence = SP)
probabilities$BPC_plus_damaging_RV <- apply_rv_probability(
  probabilities$BPC, p_damaging = p_damaging
)
probabilities$BPC_plus_protective_RV <- apply_rv_probability(
  probabilities$BPC, p_protective = p_protective
)

# Common high-impact example: APOE frequencies are illustrative only. In an AD
# application, exclude the separately modelled APOE region from the PRS and
# replace these values with justified population-matched evidence.
apoe_reference <- apoe_genotype_reference(
  case_allele_frequencies = c(e2 = 0.04, e3 = 0.66, e4 = 0.30),
  control_allele_frequencies = c(e2 = 0.08, e3 = 0.79, e4 = 0.13)
)
apoe <- rep(c("e3/e3", "e3/e4", "e4/e4"), length.out = nrow(probabilities))
probabilities$BPC_plus_APOE <- high_impact_method_probability(
  target_prs, apoe, apoe_reference,
  K = K, SP = SP, method = "BPC", r2_liability = r2_liability
)
probabilities$BPC_plus_APOE_plus_RV <- apply_rv_probability(
  probabilities$BPC_plus_APOE, p_damaging = p_damaging
)

print(round(probabilities, 4))
