# TIGER: Translatable Integrated Genetic Risk.
# Copyright (C) 2026 TIGER study authors
# Licensed under GNU GPL v3 or later; distributed WITHOUT ANY WARRANTY.

# Validate a case/control reference table for a separately modelled common
# high-impact genotype. Frequencies are conditional on disease status and must
# each sum to one across the mutually exclusive genotype categories.
prepare_high_impact_reference <- function(x, tolerance = 1e-8) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  required <- c("Genotype", "Case_freq", "Control_freq")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing column(s): ", paste(missing, collapse = ", "))
  out <- data.frame(
    Genotype = as.character(x$Genotype),
    Case_freq = as.numeric(x$Case_freq),
    Control_freq = as.numeric(x$Control_freq),
    stringsAsFactors = FALSE
  )
  if (!nrow(out) || anyNA(out) || any(!nzchar(out$Genotype)) ||
      anyDuplicated(out$Genotype) || any(!is.finite(out$Case_freq)) ||
      any(!is.finite(out$Control_freq)) || any(out$Case_freq < 0) ||
      any(out$Control_freq < 0)) {
    stop("genotypes must be unique and frequencies must be finite and non-negative")
  }
  if (abs(sum(out$Case_freq) - 1) > tolerance ||
      abs(sum(out$Control_freq) - 1) > tolerance) {
    stop("case and control genotype frequencies must each sum to one")
  }
  out$Likelihood_ratio <- out$Case_freq / out$Control_freq
  out
}

# Convert APOE e2/e3/e4 allele frequencies in cases and controls to the six
# unordered genotype frequencies under Hardy-Weinberg equilibrium (HWE).
# HWE is an explicit modelling assumption; use prepare_high_impact_reference()
# with observed genotype frequencies when those are available.
apoe_genotype_reference <- function(case_allele_frequencies,
                                    control_allele_frequencies) {
  required <- c("e2", "e3", "e4")
  check_alleles <- function(x, name) {
    if (!is.numeric(x) || !identical(sort(names(x)), sort(required)) ||
        any(!is.finite(x)) || any(x < 0) || abs(sum(x) - 1) > 1e-8) {
      stop(name, " must be named e2/e3/e4 frequencies summing to one")
    }
  }
  check_alleles(case_allele_frequencies, "case_allele_frequencies")
  check_alleles(control_allele_frequencies, "control_allele_frequencies")
  genotypes <- data.frame(
    Genotype = c("e2/e2", "e2/e3", "e2/e4", "e3/e3", "e3/e4", "e4/e4"),
    allele_1 = c("e2", "e2", "e2", "e3", "e3", "e4"),
    allele_2 = c("e2", "e3", "e4", "e3", "e4", "e4"),
    stringsAsFactors = FALSE
  )
  multiplier <- ifelse(genotypes$allele_1 == genotypes$allele_2, 1, 2)
  prepare_high_impact_reference(data.frame(
    Genotype = genotypes$Genotype,
    Case_freq = multiplier * case_allele_frequencies[genotypes$allele_1] *
      case_allele_frequencies[genotypes$allele_2],
    Control_freq = multiplier * control_allele_frequencies[genotypes$allele_1] *
      control_allele_frequencies[genotypes$allele_2]
  ))
}

# Derive a three-genotype reference for one biallelic common high-impact
# variant from effect-allele frequencies in cases and controls. HWE is assumed
# separately within each phenotype group. Prefer observed genotype frequencies
# with prepare_high_impact_reference() when they are available.
biallelic_genotype_reference <- function(
    case_effect_allele_frequency,
    control_effect_allele_frequency,
    genotype_labels = c("0", "1", "2")) {
  check_frequency <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || x > 1) {
      stop(name, " must be one finite allele frequency in [0, 1]")
    }
  }
  check_frequency(case_effect_allele_frequency,
                  "case_effect_allele_frequency")
  check_frequency(control_effect_allele_frequency,
                  "control_effect_allele_frequency")
  if (!is.character(genotype_labels) || length(genotype_labels) != 3L ||
      anyNA(genotype_labels) || any(!nzchar(genotype_labels)) ||
      anyDuplicated(genotype_labels)) {
    stop("genotype_labels must contain three unique non-empty labels for 0, 1 and 2 effect alleles")
  }
  hwe <- function(frequency) {
    c((1 - frequency)^2, 2 * frequency * (1 - frequency), frequency^2)
  }
  prepare_high_impact_reference(data.frame(
    Genotype = genotype_labels,
    Case_freq = hwe(case_effect_allele_frequency),
    Control_freq = hwe(control_effect_allele_frequency),
    stringsAsFactors = FALSE
  ))
}

# Update a PRS-derived probability using the case/control likelihood ratio for
# an observed common high-impact genotype. This is a Bayes update and assumes
# the separately modelled genotype is not already represented in the PRS.
apply_high_impact_probability <- function(p_prs, genotype, reference) {
  if (!is.numeric(p_prs) || !length(p_prs) || any(!is.finite(p_prs)) ||
      any(p_prs < 0) || any(p_prs > 1)) {
    stop("p_prs must contain finite probabilities in [0, 1]")
  }
  reference <- prepare_high_impact_reference(reference)
  target_length <- max(length(p_prs), length(genotype))
  if (!length(genotype) ||
      any(!c(length(p_prs), length(genotype)) %in% c(1L, target_length))) {
    stop("p_prs and genotype must have equal lengths or be scalar")
  }
  p_prs <- rep(p_prs, length.out = target_length)
  genotype <- rep(as.character(genotype), length.out = target_length)
  index <- match(genotype, reference$Genotype)
  if (anyNA(index)) {
    stop("Unknown genotype(s): ", paste(unique(genotype[is.na(index)]), collapse = ", "))
  }
  case_frequency <- reference$Case_freq[index]
  control_frequency <- reference$Control_freq[index]
  denominator <- p_prs * case_frequency + (1 - p_prs) * control_frequency
  if (any(denominator <= 0)) stop("genotype has zero probability in both groups")
  p_prs * case_frequency / denominator
}

# Derive genotype-specific population prevalence and target-sample SP from
# case/control genotype frequencies. These are the K_g and sample-prevalence
# quantities used by TIGER when each PRS conversion is recalculated by genotype.
high_impact_prevalence_parameters <- function(reference, K, SP = 0.5) {
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  reference <- prepare_high_impact_reference(reference)
  population_denominator <- K * reference$Case_freq +
    (1 - K) * reference$Control_freq
  sample_denominator <- SP * reference$Case_freq +
    (1 - SP) * reference$Control_freq
  if (any(population_denominator <= 0) || any(sample_denominator <= 0)) {
    stop("each genotype must occur in the population and target mixture")
  }
  reference$K_genotype <-
    K * reference$Case_freq / population_denominator
  reference$SP_genotype <-
    SP * reference$Case_freq / sample_denominator
  if (any(reference$K_genotype <= 0 | reference$K_genotype >= 1) ||
      any(reference$SP_genotype <= 0 | reference$SP_genotype >= 1)) {
    stop(
      "method-specific high-impact updates require positive case and control ",
      "frequency for every genotype"
    )
  }
  reference
}

# Recalculate a selected PRS-to-probability method using genotype-specific K_g
# and sample prevalence. This is TIGER's canonical high-impact-variant update.
high_impact_method_probability <- function(
    prs_liability, genotype, reference, K, SP = 0.5,
    method = c("PAIR (summary)", "PAIR (sample)", "BPC", "GenoPred"),
    r2_liability = NULL, reference_prs_liability = NULL,
    r2_observed = NULL, case_mean = NULL, case_sd = NULL,
    control_mean = NULL, control_sd = NULL, n_quantiles = 100) {
  if (!is.numeric(prs_liability) || !length(prs_liability) ||
      any(!is.finite(prs_liability))) {
    stop("prs_liability must contain finite numeric values")
  }
  method <- match.arg(method)
  genotype <- as.character(genotype)
  target_length <- max(length(prs_liability), length(genotype))
  if (!length(genotype) ||
      any(!c(length(prs_liability), length(genotype)) %in% c(1L, target_length))) {
    stop("prs_liability and genotype must have equal lengths or be scalar")
  }
  prs_liability <- rep(prs_liability, length.out = target_length)
  genotype <- rep(genotype, length.out = target_length)
  parameters <- high_impact_prevalence_parameters(reference, K, SP)
  index <- match(genotype, parameters$Genotype)
  if (anyNA(index)) {
    stop("Unknown genotype(s): ",
         paste(unique(genotype[is.na(index)]), collapse = ", "))
  }
  output <- numeric(target_length)
  for (g in unique(genotype)) {
    selected <- genotype == g
    row <- parameters[parameters$Genotype == g, , drop = FALSE]
    Kg <- row$K_genotype
    SP_g <- row$SP_genotype
    r2_observed_g <- if (method == "GenoPred" && !is.null(r2_liability)) {
      liability_to_observed_r2(r2_liability, Kg, SP_g)
    } else {
      r2_observed
    }
    output[selected] <- switch(
      method,
      "BPC" = bpc_probability(
        prs_liability[selected], Kg, SP_g, r2_liability
      ),
      "GenoPred" = genopred_probability(
        prs_liability[selected], reference_prs_liability,
        r2_observed_g, Kg, SP = SP_g, n_quantiles = n_quantiles
      ),
      "PAIR (summary)" = pair_probability_summary(
        prs_liability[selected], Kg, SP_g, r2_liability
      ),
      "PAIR (sample)" = pair_probability_sample(
        prs_liability[selected], case_mean, case_sd,
        control_mean, control_sd, Kg, SP_g
      )
    )
  }
  output
}
