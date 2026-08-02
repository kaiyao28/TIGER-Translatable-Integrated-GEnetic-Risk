# TIGER: Translatable Integrated Genetic Risk framework.
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
