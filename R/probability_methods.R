# PRS and rare-variant probability toolkit.
# Copyright (C) 2026 TIGER study authors
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed WITHOUT ANY WARRANTY. See LICENSE.

.check_probability <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop(name, " must be one finite value strictly between 0 and 1")
  }
}

.check_r2 <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 0 || x >= 1) {
    stop(name, " must be one finite value in [0, 1)")
  }
}

# Effective sample size for a case/control GWAS or meta-analysis cohort.
effective_sample_size <- function(n_cases, n_controls) {
  if (length(n_cases) != length(n_controls) ||
      any(!is.finite(n_cases)) || any(!is.finite(n_controls)) ||
      any(n_cases <= 0) || any(n_controls <= 0)) {
    stop("n_cases and n_controls must be equal-length positive values")
  }
  4 / (1 / n_cases + 1 / n_controls)
}

# Liability multiplier for a PRS on the standardised observed scale with
# ascertainment SP. SP = 0.5 is used when the score construction represents a
# balanced case-control setting through effective sample size.
liability_prs_scale_factor <- function(K, SP = 0.5) {
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  threshold <- -stats::qnorm(K)
  density <- stats::dnorm(threshold)
  sqrt(K^2 * (1 - K)^2 / (density^2 * SP * (1 - SP)))
}

# Convert observed-scale target and reference PRSs to liability scale and
# estimate R2_liability as the variance of the reference PRS. Scores must have
# been constructed with the same SNPs, posterior effects and reference-allele
# centring. Set center_on_reference only if scores were not already centred.
prepare_liability_prs_inputs <- function(
    target_prs_observed, reference_prs_observed,
    K, SP = 0.5, center_on_reference = FALSE) {
  if (!is.numeric(target_prs_observed) ||
      any(!is.finite(target_prs_observed)) ||
      !length(target_prs_observed)) stop("target PRS must be finite numeric")
  if (!is.numeric(reference_prs_observed) ||
      any(!is.finite(reference_prs_observed)) ||
      length(reference_prs_observed) < 2L) {
    stop("reference PRS must contain at least two finite values")
  }
  if (isTRUE(center_on_reference)) {
    reference_mean <- mean(reference_prs_observed)
    target_prs_observed <- target_prs_observed - reference_mean
    reference_prs_observed <- reference_prs_observed - reference_mean
  }
  multiplier <- liability_prs_scale_factor(K, SP)
  target_liability <- target_prs_observed * multiplier
  reference_liability <- reference_prs_observed * multiplier
  list(
    target_prs_liability = target_liability,
    reference_prs_liability = reference_liability,
    r2_liability = stats::var(reference_liability),
    liability_scale_factor = multiplier
  )
}

# Prepare log-odds GWAS effects and standard errors for SBayesR on the
# standardised observed scale with 50% ascertainment, following BPC guidance.
prepare_sbayesr_summary_statistics <- function(beta, standard_error,
                                               n_effective) {
  lengths <- c(length(beta), length(standard_error), length(n_effective))
  target_length <- max(lengths)
  if (any(!lengths %in% c(1L, target_length))) {
    stop("inputs must have equal lengths or be scalar")
  }
  if (any(!is.finite(beta)) || any(!is.finite(standard_error)) ||
      any(standard_error <= 0) || any(!is.finite(n_effective)) ||
      any(n_effective <= 0)) stop("invalid summary-statistic inputs")
  data.frame(
    beta_50_50 = beta / (standard_error * sqrt(n_effective)),
    standard_error_50_50 = 1 / sqrt(n_effective)
  )
}

# Convert observed-scale variance explained to the liability scale using the
# Lee et al. observed-to-liability transformation.
observed_to_liability_r2 <- function(r2_observed, K, SP = 0.5) {
  .check_r2(r2_observed, "r2_observed")
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  threshold <- -stats::qnorm(K)
  density <- stats::dnorm(threshold)
  r2_observed * K^2 * (1 - K)^2 /
    (density^2 * SP * (1 - SP))
}

# Exact algebraic inverse of observed_to_liability_r2().
liability_to_observed_r2 <- function(r2_liability, K, SP = 0.5) {
  .check_r2(r2_liability, "r2_liability")
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  threshold <- -stats::qnorm(K)
  density <- stats::dnorm(threshold)
  r2_liability * density^2 * SP * (1 - SP) /
    (K^2 * (1 - K)^2)
}

# Bayesian polygenic score Probability Conversion (BPC).
#
# prs_liability: PRS values on the liability scale, centred using an
# ancestry-matched population reference sample.
# K: population lifetime prevalence.
# SP: sample prevalence in the target/test context.
# r2_liability: liability-scale PRS variance explained.
bpc_probability <- function(prs_liability, K, SP, r2_liability) {
  if (!is.numeric(prs_liability) || !length(prs_liability) ||
      any(!is.finite(prs_liability))) stop("prs_liability must be finite numeric")
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  .check_r2(r2_liability, "r2_liability")
  if (r2_liability <= 0) return(rep(SP, length(prs_liability)))

  threshold <- -stats::qnorm(K)
  density <- stats::dnorm(threshold)
  i_case <- density / K
  k_case <- i_case * (i_case - threshold)
  i_control <- -density / (1 - K)
  k_control <- i_control * (i_control - threshold)
  mean_case <- i_case * r2_liability
  mean_control <- i_control * r2_liability
  variance_case <- r2_liability - k_case * r2_liability^2
  variance_control <- r2_liability - k_control * r2_liability^2
  if (variance_case <= 0 || variance_control <= 0) {
    stop("r2_liability produces an invalid theoretical PRS variance")
  }

  # Work on the log-odds scale so extreme PRS values do not cause both normal
  # densities to underflow to zero and produce NaN.
  log_case_density <- stats::dnorm(
    prs_liability, mean_case, sqrt(variance_case), log = TRUE
  )
  log_control_density <- stats::dnorm(
    prs_liability, mean_control, sqrt(variance_control), log = TRUE
  )
  stats::plogis(
    stats::qlogis(SP) + log_case_density - log_control_density
  )
}

# Pain et al. quantile conversion with an optional population-prevalence
# adaptation. Set corrected = FALSE to return the original SP-based version.
genopred_probability <- function(prs_liability, reference_prs_liability,
                                 r2_observed, K, SP = 0.5,
                                 n_quantiles = 100, corrected = TRUE) {
  if (!is.numeric(prs_liability) || !length(prs_liability) ||
      any(!is.finite(prs_liability))) stop("prs_liability must be finite numeric")
  if (!is.numeric(reference_prs_liability) ||
      sum(is.finite(reference_prs_liability)) < 2L) {
    stop("reference_prs_liability must contain at least two finite values")
  }
  .check_r2(r2_observed, "r2_observed")
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  if (length(n_quantiles) != 1L || !is.finite(n_quantiles)) {
    stop("n_quantiles must be one finite integer of at least 2")
  }
  n_quantiles <- as.integer(n_quantiles)
  if (n_quantiles < 2L) stop("n_quantiles must be at least 2")
  if (length(corrected) != 1L || is.na(corrected) || !is.logical(corrected)) {
    stop("corrected must be TRUE or FALSE")
  }
  if (r2_observed <= 1e-4) return(rep(SP, length(prs_liability)))

  ref_mean <- mean(reference_prs_liability, na.rm = TRUE)
  ref_sd <- stats::sd(reference_prs_liability, na.rm = TRUE)
  if (!is.finite(ref_sd) || ref_sd <= 0) stop("reference PRS has zero variance")
  z_prs <- (prs_liability - ref_mean) / ref_sd

  discovery_prevalence <- 0.5
  alpha <- 1 / (discovery_prevalence * (1 - discovery_prevalence))
  d <- sqrt(alpha) * sqrt(r2_observed) / sqrt(1 - r2_observed)
  bin_mass <- 1 / n_quantiles
  quantile_probabilities <- seq(
    bin_mass, 1 - bin_mass, by = bin_mass
  )
  mixture_quantile <- function(pq) {
    stats::uniroot(
      function(x) SP * stats::pnorm(x - d) +
        (1 - SP) * stats::pnorm(x) - pq,
      interval = c(-2.5, 2.5), extendInt = "yes", tol = 6e-12
    )$root
  }
  raw_breaks <- vapply(quantile_probabilities, mixture_quantile, numeric(1))
  moment_prevalence <- if (isTRUE(corrected)) K else SP
  variance_prs <- moment_prevalence *
    (1 + d^2 - (d * moment_prevalence)^2) +
    (1 - moment_prevalence) * (1 - (d * moment_prevalence)^2)
  mean_prs <- d * moment_prevalence
  raw_bounds <- c(-Inf, raw_breaks, Inf)
  standardised_bounds <- (raw_bounds - mean_prs) / sqrt(variance_prs)
  case_mass <- stats::pnorm(raw_bounds[-1], mean = d) -
    stats::pnorm(raw_bounds[-length(raw_bounds)], mean = d)
  p_case <- case_mass * SP / bin_mass
  bin <- findInterval(z_prs, standardised_bounds[-1]) + 1L
  bin <- pmin(pmax(bin, 1L), n_quantiles)
  pmin(pmax(p_case[bin], 0), 1)
}

# PAIR conversion using supplied case/control PRS moments. population_prevalence
# reconstructs population PRS moments; SP calibrates the output probability
# to the intended target context.
.pair_probability_from_moments <- function(
    prs_liability, case_mean, case_sd, control_mean, control_sd,
    population_prevalence, SP) {
  if (!is.numeric(prs_liability) || !length(prs_liability) ||
      any(!is.finite(prs_liability))) stop("prs_liability must be finite numeric")
  moments <- c(case_mean, case_sd, control_mean, control_sd)
  if (any(lengths(list(case_mean, case_sd, control_mean, control_sd)) != 1L) ||
      any(!is.finite(moments)) || case_sd <= 0 || control_sd <= 0) {
    stop("case/control means must be finite scalars and SDs must be positive")
  }
  .check_probability(population_prevalence, "population_prevalence")
  .check_probability(SP, "SP")

  mixture_mean <- population_prevalence * case_mean +
    (1 - population_prevalence) * control_mean
  mixture_variance <- population_prevalence * case_sd^2 +
    (1 - population_prevalence) * control_sd^2 +
    population_prevalence * (1 - population_prevalence) *
      (case_mean - control_mean)^2
  mean_difference_sq <- (case_mean - control_mean)^2
  r0 <- (case_sd^2 + mean_difference_sq) / control_sd^2
  r1 <- (control_sd^2 + mean_difference_sq) / case_sd^2
  beta <- (case_mean - control_mean) *
    (SP * (1 - SP) * ((r0 + r1) / 2 - 1) +
       SP * case_sd^2 / control_sd^2 +
       (1 - SP) * control_sd^2 / case_sd^2) /
    mixture_variance
  alpha <- log(SP * control_sd / ((1 - SP) * case_sd)) +
    ((r0 - 1) * SP + (1 - r1) * (1 - SP)) / 2 -
    mixture_mean * beta
  stats::plogis(alpha + beta * prs_liability)
}

# PAIR (summary): phenotype-independent conversion using theoretical
# case/control PRS moments implied by population K and liability-scale PRS R2.
pair_probability_summary <- function(prs_liability, K, SP = 0.5,
                                     r2_liability) {
  if (!is.numeric(prs_liability) || !length(prs_liability) ||
      any(!is.finite(prs_liability))) stop("prs_liability must be finite numeric")
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  .check_r2(r2_liability, "r2_liability")
  if (r2_liability <= 0) return(rep(SP, length(prs_liability)))

  threshold <- -stats::qnorm(K)
  density <- stats::dnorm(threshold)
  i_case <- density / K
  k_case <- i_case * (i_case - threshold)
  i_control <- -density / (1 - K)
  k_control <- i_control * (i_control - threshold)
  m1 <- i_case * r2_liability
  m0 <- i_control * r2_liability
  variance_case <- r2_liability - k_case * r2_liability^2
  variance_control <- r2_liability - k_control * r2_liability^2
  if (variance_case <= 0 || variance_control <= 0) {
    stop("r2_liability produces an invalid theoretical PRS variance")
  }
  s1 <- sqrt(variance_case)
  s0 <- sqrt(variance_control)

  .pair_probability_from_moments(
    prs_liability, m1, s1, m0, s0,
    population_prevalence = K, SP = SP
  )
}

# PAIR (sample): conversion using case/control PRS moments estimated in a
# labelled sample. For prospective application, estimate the moments in a
# separate representative calibration dataset, not in the people being scored.
pair_probability_sample <- function(
    prs_liability, case_mean, case_sd, control_mean, control_sd,
    K, SP = 0.5) {
  .pair_probability_from_moments(
    prs_liability, case_mean, case_sd, control_mean, control_sd,
    population_prevalence = K, SP = SP
  )
}
