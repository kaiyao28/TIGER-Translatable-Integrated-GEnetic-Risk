# TIGER framework. GNU GPL v3 or later. See LICENSE.
if (requireNamespace("TIGER", quietly = TRUE)) {
  library(TIGER)
} else {
  source("R/tiger.R")
  loaded_tiger_files <- load_tiger()
  stopifnot(length(loaded_tiger_files) == 6L)
}

K <- 0.01
SP <- 0.5
r2l <- 0.1
prs <- seq(-0.5, 0.5, length.out = 21)
reference <- seq(-1, 1, length.out = 500) * sqrt(r2l)
r2o <- liability_to_observed_r2(r2l, K, SP)
stopifnot(abs(observed_to_liability_r2(r2o, K, SP) - r2l) < 1e-12)
stopifnot(effective_sample_size(100, 100) == 200)
prepared_sbayesr <- prepare_sbayesr_summary_statistics(
  beta = 0.1, standard_error = 0.02, n_effective = 1000
)
stopifnot(
  abs(prepared_sbayesr$beta_50_50 - 0.1 / (0.02 * sqrt(1000))) < 1e-12,
  abs(prepared_sbayesr$standard_error_50_50 - 1 / sqrt(1000)) < 1e-12
)
scale_factor <- bpc_liability_scale_factor(K, SP)
prepared_bpc <- prepare_bpc_inputs(
  target_prs_observed = prs / scale_factor,
  reference_prs_observed = reference / scale_factor,
  K = K, SP = SP
)
stopifnot(
  max(abs(prepared_bpc$target_prs_liability - prs)) < 1e-12,
  abs(prepared_bpc$r2_liability - stats::var(reference)) < 1e-12
)

outputs <- list(
  bpc_probability(prs, K, SP, r2l),
  genopred_probability(prs, reference, r2o, K, SP),
  pair_probability_summary(prs, K, SP, r2l)
)
stopifnot(identical(genopred_probability, pain_probability))
stopifnot(all(vapply(outputs, function(x) {
  length(x) == length(prs) && all(is.finite(x)) &&
    all(x >= 0) && all(x <= 1)
}, logical(1))))
stopifnot(all(vapply(outputs, function(x) all(diff(x) >= 0), logical(1))))
extreme_bpc <- bpc_probability(c(-100, 100), K, SP, r2l)
stopifnot(all(is.finite(extreme_bpc)), all(extreme_bpc >= 0),
          all(extreme_bpc <= 1))

pair_summary <- pair_probability_summary(prs, K, SP, r2l)
stopifnot(identical(pair_probability(prs, K, SP, r2l), pair_summary))
threshold <- -stats::qnorm(K)
density <- stats::dnorm(threshold)
i_case <- density / K
i_control <- -density / (1 - K)
sample_case_mean <- i_case * r2l
sample_control_mean <- i_control * r2l
sample_case_sd <- sqrt(r2l - i_case * (i_case - threshold) * r2l^2)
sample_control_sd <- sqrt(r2l - i_control * (i_control - threshold) * r2l^2)
pair_sample <- pair_probability_sample(
  prs, sample_case_mean, sample_case_sd,
  sample_control_mean, sample_control_sd, K = K, SP = SP
)
stopifnot(
  length(pair_sample) == length(prs),
  all(is.finite(pair_sample)), all(diff(pair_sample) >= 0),
  max(abs(pair_sample - pair_summary)) < 1e-12
)

p_dam <- intrinsic_rv_probability(5, K)
p_pro <- intrinsic_rv_probability(0.5, K)
stopifnot(
  identical(intrinsic_rv_probability(5),
            intrinsic_rv_probability(5, prevalence = 0.5))
)
baseline <- c(0.1, 0.5, 0.9)
stopifnot(all(apply_rv_probability(baseline, p_damaging = p_dam) >= baseline))
stopifnot(all(apply_rv_probability(baseline, p_protective = p_pro) <= baseline))
stopifnot(abs(combine_independent_rv_probabilities(c(0.2, 0.3)) - 0.44) < 1e-12)

reference_table <- data.frame(
  Symbol = c("RISK1", "PRO1"), Case_freq = c(0.01, 0.002),
  Control_freq = c(0.001, 0.004), OR = c(5, 0.5)
)
clean <- prepare_rv_reference(reference_table)
stopifnot(nrow(clean) == 2L, identical(clean$Direction, c("Damaging", "Protective")))
duplicate_reference <- transform(reference_table, ID = c("RV1", "RV1"))
stopifnot(inherits(try(prepare_rv_reference(duplicate_reference), silent = TRUE),
                     "try-error"))

scz_ptv_test <- harmonise_scz_reference(data.frame(
  ID = "ENSG1", Symbol = "GENE1", CAP_freq = 0.01,
  COP_freq = 0.002, CP_OR = 5
), "PTV")
stopifnot(
  scz_ptv_test$Class == "PTV", scz_ptv_test$Case_freq == 0.01,
  scz_ptv_test$OR == 5
)
ad_lof_test <- harmonise_ad_reference(data.frame(
  Gene = "GENE2", Case_freq = 0.02,
  Control_freq = 0.005, Case_OR = 4
), "LOF")
stopifnot(
  ad_lof_test$ID == "LOF:GENE2", ad_lof_test$Class == "LOF",
  ad_lof_test$OR == 4
)
if (requireNamespace("readxl", quietly = TRUE)) {
  stopifnot(
    nrow(read_tiger_reference(
      tiger_reference_file("SCZ", "SCZ_PTV.xlsx"), "SCZ", "PTV"
    )) > 0,
    nrow(read_tiger_reference(
      tiger_reference_file("AD", "AD_LOF.xlsx"), "AD", "LOF"
    )) > 0
  )
}

apoe_reference <- apoe_genotype_reference(
  c(e2 = 0.04, e3 = 0.66, e4 = 0.30),
  c(e2 = 0.08, e3 = 0.79, e4 = 0.13)
)
stopifnot(
  nrow(apoe_reference) == 6L,
  abs(sum(apoe_reference$Case_freq) - 1) < 1e-12,
  abs(sum(apoe_reference$Control_freq) - 1) < 1e-12
)
apoe_updated <- apply_high_impact_probability(
  c(0.1, 0.1, 0.1), c("e3/e3", "e3/e4", "e4/e4"), apoe_reference
)
stopifnot(
  length(apoe_updated) == 3L,
  apoe_updated[2] > apoe_updated[1],
  apoe_updated[3] > apoe_updated[2]
)
apoe_parameters <- high_impact_prevalence_parameters(
  apoe_reference, K = 0.01, SP = 0.5
)
stopifnot(
  all(apoe_parameters$K_genotype > 0 & apoe_parameters$K_genotype < 1),
  all(apoe_parameters$SP_genotype > 0 &
        apoe_parameters$SP_genotype < 1)
)
method_updated <- high_impact_method_probability(
  c(-0.2, 0.1, 0.3), c("e3/e3", "e3/e4", "e4/e4"), apoe_reference,
  K = 0.01, SP = 0.5, method = "PAIR (summary)", r2_liability = 0.1
)
manual_updated <- vapply(seq_along(method_updated), function(i) {
  row <- apoe_parameters[
    apoe_parameters$Genotype == c("e3/e3", "e3/e4", "e4/e4")[i],
  ]
  pair_probability_summary(
    c(-0.2, 0.1, 0.3)[i], row$K_genotype,
    row$SP_genotype, 0.1
  )
}, numeric(1))
stopifnot(max(abs(method_updated - manual_updated)) < 1e-12)

# A likelihood ratio of one leaves the input probability unchanged.
neutral <- data.frame(
  Genotype = c("A", "B"),
  Case_freq = c(0.4, 0.6), Control_freq = c(0.4, 0.6)
)
stopifnot(all.equal(
  apply_high_impact_probability(c(0.2, 0.8), c("A", "B"), neutral),
  c(0.2, 0.8)
))

carrier_table <- data.frame(
  ID = c("I1", "I2", "I2", "I3"),
  Variant_ID = c("D1", "D1", "D2", "P1"),
  Carrier = 1
)
carrier_matrix <- prepare_rv_carrier_matrix(
  carrier_table,
  individual_ids = c("I1", "I2", "I3", "I4"),
  variant_ids = c("D1", "D2", "P1")
)
carrier_update <- apply_rv_carriers(
  rep(0.1, 4), carrier_matrix,
  odds_ratios = c(D1 = 4, D2 = 8, P1 = 0.4), prevalence = K
)
stopifnot(
  identical(carrier_update$RV_count, c(1, 2, 1, 0)),
  carrier_update$probability_after[2] > carrier_update$probability_after[1],
  carrier_update$probability_after[3] < 0.1,
  abs(carrier_update$probability_after[4] - 0.1) < 1e-12
)
named_probabilities <- stats::setNames(rep(0.1, 4), rownames(carrier_matrix))
stopifnot(nrow(apply_rv_carriers(
  named_probabilities, carrier_matrix,
  odds_ratios = c(D1 = 4, D2 = 8, P1 = 0.4), prevalence = K
)) == 4L)
misordered_probabilities <- named_probabilities[rev(seq_along(named_probabilities))]
stopifnot(inherits(try(apply_rv_carriers(
  misordered_probabilities, carrier_matrix,
  odds_ratios = c(D1 = 4, D2 = 8, P1 = 0.4), prevalence = K
), silent = TRUE), "try-error"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_carriers <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
  plot_before <- seq(0.1, 0.5, by = 0.1)
  plot_after <- apply_rv_probability(
    plot_before, p_damaging = ifelse(plot_carriers, 0.2, 0)
  )
  tiger_plot <- plot_tiger_before_after_rv(
    prs = seq(-1, 1, length.out = 5),
    probability_before = plot_before,
    probability_after = plot_after,
    rv_carrier = plot_carriers,
    rv_count = as.integer(plot_carriers),
    n_prs_bins = 3,
    probability_method = "PAIR (summary)"
  )
  stopifnot(inherits(tiger_plot, "ggplot"))
  tiger_plot_build <- ggplot2::ggplot_build(tiger_plot)
  tiger_point_layer <- tiger_plot_build$data[[2]]
  stopifnot(
    all(unique(tiger_point_layer$shape) %in% c(21, 24)),
    length(unique(tiger_point_layer$fill)) == 1L
  )
  tiger_apoe_plot <- plot_tiger_before_after_apoe(
    prs = seq(-1, 1, length.out = 5),
    probability_before = plot_before,
    probability_after = pmin(plot_before * c(0.8, 1, 1.2, 1.5, 2), 1),
    apoe_genotype = c("e2/e3", "e3/e3", "e3/e4", "e3/e4", "e4/e4"),
    n_prs_bins = 3
  )
  stopifnot(inherits(tiger_apoe_plot, "ggplot"))
  apoe_curve_plot <- plot_tiger_apoe_curves(
    prs = seq(-4, 4, length.out = 41),
    probability_prs = pair_probability_summary(
      seq(-4, 4, length.out = 41), K = 0.1, SP = 0.5,
      r2_liability = 0.1
    ),
    apoe_reference = apoe_reference,
    K = 0.01, SP = 0.5, r2_liability = 0.1
  )
  stopifnot(inherits(apoe_curve_plot, "ggplot"))
  carrier_point_plot <- plot_tiger_rv_carrier_points(
    prs_curve = seq(-4, 4, length.out = 41),
    probability_curve = pair_probability_summary(
      seq(-4, 4, length.out = 41), K = 0.01, SP = 0.5,
      r2_liability = 0.1
    ),
    carrier_prs = c(-1, 1),
    carrier_probability_before = c(0.1, 0.7),
    carrier_probability_after = c(0.3, 0.9),
    rv_count = c(1, 2),
    show_rv_labels = TRUE, rv_labels = c("RV_A", "RV_A; RV_B"),
    rv_point_size = 2.5, rv_point_alpha = 0.4,
    rv_shapes = c("1 RV" = 21, "2+ RVs" = 24)
  )
  stopifnot(inherits(carrier_point_plot, "ggplot"))
  carrier_point_build <- ggplot2::ggplot_build(carrier_point_plot)
  carrier_point_layer <- carrier_point_build$data[[2]]
  stopifnot(
    identical(sort(unique(carrier_point_layer$shape)), c(21, 24)),
    length(unique(carrier_point_layer$fill)) == 1L,
    identical(carrier_point_build$data[[3]]$label,
              c("RV_A", "RV_A; RV_B"))
  )
  combined_plot <- plot_tiger_apoe_rv_carrier_points(
    prs_curve = seq(-4, 4, length.out = 41),
    probability_prs = pair_probability_summary(
      seq(-4, 4, length.out = 41), K = 0.01, SP = 0.5,
      r2_liability = 0.1
    ),
    apoe_reference = apoe_reference,
    carrier_prs = c(-1, 1),
    carrier_apoe = c("e3/e4", "e4/e4"),
    carrier_probability_before_rv = c(0.2, 0.8),
    carrier_probability_after_rv = c(0.4, 0.9),
    rv_count = c(1, 2),
    carrier_group = c("Group A", "Group B"),
    K = 0.01, SP = 0.5, r2_liability = 0.1
  )
  stopifnot(inherits(combined_plot, "ggplot"))
  combined_plot_build <- ggplot2::ggplot_build(combined_plot)
  combined_point_layer <- combined_plot_build$data[[2]]
  stopifnot(
    identical(sort(unique(combined_point_layer$shape)), c(21, 24)),
    length(unique(combined_point_layer$fill)) == 2L,
    all(combined_point_layer$alpha == 1),
    identical(combined_plot$labels$fill, "Carrier group"),
    identical(combined_plot$labels$shape, "Number of RVs")
  )
}
cat("All tests passed.\n")
