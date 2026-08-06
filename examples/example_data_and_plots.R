# TIGER. GNU GPL v3 or later. See LICENSE.
# Simulation-free example using supplied synthetic input tables.
# Run from the TIGER repository root:
# Rscript examples/example_data_and_plots.R

source("R/tiger.R")
load_tiger()

individuals <- read.csv(tiger_example_file("example_individuals.csv"),
                        stringsAsFactors = FALSE)
apoe_status <- read.csv(tiger_example_file("example_apoe_status.csv"),
                        stringsAsFactors = FALSE)
rv_reference <- prepare_rv_reference(
  read.csv(tiger_example_file("example_rv_reference.csv"), stringsAsFactors = FALSE)
)
rv_carriers <- read.csv(tiger_example_file("example_rv_carriers.csv"),
                        stringsAsFactors = FALSE)
apoe_reference <- prepare_high_impact_reference(
  read.csv(tiger_example_file("example_apoe_reference.csv"), stringsAsFactors = FALSE)
)

# Join the separate APOE-status record by ID. Never assume that two files have
# the same row order.
if (anyDuplicated(individuals$ID) || anyDuplicated(apoe_status$ID)) {
  stop("individual and APOE-status IDs must be unique")
}
apoe_match <- match(individuals$ID, apoe_status$ID)
if (anyNA(apoe_match)) stop("APOE status is missing for one or more individuals")
individuals$APOE <- apoe_status$APOE[apoe_match]

# Add the presence-only RV records to the same individual table. Multiple RV
# IDs are separated by semicolons; an empty value means no carried RV.
rv_ids_by_person <- split(rv_carriers$Variant_ID, rv_carriers$ID)
individuals$RV_status <- vapply(individuals$ID, function(id) {
  carried <- rv_ids_by_person[[id]]
  if (is.null(carried)) "" else paste(unique(carried), collapse = ";")
}, character(1))

# Illustrative modelling values only; users must replace and justify these.
K <- 0.01
SP <- 0.50
r2_liability <- 0.10
SP_RV <- 0.50

# Calculate all four conditions through TIGER's main merged-table interface.
# The example PRS is assumed to have excluded the APOE region.
individuals <- tiger_probabilities(
  individuals, K = K, SP = SP, method = "PAIR (summary)",
  r2_liability = r2_liability,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_status", rv_prevalence = SP_RV,
  include_high_impact = TRUE, high_impact_col = "APOE",
  high_impact_reference = apoe_reference
)

example_output <- individuals[, c(
  "ID", "PRS_liability", "APOE", "RV_count", "Probability_PRS",
  "Probability_PRS_RV", "Probability_PRS_HIGH_IMPACT",
  "Probability_PRS_HIGH_IMPACT_RV"
)]
print(utils::head(example_output, 12))
message("Calculated probabilities for ", nrow(example_output),
        " synthetic individuals; the first 12 are shown.")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  figure_dir <- "examples/figures"
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # A dense PRS sequence is retained for the genotype-specific APOE curves.
  prs_sequence <- seq(-4, 4, length.out = 321)
  plot_K <- 0.01
  plot_SP <- 0.50
  plot_r2l <- 0.10
  plot_p_prs <- pair_probability_summary(
    prs_sequence, K = plot_K, SP = plot_SP, r2_liability = plot_r2l
  )
  reference_quantiles <- stats::qnorm(
    (seq_len(500) - 0.5) / 500
  ) * sqrt(plot_r2l)
  plot_r2o <- liability_to_observed_r2(plot_r2l, plot_K, plot_SP)
  method_plot <- plot_tiger_prs_methods(
    prs = prs_sequence,
    bpc_probability = bpc_probability(
      prs_sequence, plot_K, plot_SP, plot_r2l
    ),
    genopred_probability = genopred_probability(
      prs_sequence, reference_quantiles, plot_r2o,
      plot_K, SP = plot_SP
    ),
    pair_probability = plot_p_prs
  )
  ggplot2::ggsave(
    file.path(figure_dir, "prs_method_comparison_v1.png"), method_plot,
    width = 7.5, height = 5.2, dpi = 160
  )

  rv_reference_plot <- plot_tiger_rv_reference(
    rv_reference, prevalence = SP_RV
  )
  ggplot2::ggsave(
    file.path(figure_dir, "rv_probability_reference_v1.png"),
    rv_reference_plot, width = 7.5, height = 4.8, dpi = 160
  )
  if (anyNA(individuals$Group) || any(!nzchar(individuals$Group))) {
    stop("Example group is missing for one or more individuals")
  }
  plotted_prs_group <- individuals$Group
  carrier_index <- individuals$RV_count > 0
  plotted_carrier_group <- plotted_prs_group[carrier_index]
  rv_labels_by_id <- tapply(
    rv_carriers$Variant_ID, rv_carriers$ID,
    function(x) paste(unique(x), collapse = "; ")
  )
  plotted_rv_labels <- unname(
    rv_labels_by_id[individuals$ID[carrier_index]]
  )
  plotted_group_colours <- c(
    "Control" = "#3182BD", "Case" = "#E31A1C"
  )
  rv_plot <- plot_tiger_rv_carrier_points(
    prs_curve = individuals$PRS_liability,
    probability_curve = individuals$Probability_PRS,
    prs_group = plotted_prs_group,
    carrier_prs = individuals$PRS_liability[carrier_index],
    carrier_probability_before = individuals$Probability_PRS[carrier_index],
    carrier_probability_after = individuals$Probability_PRS_RV[carrier_index],
    rv_count = individuals$RV_count[carrier_index],
    rv_labels = plotted_rv_labels,
    carrier_group = plotted_carrier_group,
    group_colours = plotted_group_colours,
    probability_method = "PAIR (summary)"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "rv_carrier_points_v4.png"), rv_plot,
    width = 7.5, height = 5.2, dpi = 160
  )

  apoe_plot <- plot_tiger_apoe_curves(
    prs_sequence,
    plot_p_prs,
    apoe_reference,
    probability_method = "PAIR (summary)",
    K = K, SP = SP, r2_liability = r2_liability
  )
  ggplot2::ggsave(
    file.path(figure_dir, "apoe_genotype_curves_v2.png"), apoe_plot,
    width = 7.5, height = 5.2, dpi = 160
  )

  single_variant_reference <- biallelic_genotype_reference(
    case_effect_allele_frequency = 0.20,
    control_effect_allele_frequency = 0.10
  )
  single_variant_plot <- plot_tiger_high_impact_curves(
    prs_sequence,
    plot_p_prs,
    single_variant_reference,
    probability_method = "PAIR (summary)",
    K = K, SP = SP, r2_liability = r2_liability,
    genotype_colours = c("0" = "#3182BD", "1" = "#737373", "2" = "#E31A1C"),
    plot_title = "Single common high-impact variant"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "single_high_impact_variant_curves_v1.png"),
    single_variant_plot, width = 7.5, height = 5.2, dpi = 160
  )

  # Three genotype-specific sample distributions with the same labelled RV
  # overlay used in the six-genotype APOE figure.
  single_variant_individuals <- individuals
  single_variant_individuals$Variant_genotype <- rep(
    c("0", "1", "2"), each = nrow(single_variant_individuals) / 3
  )
  single_variant_individuals <- tiger_probabilities(
    single_variant_individuals,
    K = K, SP = SP, method = "PAIR (summary)",
    r2_liability = r2_liability,
    include_rv = TRUE, rv_reference = rv_reference,
    rv_status_col = "RV_status", rv_prevalence = SP_RV,
    include_high_impact = TRUE,
    high_impact_col = "Variant_genotype",
    high_impact_reference = single_variant_reference
  )
  single_variant_rv_plot <- plot_tiger_high_impact_rv_carrier_points(
    prs_curve = single_variant_individuals$PRS_liability,
    probability_prs = single_variant_individuals$Probability_PRS,
    high_impact_reference = single_variant_reference,
    carrier_prs = single_variant_individuals$PRS_liability[carrier_index],
    carrier_genotype =
      single_variant_individuals$Variant_genotype[carrier_index],
    carrier_probability_before_rv =
      single_variant_individuals$Probability_PRS_HIGH_IMPACT[carrier_index],
    carrier_probability_after_rv =
      single_variant_individuals$Probability_PRS_HIGH_IMPACT_RV[carrier_index],
    rv_count = single_variant_individuals$RV_count[carrier_index],
    prs_group = single_variant_individuals$Group,
    prs_genotype = single_variant_individuals$Variant_genotype,
    carrier_group = single_variant_individuals$Group[carrier_index],
    group_colours = plotted_group_colours,
    rv_labels = plotted_rv_labels,
    probability_method = "PAIR (summary)",
    K = K, SP = SP, r2_liability = r2_liability,
    genotype_colours = c("0" = "#3182BD", "1" = "#737373", "2" = "#E31A1C")
  )
  ggplot2::ggsave(
    file.path(
      figure_dir, "single_high_impact_variant_rv_carrier_points_v1.png"
    ),
    single_variant_rv_plot, width = 11.5, height = 7.2, dpi = 160
  )

  apoe_rv_plot <- plot_tiger_apoe_rv_carrier_points(
    prs_curve = individuals$PRS_liability,
    probability_prs = individuals$Probability_PRS,
    apoe_reference = apoe_reference,
    carrier_prs = individuals$PRS_liability[carrier_index],
    carrier_apoe = individuals$APOE[carrier_index],
    carrier_probability_before_rv =
      individuals$Probability_PRS_HIGH_IMPACT[carrier_index],
    carrier_probability_after_rv =
      individuals$Probability_PRS_HIGH_IMPACT_RV[carrier_index],
    rv_count = individuals$RV_count[carrier_index],
    prs_group = plotted_prs_group,
    prs_apoe = individuals$APOE,
    rv_labels = plotted_rv_labels,
    carrier_group = plotted_carrier_group,
    group_colours = plotted_group_colours,
    probability_method = "PAIR (summary)",
    K = K, SP = SP, r2_liability = r2_liability
  )
  ggplot2::ggsave(
    file.path(figure_dir, "apoe_rv_carrier_points_v1.png"), apoe_rv_plot,
    width = 11.5, height = 7.2, dpi = 160
  )
}
