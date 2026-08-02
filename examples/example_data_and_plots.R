# TIGER framework. GNU GPL v3 or later. See LICENSE.
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
carrier_groups <- read.csv(tiger_example_file("example_carrier_groups.csv"),
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

# Illustrative modelling values only; users must replace and justify these.
K <- 0.01
SP <- 0.50
r2_liability <- 0.10
SP_RV <- 0.50

# 1. PRS probability. PAIR (summary) is the worked-example default.
individuals$Probability_PRS <- pair_probability_summary(
  individuals$PRS_liability, K = K, SP = SP,
  r2_liability = r2_liability
)

# 2. APOE update. The example PRS is assumed to have excluded the APOE region.
individuals$Probability_PRS_APOE <- high_impact_method_probability(
  individuals$PRS_liability, individuals$APOE, apoe_reference,
  K = K, SP = SP, method = "PAIR (summary)",
  r2_liability = r2_liability
)

# 3. RV update from a separate presence-only carrier table. Because the full ID
# and variant sets are supplied, omitted pairs are interpreted as non-carriers.
carrier_matrix <- prepare_rv_carrier_matrix(
  rv_carriers,
  individual_ids = individuals$ID,
  variant_ids = rv_reference$ID
)
odds_ratios <- stats::setNames(rv_reference$OR, rv_reference$ID)
rv_update <- apply_rv_carriers(
  individuals$Probability_PRS,
  carrier_matrix = carrier_matrix,
  odds_ratios = odds_ratios,
  prevalence = SP_RV
)
individuals$RV_count <- rv_update$RV_count
individuals$Probability_PRS_RV <- rv_update$probability_after

# 4. Apply the same RV layer after APOE.
combined_update <- apply_rv_carriers(
  individuals$Probability_PRS_APOE,
  carrier_matrix = carrier_matrix,
  odds_ratios = odds_ratios,
  prevalence = SP_RV
)
individuals$Probability_PRS_APOE_RV <- combined_update$probability_after

print(individuals[, c(
  "ID", "PRS_liability", "APOE", "RV_count", "Probability_PRS",
  "Probability_PRS_RV", "Probability_PRS_APOE",
  "Probability_PRS_APOE_RV"
)])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  figure_dir <- "examples/figures"
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Instructional distributions are separate from the individual carrier file.
  # A dense standard liability-PRS sequence produces smooth, interpretable
  # curves; it is not a simulated cohort.
  prs_sequence <- seq(-4, 4, by = 0.05)
  plot_K <- 0.01
  plot_SP <- 0.50
  plot_r2l <- 0.10
  plot_p_prs <- pair_probability_summary(
    prs_sequence, K = plot_K, SP = plot_SP, r2_liability = plot_r2l
  )
  carrier_index <- rv_update$RV_count > 0
  carrier_group_match <- match(individuals$ID[carrier_index], carrier_groups$ID)
  if (anyNA(carrier_group_match)) {
    stop("Example carrier group is missing for one or more RV carriers")
  }
  plotted_carrier_group <-
    carrier_groups$Carrier_group[carrier_group_match]
  rv_labels_by_id <- tapply(
    rv_carriers$Variant_ID, rv_carriers$ID,
    function(x) paste(unique(x), collapse = "; ")
  )
  plotted_rv_labels <- unname(
    rv_labels_by_id[individuals$ID[carrier_index]]
  )
  plotted_group_colours <- c(
    "Group A" = "#3182BD", "Group B" = "#E31A1C"
  )
  rv_plot <- plot_tiger_rv_carrier_points(
    prs_curve = prs_sequence,
    probability_curve = plot_p_prs,
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

  apoe_rv_plot <- plot_tiger_apoe_rv_carrier_points(
    prs_curve = prs_sequence,
    probability_prs = plot_p_prs,
    apoe_reference = apoe_reference,
    carrier_prs = individuals$PRS_liability[carrier_index],
    carrier_apoe = individuals$APOE[carrier_index],
    carrier_probability_before_rv =
      individuals$Probability_PRS_APOE[carrier_index],
    carrier_probability_after_rv =
      individuals$Probability_PRS_APOE_RV[carrier_index],
    rv_count = individuals$RV_count[carrier_index],
    rv_labels = plotted_rv_labels,
    carrier_group = plotted_carrier_group,
    group_colours = plotted_group_colours,
    probability_method = "PAIR (summary)",
    K = K, SP = SP, r2_liability = r2_liability
  )
  ggplot2::ggsave(
    file.path(figure_dir, "apoe_rv_carrier_points_v1.png"), apoe_rv_plot,
    width = 10.5, height = 6.5, dpi = 160
  )
}
