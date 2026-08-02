# TIGER: Translatable Integrated Genetic Risk framework.
# Copyright (C) 2026 TIGER study authors
# Licensed under GNU GPL v3 or later; distributed WITHOUT ANY WARRANTY.

.check_plot_probability_vector <- function(x, name, n = NULL) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(x < 0) || any(x > 1)) {
    stop(name, " must contain finite probabilities in [0, 1]")
  }
  if (!is.null(n) && length(x) != n) stop(name, " must have length ", n)
}

# Shared TIGER visual theme. All exported plotting helpers use this specification
# so figures have consistent typography, axes, panel treatment and legends.
tiger_plot_theme <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("tiger_plot_theme() requires the ggplot2 package")
  }
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "grey25", linewidth = 0.6),
      axis.title = ggplot2::element_text(colour = "grey10"),
      axis.text = ggplot2::element_text(colour = "grey25"),
      plot.title = ggplot2::element_text(face = "bold", colour = "grey10"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

# Plot PRS-based probabilities before and after an RV update.
#
# The PRS itself does not change when an RV is added. Lines show equal-count
# PRS-bin means; optional points show individual RV-adjusted probabilities.
plot_tiger_before_after_rv <- function(
    prs, probability_before, probability_after,
    rv_carrier = NULL, rv_count = NULL, show_rv_points = TRUE,
    point_selection = c("carriers", "all"), n_prs_bins = 50,
    probability_method = "PAIR (summary)", y_limits = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_before_after_rv() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || !length(prs) || any(!is.finite(prs))) {
    stop("prs must be a non-empty numeric vector of finite values")
  }
  n <- length(prs)
  .check_plot_probability_vector(probability_before, "probability_before", n)
  .check_plot_probability_vector(probability_after, "probability_after", n)
  point_selection <- match.arg(point_selection)
  if (!is.numeric(n_prs_bins) || length(n_prs_bins) != 1L ||
      !is.finite(n_prs_bins) || n_prs_bins < 2) {
    stop("n_prs_bins must be one finite value of at least 2")
  }
  n_prs_bins <- min(as.integer(n_prs_bins), n)
  if (n_prs_bins < 2L) stop("at least two observations are required")
  if (is.null(y_limits)) {
    upper <- max(c(probability_before, probability_after)) * 1.08
    y_limits <- c(0, min(1, max(upper, 0.05)))
  } else if (!is.numeric(y_limits) || length(y_limits) != 2L ||
             any(!is.finite(y_limits)) || y_limits[1] < 0 ||
             y_limits[2] > 1 || y_limits[1] >= y_limits[2]) {
    stop("y_limits must be NULL or increasing limits within [0, 1]")
  }

  if (is.null(rv_carrier)) {
    # A changed probability is a useful default when an explicit carrier flag
    # is unavailable. Supplying rv_carrier is preferred for transparent plots.
    rv_carrier <- abs(probability_after - probability_before) >
      sqrt(.Machine$double.eps)
  }
  if (!is.logical(rv_carrier) || length(rv_carrier) != n || anyNA(rv_carrier)) {
    stop("rv_carrier must be a logical vector with one value per person")
  }
  if (is.null(rv_count)) rv_count <- as.integer(rv_carrier)
  if (!is.numeric(rv_count) || length(rv_count) != n || any(!is.finite(rv_count)) ||
      any(rv_count < 0) || any(rv_count != as.integer(rv_count)) ||
      any((rv_count > 0) != rv_carrier)) {
    stop("rv_count must be non-negative integers consistent with rv_carrier")
  }

  rank_bin <- pmin(
    ceiling(rank(prs, ties.method = "first") / n * n_prs_bins),
    n_prs_bins
  )
  long <- rbind(
    data.frame(PRS = prs, Probability = probability_before,
               Condition = "Before RV", Bin = rank_bin),
    data.frame(PRS = prs, Probability = probability_after,
               Condition = "After RV", Bin = rank_bin)
  )
  long$Condition <- factor(long$Condition,
                           levels = c("Before RV", "After RV"))
  curves <- stats::aggregate(
    cbind(PRS, Probability) ~ Condition + Bin,
    data = long, FUN = mean
  )
  curves <- curves[order(curves$Condition, curves$Bin), , drop = FALSE]

  point_index <- if (point_selection == "carriers") rv_carrier else rep(TRUE, n)
  points <- data.frame(
    PRS = prs[point_index], Probability = probability_after[point_index],
    RV_group = ifelse(rv_count[point_index] >= 2, "2+ RVs",
                      ifelse(rv_count[point_index] == 1, "1 RV", "No RV"))
  )
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) {
    ""
  } else {
    paste0(": ", trimws(probability_method))
  }

  plot <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = PRS, y = Probability, colour = Condition,
                 linetype = Condition, group = Condition)
  ) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = c("Before RV" = "#3182BD", "After RV" = "#E31A1C")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Before RV" = "dashed", "After RV" = "solid")
    ) +
    ggplot2::scale_y_continuous() +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      x = "Liability PRS",
      y = "Estimated disorder probability",
      colour = NULL, linetype = NULL,
      title = paste0("PRS-based probability before and after RV", method_text),
      subtitle = paste0(
        "Lines are means across ", n_prs_bins,
        " equal-count PRS bins; the PRS itself is unchanged"
      )
    ) +
    tiger_plot_theme()

  if (isTRUE(show_rv_points) && nrow(points)) {
    plot <- plot + ggplot2::geom_point(
      data = points,
      ggplot2::aes(x = PRS, y = Probability, fill = RV_group),
      inherit.aes = FALSE, shape = 21, colour = "grey25",
      alpha = 0.55, size = 1.8, stroke = 0.25
    ) +
      ggplot2::scale_fill_manual(
        values = c("No RV" = "grey75", "1 RV" = "#FDB863", "2+ RVs" = "#B2182B"),
        drop = FALSE
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Carrier status"))
  }
  plot
}

# Plot PRS-based probability before and after a separate APOE update. Lines are
# equal-count PRS-bin means and optional points show genotype-specific adjusted
# probabilities for individuals.
plot_tiger_before_after_apoe <- function(
    prs, probability_before, probability_after, apoe_genotype,
    show_apoe_points = TRUE, n_prs_bins = 50,
    probability_method = "PAIR (summary)", y_limits = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_before_after_apoe() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || length(prs) < 2L || any(!is.finite(prs))) {
    stop("prs must contain at least two finite numeric values")
  }
  n <- length(prs)
  .check_plot_probability_vector(probability_before, "probability_before", n)
  .check_plot_probability_vector(probability_after, "probability_after", n)
  apoe_genotype <- as.character(apoe_genotype)
  if (length(apoe_genotype) != n || anyNA(apoe_genotype) ||
      any(!nzchar(apoe_genotype))) stop("apoe_genotype must have one label per person")
  n_prs_bins <- min(as.integer(n_prs_bins), n)
  if (!is.finite(n_prs_bins) || n_prs_bins < 2L) {
    stop("n_prs_bins must be at least 2")
  }
  if (is.null(y_limits)) {
    upper <- max(c(probability_before, probability_after)) * 1.08
    y_limits <- c(0, min(1, max(upper, 0.05)))
  } else if (!is.numeric(y_limits) || length(y_limits) != 2L ||
             any(!is.finite(y_limits)) || y_limits[1] < 0 ||
             y_limits[2] > 1 || y_limits[1] >= y_limits[2]) {
    stop("y_limits must be NULL or increasing limits within [0, 1]")
  }
  rank_bin <- pmin(ceiling(rank(prs, ties.method = "first") / n * n_prs_bins),
                   n_prs_bins)
  long <- rbind(
    data.frame(PRS = prs, Probability = probability_before,
               Condition = "PRS", Bin = rank_bin),
    data.frame(PRS = prs, Probability = probability_after,
               Condition = "PRS + APOE", Bin = rank_bin)
  )
  long$Condition <- factor(long$Condition, levels = c("PRS", "PRS + APOE"))
  curves <- stats::aggregate(cbind(PRS, Probability) ~ Condition + Bin,
                             long, mean)
  curves <- curves[order(curves$Condition, curves$Bin), , drop = FALSE]
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  plot <- ggplot2::ggplot(
    curves,
    ggplot2::aes(PRS, Probability, colour = Condition,
                 linetype = Condition, group = Condition)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = c("PRS" = "#3182BD",
                                            "PRS + APOE" = "#E66101")) +
    ggplot2::scale_linetype_manual(values = c("PRS" = "dashed",
                                              "PRS + APOE" = "solid")) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = NULL, linetype = NULL,
      title = paste0("PRS-based probability before and after APOE", method_text),
      subtitle = "APOE is modelled separately; the PRS itself is unchanged"
    ) +
    tiger_plot_theme()
  if (isTRUE(show_apoe_points)) {
    points <- data.frame(PRS = prs, Probability = probability_after,
                         APOE = apoe_genotype)
    plot <- plot + ggplot2::geom_point(
      data = points, ggplot2::aes(PRS, Probability, fill = APOE),
      inherit.aes = FALSE, shape = 21, colour = "grey25",
      alpha = 0.65, size = 1.8, stroke = 0.25
    ) + ggplot2::guides(
      fill = ggplot2::guide_legend(title = "APOE genotype", nrow = 1, order = 1),
      colour = ggplot2::guide_legend(order = 2),
      linetype = ggplot2::guide_legend(order = 2)
    )
  }
  plot
}

# Instructional APOE genotype curves across a supplied PRS sequence, following
# the visual structure used in APOE6_Freq.R. The case/control reference is used
# to update the same PRS probability curve for all six genotypes.
plot_tiger_apoe_curves <- function(
    prs, probability_prs, apoe_reference,
    probability_method = "PAIR (summary)", y_limits = c(0, 1)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_apoe_curves() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || length(prs) < 2L || any(!is.finite(prs))) {
    stop("prs must contain at least two finite numeric values")
  }
  .check_plot_probability_vector(probability_prs, "probability_prs", length(prs))
  reference <- prepare_high_impact_reference(apoe_reference)
  preferred <- c("e4/e4", "e3/e4", "e2/e4", "e3/e3", "e2/e3", "e2/e2")
  genotype_order <- c(preferred[preferred %in% reference$Genotype],
                      setdiff(reference$Genotype, preferred))
  rows <- lapply(genotype_order, function(genotype) {
    data.frame(
      PRS = prs,
      Probability = apply_high_impact_probability(
        probability_prs, genotype, reference
      ),
      Genotype = genotype
    )
  })
  long <- do.call(rbind, rows)
  long$Genotype <- factor(long$Genotype, levels = genotype_order)
  genotype_colours <- c(
    "e4/e4" = "#7570B3", "e3/e4" = "#E6AB02", "e2/e4" = "#D95F02",
    "e3/e3" = "grey50", "e2/e3" = "#1B9E77", "e2/e2" = "#E7298A"
  )
  missing_colours <- setdiff(genotype_order, names(genotype_colours))
  if (length(missing_colours)) {
    genotype_colours <- c(genotype_colours,
      stats::setNames(grDevices::hcl.colors(length(missing_colours), "Dark 3"),
                      missing_colours))
  }
  genotype_lines <- stats::setNames(rep("solid", length(genotype_order)),
                                    genotype_order)
  if ("e3/e3" %in% genotype_order) genotype_lines["e3/e3"] <- "dotted"
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  ggplot2::ggplot(
    long,
    ggplot2::aes(PRS, Probability, colour = Genotype, linetype = Genotype)
  ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                       colour = "grey75", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_linetype_manual(values = genotype_lines) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::coord_cartesian(xlim = range(prs), ylim = y_limits) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = "APOE genotype", linetype = "APOE genotype",
      title = paste0("APOE genotype probability curves", method_text),
      subtitle = "APOE is modelled separately from the PRS"
    ) +
    tiger_plot_theme()
}

# Instructional/observed RV-carrier overlay. The line is the unadjusted PRS
# probability curve. Only people carrying at least one RV are plotted, at their
# RV-adjusted probabilities. Shape distinguishes one from multiple carried RVs.
plot_tiger_rv_carrier_points <- function(
    prs_curve, probability_curve,
    carrier_prs, carrier_probability_before, carrier_probability_after,
    rv_count, rv_effect = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    show_shift_segments = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_rv_carrier_points() requires the ggplot2 package")
  }
  if (!is.numeric(prs_curve) || length(prs_curve) < 2L ||
      any(!is.finite(prs_curve))) {
    stop("prs_curve must contain at least two finite values")
  }
  .check_plot_probability_vector(
    probability_curve, "probability_curve", length(prs_curve)
  )
  n <- length(carrier_prs)
  if (!is.numeric(carrier_prs) || !n || any(!is.finite(carrier_prs))) {
    stop("carrier_prs must contain finite PRSs for RV carriers")
  }
  .check_plot_probability_vector(
    carrier_probability_before, "carrier_probability_before", n
  )
  .check_plot_probability_vector(
    carrier_probability_after, "carrier_probability_after", n
  )
  if (!is.numeric(rv_count) || length(rv_count) != n ||
      any(!is.finite(rv_count)) || any(rv_count < 1) ||
      any(rv_count != as.integer(rv_count))) {
    stop("rv_count must contain positive integer counts for RV carriers")
  }
  if (is.null(rv_effect)) rv_effect <- rep("RV-adjusted carrier", n)
  rv_effect <- as.character(rv_effect)
  if (length(rv_effect) != n || anyNA(rv_effect) || any(!nzchar(rv_effect))) {
    stop("rv_effect must contain one non-empty label per carrier")
  }
  carriers <- data.frame(
    PRS = carrier_prs,
    Probability_before = carrier_probability_before,
    Probability_after = carrier_probability_after,
    RV_count_group = factor(ifelse(rv_count == 1, "1 RV", "2+ RVs"),
                            levels = c("1 RV", "2+ RVs")),
    RV_effect = rv_effect
  )
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  plot <- ggplot2::ggplot(
    data.frame(PRS = prs_curve, Probability = probability_curve),
    ggplot2::aes(PRS, Probability)
  ) +
    ggplot2::geom_line(colour = "#3182BD", linewidth = 1.1) +
    ggplot2::coord_cartesian(xlim = range(prs_curve), ylim = y_limits) +
    ggplot2::scale_shape_manual(values = c("1 RV" = 24, "2+ RVs" = 22)) +
    ggplot2::scale_fill_manual(values = c(
      "Damaging" = "#F8766D", "Mixed" = "#00BA38",
      "Protective" = "#619CFF", "RV-adjusted carrier" = "#FDB863"
    )) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1, override.aes = list(shape = 21, size = 3)
      ),
      shape = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      shape = "Number carried", fill = "RV effect",
      title = paste0("RV carrier probabilities", method_text),
      subtitle = "The curve is PRS only; points are shown only for RV carriers"
    ) +
    tiger_plot_theme()
  if (isTRUE(show_shift_segments)) {
    plot <- plot + ggplot2::geom_segment(
      data = carriers,
      ggplot2::aes(x = PRS, xend = PRS,
                   y = Probability_before, yend = Probability_after),
      inherit.aes = FALSE, colour = "grey55", linewidth = 0.45,
      linetype = "dashed"
    )
  }
  plot + ggplot2::geom_point(
    data = carriers,
    ggplot2::aes(PRS, Probability_after,
                 shape = RV_count_group, fill = RV_effect),
    inherit.aes = FALSE, colour = "grey20", size = 3, stroke = 0.55
  )
}

# Combined APOE + RV figure. Lines show PRS + APOE probability curves for all
# high-impact genotypes. Only RV carriers are overlaid, at their final
# PRS + APOE + RV probabilities.
plot_tiger_apoe_rv_carrier_points <- function(
    prs_curve, probability_prs, apoe_reference,
    carrier_prs, carrier_apoe,
    carrier_probability_before_rv, carrier_probability_after_rv,
    rv_count, rv_effect = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_apoe_rv_carrier_points() requires the ggplot2 package")
  }
  if (!is.numeric(prs_curve) || length(prs_curve) < 2L ||
      any(!is.finite(prs_curve))) {
    stop("prs_curve must contain at least two finite values")
  }
  .check_plot_probability_vector(
    probability_prs, "probability_prs", length(prs_curve)
  )
  reference <- prepare_high_impact_reference(apoe_reference)
  preferred <- c("e4/e4", "e3/e4", "e2/e4", "e3/e3", "e2/e3", "e2/e2")
  genotype_order <- c(preferred[preferred %in% reference$Genotype],
                      setdiff(reference$Genotype, preferred))
  curve_rows <- lapply(genotype_order, function(genotype) {
    data.frame(
      PRS = prs_curve,
      Probability = apply_high_impact_probability(
        probability_prs, genotype, reference
      ),
      Genotype = genotype
    )
  })
  curves <- do.call(rbind, curve_rows)
  curves$Genotype <- factor(curves$Genotype, levels = genotype_order)

  n <- length(carrier_prs)
  if (!is.numeric(carrier_prs) || !n || any(!is.finite(carrier_prs))) {
    stop("carrier_prs must contain finite PRSs for RV carriers")
  }
  carrier_apoe <- as.character(carrier_apoe)
  if (length(carrier_apoe) != n || anyNA(carrier_apoe) ||
      any(!carrier_apoe %in% reference$Genotype)) {
    stop("carrier_apoe must contain one known genotype per RV carrier")
  }
  .check_plot_probability_vector(
    carrier_probability_before_rv, "carrier_probability_before_rv", n
  )
  .check_plot_probability_vector(
    carrier_probability_after_rv, "carrier_probability_after_rv", n
  )
  if (!is.numeric(rv_count) || length(rv_count) != n ||
      any(!is.finite(rv_count)) || any(rv_count < 1) ||
      any(rv_count != as.integer(rv_count))) {
    stop("rv_count must contain positive integer counts for RV carriers")
  }
  if (is.null(rv_effect)) rv_effect <- rep("RV-adjusted carrier", n)
  rv_effect <- as.character(rv_effect)
  if (length(rv_effect) != n || anyNA(rv_effect) || any(!nzchar(rv_effect))) {
    stop("rv_effect must contain one non-empty label per carrier")
  }
  carriers <- data.frame(
    PRS = carrier_prs,
    APOE = carrier_apoe,
    Probability_before_RV = carrier_probability_before_rv,
    Probability_after_RV = carrier_probability_after_rv,
    RV_count_group = factor(ifelse(rv_count == 1, "1 RV", "2+ RVs"),
                            levels = c("1 RV", "2+ RVs")),
    RV_effect = rv_effect
  )
  genotype_colours <- c(
    "e4/e4" = "#7570B3", "e3/e4" = "#E6AB02", "e2/e4" = "#D95F02",
    "e3/e3" = "grey50", "e2/e3" = "#1B9E77", "e2/e2" = "#E7298A"
  )
  missing_colours <- setdiff(genotype_order, names(genotype_colours))
  if (length(missing_colours)) {
    genotype_colours <- c(genotype_colours,
      stats::setNames(grDevices::hcl.colors(length(missing_colours), "Dark 3"),
                      missing_colours))
  }
  genotype_lines <- stats::setNames(rep("solid", length(genotype_order)),
                                    genotype_order)
  if ("e3/e3" %in% genotype_order) genotype_lines["e3/e3"] <- "dotted"
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  ggplot2::ggplot(
    curves,
    ggplot2::aes(PRS, Probability, colour = Genotype, linetype = Genotype)
  ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(
      data = carriers,
      ggplot2::aes(PRS, Probability_after_RV,
                   shape = RV_count_group, fill = RV_effect),
      inherit.aes = FALSE, colour = "grey15", size = 3.1, stroke = 0.6
    ) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_linetype_manual(values = genotype_lines) +
    ggplot2::scale_shape_manual(values = c("1 RV" = 24, "2+ RVs" = 22)) +
    ggplot2::scale_fill_manual(values = c(
      "Damaging" = "#F8766D", "Mixed" = "#00BA38",
      "Protective" = "#619CFF", "RV-adjusted carrier" = "#FDB863"
    )) +
    ggplot2::coord_cartesian(xlim = range(prs_curve), ylim = y_limits) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1),
      fill = ggplot2::guide_legend(
        order = 2, override.aes = list(shape = 21, size = 3)
      ),
      shape = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = "APOE genotype", linetype = "APOE genotype",
      fill = "RV effect", shape = "Number carried",
      title = paste0("APOE curves with RV carriers", method_text),
      subtitle = paste0(
        "Curves: PRS + APOE; points: PRS + APOE + RV (carriers only)"
      )
    ) +
    tiger_plot_theme()
}
