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

# Shared RV-carrier point specification. Fill is reserved for optional external
# groups. Shape alone identifies whether one or multiple RVs are carried.
.tiger_rv_shapes <- c("No RV" = 22, "1 RV" = 21, "2+ RVs" = 24)
.tiger_rv_point_alpha <- 1
.tiger_rv_point_size <- 1.25
.tiger_rv_shape_legend <- list(
  alpha = 0.90, size = 3.2, fill = "grey65",
  colour = "grey25", stroke = 0.45
)
.tiger_group_legend <- list(
  alpha = 1.00, size = 3.2, shape = 21,
  colour = "grey25", stroke = 0.45
)

# Add optional RV-name labels. ggrepel provides guide lines when labels would
# otherwise collide. Keeping it optional avoids making plotting dependencies
# mandatory for users who do not need labels.
.add_tiger_rv_labels <- function(plot, data, y_column) {
  mapping <- ggplot2::aes(x = PRS)
  mapping$y <- as.name(y_column)
  mapping$label <- as.name("RV_label")
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot + ggrepel::geom_text_repel(
      data = data, mapping = mapping, inherit.aes = FALSE,
      size = 3, colour = "grey15", box.padding = 0.25,
      point.padding = 0.15, min.segment.length = 0,
      segment.colour = "grey55", segment.size = 0.3,
      max.overlaps = Inf
    )
  } else {
    plot + ggplot2::geom_text(
      data = data, mapping = mapping, inherit.aes = FALSE,
      nudge_y = 0.02, size = 3, colour = "grey15",
      check_overlap = TRUE
    )
  }
}

.check_tiger_point_style <- function(size, alpha, border, stroke, shapes,
                                     required_shapes) {
  if (length(size) != 1L || !is.finite(size) || size <= 0) {
    stop("rv_point_size must be one positive finite value")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("rv_point_alpha must be one finite value in [0, 1]")
  }
  if (length(border) != 1L || is.na(border) || !nzchar(border)) {
    stop("rv_point_border must be one colour")
  }
  if (length(stroke) != 1L || !is.finite(stroke) || stroke < 0) {
    stop("rv_point_stroke must be one non-negative finite value")
  }
  if (is.null(names(shapes)) || !all(required_shapes %in% names(shapes)) ||
      any(!is.finite(shapes[required_shapes]))) {
    stop("rv_shapes must be a named vector containing: ",
         paste(required_shapes, collapse = ", "))
  }
  invisible(TRUE)
}

# Prepare an optional user-defined grouping for carrier-point colours. RV effect
# direction is deliberately not inferred here. If no group is supplied, all
# carriers receive one neutral fill and no colour legend is drawn.
.prepare_tiger_carrier_group <- function(carrier_group, n,
                                         group_colours = NULL) {
  supplied <- !is.null(carrier_group)
  if (!supplied) carrier_group <- rep("RV carrier", n)
  carrier_group <- as.character(carrier_group)
  if (length(carrier_group) != n || anyNA(carrier_group) ||
      any(!nzchar(trimws(carrier_group)))) {
    stop("carrier_group must contain one non-empty group label per person")
  }
  observed_levels <- unique(carrier_group)
  if (is.null(group_colours)) {
    levels <- observed_levels
    colours <- if (length(levels) == 1L) {
      stats::setNames("#7F7F7F", levels)
    } else {
      stats::setNames(grDevices::hcl.colors(length(levels), "Dark 3"), levels)
    }
  } else {
    if (is.null(names(group_colours)) ||
        !all(observed_levels %in% names(group_colours))) {
      stop("group_colours must be a named vector covering every carrier group")
    }
    levels <- names(group_colours)[names(group_colours) %in% observed_levels]
    colours <- group_colours[levels]
  }
  list(
    group = factor(carrier_group, levels = levels),
    colours = colours,
    show_legend = supplied
  )
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
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = grid::unit(5, "pt"),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2),
      legend.title = ggplot2::element_text(face = "bold", size = 9.5),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

# Plot PRS-based probabilities before and after an RV update.
#
# The PRS itself does not change when an RV is added. Lines show equal-count
# PRS-bin means; optional points show individual RV-adjusted probabilities.
plot_tiger_before_after_rv <- function(
    prs, probability_before, probability_after,
    rv_carrier = NULL, rv_count = NULL, carrier_group = NULL,
    group_colours = NULL, show_rv_points = TRUE,
    show_rv_labels = TRUE, rv_labels = NULL,
    point_selection = c("carriers", "all"), n_prs_bins = 50,
    probability_method = "PAIR (summary)", y_limits = NULL,
    rv_point_size = .tiger_rv_point_size,
    rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey25", rv_point_stroke = 0.18,
    rv_shapes = .tiger_rv_shapes) {
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
  if (isTRUE(show_rv_labels) && !is.null(rv_labels) &&
      (length(rv_labels) != n || anyNA(rv_labels))) {
    stop("rv_labels must contain one non-missing label per person")
  }
  groups <- .prepare_tiger_carrier_group(carrier_group, n, group_colours)
  .check_tiger_point_style(
    rv_point_size, rv_point_alpha, rv_point_border, rv_point_stroke,
    rv_shapes, c("No RV", "1 RV", "2+ RVs")
  )

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
                      ifelse(rv_count[point_index] == 1, "1 RV", "No RV")),
    Carrier_group = groups$group[point_index],
    RV_label = if (is.null(rv_labels)) "" else as.character(rv_labels[point_index])
  )
  points$RV_group <- factor(
    points$RV_group, levels = c("No RV", "1 RV", "2+ RVs")
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
      ggplot2::aes(
        x = PRS, y = Probability, shape = RV_group, fill = Carrier_group
      ),
      inherit.aes = FALSE, colour = rv_point_border,
      alpha = rv_point_alpha,
      size = rv_point_size, stroke = rv_point_stroke
    ) +
      ggplot2::scale_shape_manual(
        values = rv_shapes, drop = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = groups$colours,
        guide = if (groups$show_legend) "legend" else "none"
      ) +
      ggplot2::labs(
        shape = "Number of RVs",
        fill = if (groups$show_legend) "Carrier group" else NULL
      ) +
      ggplot2::guides(
        shape = ggplot2::guide_legend(
          order = 2, nrow = 1, byrow = TRUE,
          override.aes = .tiger_rv_shape_legend
        ),
        fill = ggplot2::guide_legend(
          order = 3, nrow = 1, byrow = TRUE,
          override.aes = .tiger_group_legend
        )
      )
    if (isTRUE(show_rv_labels) && !is.null(rv_labels)) {
      plot <- .add_tiger_rv_labels(plot, points, "Probability")
    }
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
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    K = NULL, SP = 0.5, r2_liability = NULL,
    apoe_update = c("method-specific", "direct")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_apoe_curves() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || length(prs) < 2L || any(!is.finite(prs))) {
    stop("prs must contain at least two finite numeric values")
  }
  .check_plot_probability_vector(probability_prs, "probability_prs", length(prs))
  apoe_update <- match.arg(apoe_update)
  if (apoe_update == "method-specific" && is.null(K)) {
    stop("K is required for TIGER's method-specific APOE update")
  }
  reference <- prepare_high_impact_reference(apoe_reference)
  preferred <- c("e4/e4", "e3/e4", "e2/e4", "e3/e3", "e2/e3", "e2/e2")
  genotype_order <- c(preferred[preferred %in% reference$Genotype],
                      setdiff(reference$Genotype, preferred))
  rows <- lapply(genotype_order, function(genotype) {
    probability <- if (apoe_update == "direct") {
      apply_high_impact_probability(probability_prs, genotype, reference)
    } else {
      high_impact_method_probability(
        prs, genotype, reference, K = K, SP = SP,
        method = probability_method, r2_liability = r2_liability
      )
    }
    data.frame(
      PRS = prs,
      Probability = probability,
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
    rv_count, carrier_group = NULL, group_colours = NULL,
    show_rv_labels = TRUE, rv_labels = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    show_shift_segments = FALSE, rv_point_size = 2.1,
    rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey25", rv_point_stroke = 0.35,
    rv_shapes = .tiger_rv_shapes[c("1 RV", "2+ RVs")]) {
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
  if (isTRUE(show_rv_labels) && !is.null(rv_labels) &&
      (length(rv_labels) != n || anyNA(rv_labels))) {
    stop("rv_labels must contain one non-missing label per carrier")
  }
  groups <- .prepare_tiger_carrier_group(carrier_group, n, group_colours)
  .check_tiger_point_style(
    rv_point_size, rv_point_alpha, rv_point_border, rv_point_stroke,
    rv_shapes, c("1 RV", "2+ RVs")
  )
  carriers <- data.frame(
    PRS = carrier_prs,
    Probability_before = carrier_probability_before,
    Probability_after = carrier_probability_after,
    RV_count_group = factor(ifelse(rv_count == 1, "1 RV", "2+ RVs"),
                            levels = c("1 RV", "2+ RVs")),
    Carrier_group = groups$group,
    RV_label = if (is.null(rv_labels)) "" else as.character(rv_labels)
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
    ggplot2::scale_shape_manual(values = rv_shapes) +
    ggplot2::scale_fill_manual(
      values = groups$colours,
      guide = if (groups$show_legend) "legend" else "none"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        override.aes = .tiger_rv_shape_legend
      ),
      fill = ggplot2::guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        override.aes = .tiger_group_legend
      )
    ) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      shape = "Number of RVs",
      fill = if (groups$show_legend) "Carrier group" else NULL,
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
  plot <- plot + ggplot2::geom_point(
    data = carriers,
    ggplot2::aes(PRS, Probability_after,
                 shape = RV_count_group, fill = Carrier_group),
    inherit.aes = FALSE, colour = rv_point_border,
    alpha = rv_point_alpha,
    size = rv_point_size, stroke = rv_point_stroke
  )
  if (isTRUE(show_rv_labels) && !is.null(rv_labels)) {
    plot <- .add_tiger_rv_labels(plot, carriers, "Probability_after")
  }
  plot
}

# Combined APOE + RV figure. Lines show PRS + APOE probability curves for all
# high-impact genotypes. Only RV carriers are overlaid, at their final
# PRS + APOE + RV probabilities.
plot_tiger_apoe_rv_carrier_points <- function(
    prs_curve, probability_prs, apoe_reference,
    carrier_prs, carrier_apoe,
    carrier_probability_before_rv, carrier_probability_after_rv,
    rv_count, carrier_group = NULL, group_colours = NULL,
    show_rv_labels = TRUE, rv_labels = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    K = NULL, SP = 0.5, r2_liability = NULL,
    apoe_update = c("method-specific", "direct"),
    rv_point_size = 2.1, rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey25", rv_point_stroke = 0.35,
    rv_shapes = .tiger_rv_shapes[c("1 RV", "2+ RVs")]) {
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
  apoe_update <- match.arg(apoe_update)
  if (apoe_update == "method-specific" && is.null(K)) {
    stop("K is required for TIGER's method-specific APOE update")
  }
  reference <- prepare_high_impact_reference(apoe_reference)
  preferred <- c("e4/e4", "e3/e4", "e2/e4", "e3/e3", "e2/e3", "e2/e2")
  genotype_order <- c(preferred[preferred %in% reference$Genotype],
                      setdiff(reference$Genotype, preferred))
  curve_rows <- lapply(genotype_order, function(genotype) {
    probability <- if (apoe_update == "direct") {
      apply_high_impact_probability(probability_prs, genotype, reference)
    } else {
      high_impact_method_probability(
        prs_curve, genotype, reference, K = K, SP = SP,
        method = probability_method, r2_liability = r2_liability
      )
    }
    data.frame(
      PRS = prs_curve,
      Probability = probability,
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
  if (isTRUE(show_rv_labels) && !is.null(rv_labels) &&
      (length(rv_labels) != n || anyNA(rv_labels))) {
    stop("rv_labels must contain one non-missing label per carrier")
  }
  groups <- .prepare_tiger_carrier_group(carrier_group, n, group_colours)
  .check_tiger_point_style(
    rv_point_size, rv_point_alpha, rv_point_border, rv_point_stroke,
    rv_shapes, c("1 RV", "2+ RVs")
  )
  carriers <- data.frame(
    PRS = carrier_prs,
    APOE = carrier_apoe,
    Probability_before_RV = carrier_probability_before_rv,
    Probability_after_RV = carrier_probability_after_rv,
    RV_count_group = factor(ifelse(rv_count == 1, "1 RV", "2+ RVs"),
                            levels = c("1 RV", "2+ RVs")),
    Carrier_group = groups$group,
    RV_label = if (is.null(rv_labels)) "" else as.character(rv_labels)
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
  plot <- ggplot2::ggplot(
    curves,
    ggplot2::aes(PRS, Probability, colour = Genotype, linetype = Genotype)
  ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(
      data = carriers,
      ggplot2::aes(PRS, Probability_after_RV,
                   shape = RV_count_group, fill = Carrier_group),
      inherit.aes = FALSE, colour = rv_point_border,
      alpha = rv_point_alpha,
      size = rv_point_size, stroke = rv_point_stroke
    ) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_linetype_manual(values = genotype_lines) +
    ggplot2::scale_shape_manual(values = rv_shapes) +
    ggplot2::scale_fill_manual(
      values = groups$colours,
      guide = if (groups$show_legend) "legend" else "none"
    ) +
    ggplot2::coord_cartesian(xlim = range(prs_curve), ylim = y_limits) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1),
      shape = ggplot2::guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        override.aes = .tiger_rv_shape_legend
      ),
      fill = ggplot2::guide_legend(
        order = 3, nrow = 1, byrow = TRUE,
        override.aes = .tiger_group_legend
      )
    ) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = "APOE genotype", linetype = "APOE genotype",
      fill = if (groups$show_legend) "Carrier group" else NULL,
      shape = "Number of RVs",
      title = paste0("APOE curves with RV carriers", method_text),
      subtitle = paste0(
        "Curves: PRS + APOE; points: PRS + APOE + RV (carriers only)"
      )
    ) +
    tiger_plot_theme()
  if (isTRUE(show_rv_labels) && !is.null(rv_labels)) {
    plot <- .add_tiger_rv_labels(plot, carriers, "Probability_after_RV")
  }
  plot
}
