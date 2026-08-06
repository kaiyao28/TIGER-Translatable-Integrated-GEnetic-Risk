# TIGER: Translatable Integrated Genetic Risk.
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
.tiger_prs_point_size <- 1.25
.tiger_prs_point_stroke <- 0.65
.tiger_apoe_colours <- c(
  "e4/e4" = "#B2182B", "e3/e4" = "#EF8A62", "e2/e4" = "#F4A582",
  "e3/e3" = "#737373", "e2/e3" = "#67A9CF", "e2/e2" = "#2166AC"
)
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
  x_span <- range(data$PRS, finite = TRUE)
  x_margin <- diff(x_span) * 0.12
  if (!is.finite(x_margin) || x_margin == 0) x_margin <- 0.1
  data$Label_hjust <- ifelse(
    data$PRS >= x_span[2] - x_margin, 1,
    ifelse(data$PRS <= x_span[1] + x_margin, 0, 0.5)
  )
  label_y <- data[[y_column]]
  alternate <- seq_len(nrow(data)) %% 2L
  data$Label_y <- ifelse(
    label_y >= 0.94,
    label_y - ifelse(alternate == 0L, 0.035, 0.075),
    ifelse(
      label_y <= 0.06,
      label_y + ifelse(alternate == 0L, 0.035, 0.075),
      label_y + ifelse(alternate == 0L, -0.04, 0.04)
    )
  )
  mapping <- ggplot2::aes(x = PRS)
  mapping$y <- as.name(y_column)
  mapping$label <- as.name("RV_label")
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot + ggrepel::geom_label_repel(
      data = data, mapping = mapping, inherit.aes = FALSE,
      size = 2.9, colour = "grey15", fill = "white",
      label.size = 0.18, label.padding = grid::unit(0.12, "lines"),
      box.padding = 0.45, point.padding = 0.25,
      force = 2, force_pull = 0.25,
      max.time = 2, max.iter = 50000, seed = 1049,
      min.segment.length = 0,
      segment.colour = "grey55", segment.size = 0.3,
      segment.alpha = 0.8, max.overlaps = Inf
    )
  } else {
    mapping$y <- as.name("Label_y")
    mapping$hjust <- as.name("Label_hjust")
    plot + ggplot2::geom_label(
      data = data, mapping = mapping, inherit.aes = FALSE,
      size = 2.9, colour = "grey15", fill = "white",
      linewidth = 0.18, label.padding = grid::unit(0.12, "lines")
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

# Compare the three PRS-to-probability approaches on the same liability-PRS
# values. Calculation stays in probability_methods.R; this helper only reshapes
# supplied results for a consistent plot.
plot_tiger_prs_methods <- function(
    prs, bpc_probability, genopred_probability, pair_summary_probability,
    y_limits = c(0, 1)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_prs_methods() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || length(prs) < 2L || any(!is.finite(prs))) {
    stop("prs must contain at least two finite values")
  }
  n <- length(prs)
  .check_plot_probability_vector(bpc_probability, "bpc_probability", n)
  .check_plot_probability_vector(
    genopred_probability, "genopred_probability", n
  )
  .check_plot_probability_vector(
    pair_summary_probability, "pair_summary_probability", n
  )
  method_levels <- c("BPC", "GenoPred", "PAIR (summary)")
  long <- data.frame(
    PRS = rep(prs, 3L),
    Probability = c(
      bpc_probability, genopred_probability, pair_summary_probability
    ),
    Method = factor(rep(method_levels, each = n), levels = method_levels)
  )
  long <- long[order(long$Method, long$PRS), , drop = FALSE]
  ggplot2::ggplot(
    long, ggplot2::aes(PRS, Probability, colour = Method, linetype = Method)
  ) +
    ggplot2::geom_point(
      shape = 21, fill = "white", size = .tiger_prs_point_size,
      stroke = 0.55, alpha = 0.95
    ) +
    ggplot2::scale_colour_manual(values = c(
      "BPC" = "#D73027", "GenoPred" = "#1A9850",
      "PAIR (summary)" = "#2C7BB6"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "BPC" = "dotted", "GenoPred" = "dashed",
      "PAIR (summary)" = "solid"
    )) +
    ggplot2::coord_cartesian(xlim = range(prs), ylim = y_limits) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = "Method", linetype = "Method",
      title = "PRS-to-probability conversion methods",
      subtitle = "The same liability-scale PRS values are converted by each method"
    ) +
    tiger_plot_theme()
}

# Show the OR-to-intrinsic-probability step that precedes individual RV updates.
plot_tiger_rv_reference <- function(rv_reference, prevalence = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_rv_reference() requires the ggplot2 package")
  }
  reference <- prepare_rv_reference(rv_reference)
  reference$Intrinsic_probability <- intrinsic_rv_probability(
    reference$OR, prevalence
  )
  reference$Effect <- factor(
    reference$Direction, levels = c("Protective", "Damaging")
  )
  reference$Display <- paste0(
    reference$Symbol, "  (OR ", format(round(reference$OR, 2), trim = TRUE), ")"
  )
  reference$Display <- factor(
    reference$Display,
    levels = reference$Display[order(reference$Intrinsic_probability)]
  )
  ggplot2::ggplot(
    reference,
    ggplot2::aes(Intrinsic_probability, Display, colour = Effect)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Intrinsic_probability,
                   y = Display, yend = Display),
      linewidth = 0.8, alpha = 0.65
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("p = %.3f", Intrinsic_probability)),
      hjust = -0.12, size = 3.2, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Protective" = "#3182BD", "Damaging" = "#E31A1C"
    )) +
    ggplot2::scale_x_continuous(
      limits = c(0, min(1, max(reference$Intrinsic_probability) * 1.22)),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::labs(
      x = "Intrinsic RV probability magnitude", y = NULL, colour = "RV effect",
      title = "Rare-variant probability components",
      subtitle = paste0(
        "Background prevalence = ", format(prevalence, trim = TRUE),
        ". Damaging raises probability; protective lowers it."
      )
    ) +
    tiger_plot_theme()
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
    apoe_update = c("method-specific", "direct"),
    apoe_colours = .tiger_apoe_colours,
    sample_prs = NULL, sample_probability = NULL,
    sample_genotype = NULL, sample_group = NULL,
    group_colours = c("Control" = "#3182BD", "Case" = "#E31A1C"),
    show_genotype_labels = TRUE, genotype_label_size = 3.1,
    plot_title = "APOE genotype probability distributions",
    plot_subtitle = "APOE is modelled separately from the PRS",
    legend_title = "APOE genotype") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_tiger_apoe_curves() requires the ggplot2 package")
  }
  if (!is.numeric(prs) || length(prs) < 2L || any(!is.finite(prs))) {
    stop("prs must contain at least two finite numeric values")
  }
  .check_plot_probability_vector(probability_prs, "probability_prs", length(prs))
  apoe_update <- match.arg(apoe_update)
  if (apoe_update == "method-specific" && is.null(K)) {
    stop("K is required for TIGER's method-specific high-impact update")
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
  if (is.null(names(apoe_colours)) || anyNA(apoe_colours) ||
      any(!nzchar(apoe_colours))) {
    stop("genotype colours must be a named vector of valid colour values")
  }
  genotype_colours <- apoe_colours
  missing_colours <- setdiff(genotype_order, names(apoe_colours))
  if (length(missing_colours)) {
    genotype_colours <- c(genotype_colours,
      stats::setNames(grDevices::hcl.colors(length(missing_colours), "Dark 3"),
                      missing_colours))
  }
  sample_inputs <- list(sample_prs, sample_probability,
                        sample_genotype, sample_group)
  supplied_samples <- !vapply(sample_inputs, is.null, logical(1))
  if (any(supplied_samples) && !all(supplied_samples)) {
    stop("sample_prs, sample_probability, sample_genotype, and sample_group must be supplied together")
  }
  show_samples <- all(supplied_samples)
  if (show_samples) {
    n_sample <- length(sample_prs)
    if (!is.numeric(sample_prs) || !n_sample || any(!is.finite(sample_prs))) {
      stop("sample_prs must contain finite values")
    }
    .check_plot_probability_vector(
      sample_probability, "sample_probability", n_sample
    )
    sample_genotype <- as.character(sample_genotype)
    sample_group <- as.character(sample_group)
    if (length(sample_genotype) != n_sample || anyNA(sample_genotype) ||
        any(!sample_genotype %in% genotype_order)) {
      stop("sample_genotype must contain one known genotype per sample")
    }
    if (length(sample_group) != n_sample || anyNA(sample_group) ||
        any(!nzchar(sample_group))) {
      stop("sample_group must contain one non-empty group per sample")
    }
    missing_groups <- setdiff(unique(sample_group), names(group_colours))
    if (is.null(names(group_colours)) || length(missing_groups)) {
      stop("group_colours must provide a named colour for every sample group")
    }
    samples <- data.frame(
      PRS = sample_prs, Probability = sample_probability,
      Genotype = factor(sample_genotype, levels = genotype_order),
      Group = factor(sample_group, levels = unique(sample_group))
    )
  }
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  plot <- ggplot2::ggplot(long, ggplot2::aes(PRS, Probability)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                       colour = "grey75", linewidth = 0.4)
  if (show_samples) {
    plot <- plot + ggplot2::geom_point(
      data = samples,
      ggplot2::aes(colour = Group),
      shape = 21, fill = "white", size = .tiger_prs_point_size,
      stroke = .tiger_prs_point_stroke, alpha = 0.75
    ) + ggplot2::scale_colour_manual(values = group_colours)
    if (isTRUE(show_genotype_labels)) {
      label_targets <- stats::setNames(
        seq(0.72, 0.32, length.out = length(genotype_order)), genotype_order
      )
      labels <- do.call(rbind, lapply(genotype_order, function(genotype) {
        curve <- long[long$Genotype == genotype, , drop = FALSE]
        curve[which.min(abs(curve$Probability - label_targets[[genotype]])), ]
      }))
      labels$Label_y <- pmin(y_limits[2] - 0.02,
                             labels$Probability + 0.055 * diff(y_limits))
      plot <- plot + ggplot2::geom_label(
        data = labels,
        ggplot2::aes(PRS, Label_y, label = Genotype),
        fill = unname(genotype_colours[as.character(labels$Genotype)]),
        colour = "white", fontface = "bold", size = genotype_label_size,
        linewidth = 0.25, label.padding = grid::unit(0.16, "lines"),
        show.legend = FALSE
      )
    }
  } else {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(colour = Genotype),
      shape = 21, fill = "white", size = .tiger_prs_point_size,
      stroke = 0.55, alpha = 0.95
    ) + ggplot2::scale_colour_manual(values = genotype_colours) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(shape = 21, fill = "white", size = 2.5, alpha = 1)
      ))
  }
  plot +
    ggplot2::coord_cartesian(xlim = range(prs), ylim = y_limits) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      colour = if (show_samples) "Group" else legend_title,
      title = paste0(plot_title, method_text),
      subtitle = plot_subtitle
    ) +
    tiger_plot_theme()
}

# Generic genotype curves for one common high-impact variant or another
# mutually exclusive genotype reference. This delegates to the same plotting
# and probability logic as the APOE-specific convenience function.
plot_tiger_high_impact_curves <- function(
    prs, probability_prs, high_impact_reference,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    K = NULL, SP = 0.5, r2_liability = NULL,
    high_impact_update = c("method-specific", "direct"),
    genotype_colours = NULL,
    plot_title = "Common high-impact genotype probability distributions") {
  reference <- prepare_high_impact_reference(high_impact_reference)
  if (is.null(genotype_colours)) {
    genotype_colours <- stats::setNames(
      grDevices::hcl.colors(nrow(reference), "Dark 3"), reference$Genotype
    )
  }
  plot_tiger_apoe_curves(
    prs = prs, probability_prs = probability_prs,
    apoe_reference = reference,
    probability_method = probability_method, y_limits = y_limits,
    K = K, SP = SP, r2_liability = r2_liability,
    apoe_update = match.arg(high_impact_update),
    apoe_colours = genotype_colours,
    plot_title = plot_title,
    plot_subtitle = "The high-impact component is modelled separately from the PRS",
    legend_title = "Genotype"
  ) + ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, fill = "white", size = 2.5, alpha = 1)
    )
  )
}

# Instructional/observed RV-carrier overlay. Optional prs_group values draw
# group-specific unadjusted PRS probability curves and observations. RV carrier
# points use the matching group colour. Shape distinguishes RV count.
plot_tiger_rv_carrier_points <- function(
    prs_curve, probability_curve,
    carrier_prs, carrier_probability_before, carrier_probability_after,
    rv_count, prs_group = NULL, carrier_group = NULL, group_colours = NULL,
    show_rv_labels = TRUE, rv_labels = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    show_prs_points = TRUE, prs_point_size = .tiger_prs_point_size,
    prs_point_alpha = 0.80,
    show_shift_segments = FALSE, rv_point_size = 2.6,
    rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey15", rv_point_stroke = 0.55,
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
  grouped_curves <- !is.null(prs_group)
  if (grouped_curves) {
    prs_group <- as.character(prs_group)
    if (length(prs_group) != length(prs_curve) || anyNA(prs_group) ||
        any(!nzchar(trimws(prs_group)))) {
      stop("prs_group must contain one non-empty group label per PRS value")
    }
    if (is.null(carrier_group)) {
      stop("carrier_group is required when prs_group is supplied")
    }
  }
  if (length(prs_point_size) != 1L || !is.finite(prs_point_size) ||
      prs_point_size <= 0 || length(prs_point_alpha) != 1L ||
      !is.finite(prs_point_alpha) || prs_point_alpha < 0 ||
      prs_point_alpha > 1) {
    stop("PRS point size must be positive and alpha must be in [0, 1]")
  }
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
  group_input <- if (grouped_curves) {
    c(prs_group, as.character(carrier_group))
  } else {
    carrier_group
  }
  groups <- .prepare_tiger_carrier_group(
    group_input,
    if (grouped_curves) length(prs_curve) + n else n,
    group_colours
  )
  curve_group <- if (grouped_curves) {
    groups$group[seq_along(prs_curve)]
  } else {
    factor(rep("PRS", length(prs_curve)))
  }
  carrier_group_factor <- if (grouped_curves) {
    groups$group[length(prs_curve) + seq_len(n)]
  } else {
    groups$group
  }
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
    Carrier_group = carrier_group_factor,
    RV_label = if (is.null(rv_labels)) "" else as.character(rv_labels)
  )
  method_text <- if (is.null(probability_method) ||
                     !nzchar(trimws(probability_method))) "" else
    paste0(": ", trimws(probability_method))
  curve_data <- data.frame(
    PRS = prs_curve, Probability = probability_curve,
    PRS_group = curve_group
  )
  curve_data <- curve_data[
    order(curve_data$PRS_group, curve_data$PRS), , drop = FALSE
  ]
  plot <- ggplot2::ggplot(curve_data, ggplot2::aes(PRS, Probability))
  if (grouped_curves) {
    if (isTRUE(show_prs_points)) {
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(colour = PRS_group),
        shape = 21, fill = NA, stroke = .tiger_prs_point_stroke,
        size = prs_point_size, alpha = prs_point_alpha
      )
    }
    plot <- plot + ggplot2::scale_colour_manual(values = groups$colours)
  } else {
    plot <- plot + ggplot2::geom_line(colour = "#3182BD", linewidth = 1.1)
  }
  plot <- plot +
    ggplot2::coord_cartesian(xlim = range(prs_curve), ylim = y_limits) +
    ggplot2::scale_shape_manual(values = rv_shapes) +
    ggplot2::scale_fill_manual(
      values = groups$colours,
      guide = if (groups$show_legend && !grouped_curves) "legend" else "none"
    ) +
    ggplot2::guides(
      colour = if (grouped_curves) {
        ggplot2::guide_legend(
          order = 1, nrow = 1, byrow = TRUE,
          override.aes = list(shape = 21, fill = NA, alpha = 1, size = 3)
        )
      } else {
        "none"
      },
      shape = ggplot2::guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        override.aes = .tiger_rv_shape_legend
      ),
      fill = if (grouped_curves) {
        "none"
      } else {
        ggplot2::guide_legend(
          order = 3, nrow = 1, byrow = TRUE,
          override.aes = .tiger_group_legend
        )
      }
    ) +
    ggplot2::labs(
      x = "Liability PRS", y = "Estimated disorder probability",
      shape = "Number of RVs",
      colour = if (grouped_curves) "Group" else NULL,
      fill = if (groups$show_legend && !grouped_curves) "Carrier group" else NULL,
      title = paste0("RV carrier probabilities", method_text),
      subtitle = if (grouped_curves) {
        "Hollow circles: PRS only by group; filled points: RV-adjusted carriers"
      } else {
        "The curve is PRS only; points are shown only for RV carriers"
      }
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
    rv_count, prs_group = NULL, prs_apoe = NULL,
    carrier_group = NULL, group_colours = NULL,
    show_rv_labels = TRUE, rv_labels = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    K = NULL, SP = 0.5, r2_liability = NULL,
    apoe_update = c("method-specific", "direct"),
    apoe_colours = .tiger_apoe_colours,
    show_genotype_labels = TRUE, genotype_label_size = 3.1,
    show_prs_points = TRUE, prs_point_size = .tiger_prs_point_size,
    prs_point_alpha = 0.65,
    rv_point_size = 2.1, rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey25", rv_point_stroke = 0.35,
    rv_shapes = .tiger_rv_shapes[c("1 RV", "2+ RVs")],
    grouped_title = "APOE distributions with RV carriers",
    grouped_subtitle = "Hollow circles: PRS + APOE; solid points: PRS + APOE + RV",
    curve_title = "APOE curves with RV carriers",
    curve_subtitle = "Curves: PRS + APOE; points: PRS + APOE + RV (carriers only)",
    genotype_legend_title = "APOE genotype") {
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
  if (!is.logical(show_genotype_labels) || length(show_genotype_labels) != 1L ||
      is.na(show_genotype_labels)) {
    stop("show_genotype_labels must be TRUE or FALSE")
  }
  if (!is.numeric(genotype_label_size) || length(genotype_label_size) != 1L ||
      !is.finite(genotype_label_size) || genotype_label_size <= 0) {
    stop("genotype_label_size must be one positive finite value")
  }
  grouped_points <- !is.null(prs_group)
  if (grouped_points) {
    prs_group <- as.character(prs_group)
    if (length(prs_group) != length(prs_curve) || anyNA(prs_group) ||
        any(!nzchar(trimws(prs_group)))) {
      stop("prs_group must contain one non-empty group label per PRS value")
    }
    if (is.null(carrier_group)) {
      stop("carrier_group is required when prs_group is supplied")
    }
    prs_apoe <- as.character(prs_apoe)
    if (length(prs_apoe) != length(prs_curve) || anyNA(prs_apoe) ||
        any(!nzchar(trimws(prs_apoe)))) {
      stop("prs_apoe must contain one non-empty APOE genotype per PRS value")
    }
  }
  preferred <- c("e4/e4", "e3/e4", "e2/e4", "e3/e3", "e2/e3", "e2/e2")
  genotype_order <- c(preferred[preferred %in% reference$Genotype],
                      setdiff(reference$Genotype, preferred))
  if (grouped_points && any(!prs_apoe %in% genotype_order)) {
    stop("prs_apoe contains a genotype absent from apoe_reference")
  }
  curve_rows <- lapply(genotype_order, function(genotype) {
    keep <- if (grouped_points) prs_apoe == genotype else
      rep(TRUE, length(prs_curve))
    probability <- if (apoe_update == "direct") {
      apply_high_impact_probability(probability_prs, genotype, reference)
    } else {
      high_impact_method_probability(
        prs_curve, genotype, reference, K = K, SP = SP,
        method = probability_method, r2_liability = r2_liability
      )
    }
    data.frame(
      PRS = prs_curve[keep],
      Probability = probability[keep],
      Genotype = rep(genotype, sum(keep)),
      PRS_group = if (grouped_points) prs_group[keep] else
        rep(genotype, sum(keep))
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
  group_input <- if (grouped_points) {
    c(prs_group, as.character(carrier_group))
  } else {
    carrier_group
  }
  groups <- .prepare_tiger_carrier_group(
    group_input, if (grouped_points) length(prs_curve) + n else n,
    group_colours
  )
  carrier_groups <- if (grouped_points) {
    groups$group[length(prs_curve) + seq_len(n)]
  } else {
    groups$group
  }
  .check_tiger_point_style(
    rv_point_size, rv_point_alpha, rv_point_border, rv_point_stroke,
    rv_shapes, c("1 RV", "2+ RVs")
  )
  carriers <- data.frame(
    PRS = carrier_prs,
    Genotype = factor(carrier_apoe, levels = genotype_order),
    Probability_before_RV = carrier_probability_before_rv,
    Probability_after_RV = carrier_probability_after_rv,
    RV_count_group = factor(ifelse(rv_count == 1, "1 RV", "2+ RVs"),
                            levels = c("1 RV", "2+ RVs")),
    Carrier_group = carrier_groups,
    RV_label = if (is.null(rv_labels)) "" else as.character(rv_labels)
  )
  if (is.null(names(apoe_colours)) || anyNA(apoe_colours) ||
      any(!nzchar(apoe_colours))) {
    stop("apoe_colours must be a named vector of valid colour values")
  }
  genotype_colours <- apoe_colours
  missing_colours <- setdiff(genotype_order, names(apoe_colours))
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
  if (grouped_points) {
    reference_prs <- seq(min(prs_curve), max(prs_curve), length.out = 200L)
    reference_probability <- if (apoe_update == "direct") {
      stats::approx(prs_curve, probability_prs, xout = reference_prs,
                    ties = mean, rule = 2)$y
    } else {
      NULL
    }
    reference_curves <- do.call(rbind, lapply(genotype_order, function(genotype) {
      probability <- if (apoe_update == "direct") {
        apply_high_impact_probability(reference_probability, genotype, reference)
      } else {
        high_impact_method_probability(
          reference_prs, genotype, reference, K = K, SP = SP,
          method = probability_method, r2_liability = r2_liability
        )
      }
      data.frame(
        PRS = reference_prs, Probability = probability,
        Genotype = factor(genotype, levels = genotype_order)
      )
    }))
    label_targets <- stats::setNames(
      seq(0.65, 0.35, length.out = length(genotype_order)), genotype_order
    )
    genotype_labels <- do.call(rbind, lapply(genotype_order, function(genotype) {
      genotype_curve <- reference_curves[
        reference_curves$Genotype == genotype, , drop = FALSE
      ]
      genotype_curve[
        which.min(abs(genotype_curve$Probability - label_targets[[genotype]])),
        , drop = FALSE
      ]
    }))
    probability_span <- diff(y_limits)
    genotype_labels$Label_y <- pmin(
      y_limits[2] - 0.02 * probability_span,
      genotype_labels$Probability + 0.055 * probability_span
    )
    genotype_label_fills <- unname(
      genotype_colours[as.character(genotype_labels$Genotype)]
    )
    plot <- ggplot2::ggplot(curves, ggplot2::aes(PRS, Probability))
    if (isTRUE(show_prs_points)) {
      plot <- plot + ggplot2::geom_point(
        shape = 21, fill = "white",
        colour = unname(groups$colours[as.character(curves$PRS_group)]),
        size = prs_point_size, alpha = prs_point_alpha,
        stroke = .tiger_prs_point_stroke
      )
    }
    if (isTRUE(show_genotype_labels)) {
      plot <- plot + ggplot2::geom_label(
        data = genotype_labels,
        ggplot2::aes(x = PRS, y = Label_y, label = Genotype),
        inherit.aes = FALSE, fill = genotype_label_fills, colour = "white",
        linewidth = 0.25, label.padding = grid::unit(0.16, "lines"),
        size = genotype_label_size, fontface = "bold", show.legend = FALSE
      )
    }
    plot <- plot + ggplot2::geom_point(
      data = carriers,
      ggplot2::aes(PRS, Probability_after_RV,
                   fill = Carrier_group,
                   shape = RV_count_group),
      inherit.aes = FALSE,
      colour = rv_point_border,
      alpha = rv_point_alpha, size = rv_point_size,
      stroke = rv_point_stroke
    ) +
      ggplot2::scale_fill_manual(values = groups$colours) +
      ggplot2::scale_shape_manual(values = rv_shapes) +
      ggplot2::coord_cartesian(xlim = range(prs_curve), ylim = y_limits) +
      ggplot2::guides(
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
        colour = NULL, fill = "Carrier group",
        shape = "Number of RVs",
        title = paste0(grouped_title, method_text),
        subtitle = grouped_subtitle
      ) +
      tiger_plot_theme()
    if (isTRUE(show_rv_labels) && !is.null(rv_labels)) {
      plot <- .add_tiger_rv_labels(plot, carriers, "Probability_after_RV")
    }
    return(plot)
  }
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
      alpha = rv_point_alpha, size = rv_point_size, stroke = rv_point_stroke
    ) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_linetype_manual(values = genotype_lines) +
    ggplot2::scale_shape_manual(values = rv_shapes) +
    ggplot2::scale_fill_manual(values = groups$colours,
                               guide = if (groups$show_legend) "legend" else "none") +
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
      colour = genotype_legend_title, linetype = genotype_legend_title,
      fill = if (groups$show_legend) "Carrier group" else NULL,
      shape = "Number of RVs",
      title = paste0(curve_title, method_text),
      subtitle = curve_subtitle
    ) +
    tiger_plot_theme()
  if (isTRUE(show_rv_labels) && !is.null(rv_labels)) {
    plot <- .add_tiger_rv_labels(plot, carriers, "Probability_after_RV")
  }
  plot
}

# Generic high-impact genotype + RV plot. The point and probability logic is
# identical to the APOE-specific plot, while public argument names and labels
# apply to one biallelic 0/1/2 variant or another mutually exclusive reference.
plot_tiger_high_impact_rv_carrier_points <- function(
    prs_curve, probability_prs, high_impact_reference,
    carrier_prs, carrier_genotype,
    carrier_probability_before_rv, carrier_probability_after_rv,
    rv_count, prs_group = NULL, prs_genotype = NULL,
    carrier_group = NULL, group_colours = NULL,
    show_rv_labels = TRUE, rv_labels = NULL,
    probability_method = "PAIR (summary)", y_limits = c(0, 1),
    K = NULL, SP = 0.5, r2_liability = NULL,
    high_impact_update = c("method-specific", "direct"),
    genotype_colours = NULL,
    show_genotype_labels = TRUE, genotype_label_size = 3.1,
    show_prs_points = TRUE, prs_point_size = .tiger_prs_point_size,
    prs_point_alpha = 0.65,
    rv_point_size = 2.1, rv_point_alpha = .tiger_rv_point_alpha,
    rv_point_border = "grey25", rv_point_stroke = 0.35,
    rv_shapes = .tiger_rv_shapes[c("1 RV", "2+ RVs")]) {
  reference <- prepare_high_impact_reference(high_impact_reference)
  if (is.null(genotype_colours)) {
    genotype_colours <- stats::setNames(
      grDevices::hcl.colors(nrow(reference), "Dark 3"), reference$Genotype
    )
  }
  plot_tiger_apoe_rv_carrier_points(
    prs_curve = prs_curve, probability_prs = probability_prs,
    apoe_reference = reference,
    carrier_prs = carrier_prs, carrier_apoe = carrier_genotype,
    carrier_probability_before_rv = carrier_probability_before_rv,
    carrier_probability_after_rv = carrier_probability_after_rv,
    rv_count = rv_count, prs_group = prs_group, prs_apoe = prs_genotype,
    carrier_group = carrier_group, group_colours = group_colours,
    show_rv_labels = show_rv_labels, rv_labels = rv_labels,
    probability_method = probability_method, y_limits = y_limits,
    K = K, SP = SP, r2_liability = r2_liability,
    apoe_update = match.arg(high_impact_update), apoe_colours = genotype_colours,
    show_genotype_labels = show_genotype_labels,
    genotype_label_size = genotype_label_size,
    show_prs_points = show_prs_points, prs_point_size = prs_point_size,
    prs_point_alpha = prs_point_alpha, rv_point_size = rv_point_size,
    rv_point_alpha = rv_point_alpha, rv_point_border = rv_point_border,
    rv_point_stroke = rv_point_stroke, rv_shapes = rv_shapes,
    grouped_title = "High-impact genotype distributions with RV carriers",
    grouped_subtitle = paste0(
      "Hollow circles: PRS + high-impact genotype; solid points: ",
      "PRS + high-impact genotype + RV"
    ),
    curve_title = "High-impact genotype curves with RV carriers",
    curve_subtitle = paste0(
      "Curves: PRS + high-impact genotype; points: ",
      "PRS + high-impact genotype + RV (carriers only)"
    ),
    genotype_legend_title = "Genotype"
  )
}
