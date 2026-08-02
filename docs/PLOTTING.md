# Customising TIGER figures

TIGER plotting functions return ordinary `ggplot` objects. Figures can be
modified through function arguments or by adding standard `ggplot2` layers.
The complete implementation is available in [`R/plotting.R`](../R/plotting.R).

## RV-carrier points

The RV distribution functions expose the main point controls directly:

```r
p <- plot_tiger_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_curve = p_sequence,
  carrier_prs = carrier_prs,
  carrier_probability_before = probability_before,
  carrier_probability_after = probability_after,
  rv_count = rv_count,
  rv_point_size = 2.6,
  rv_point_alpha = 0.40,
  rv_point_border = "grey15",
  rv_point_stroke = 0.45,
  rv_shapes = c("1 RV" = 21, "2+ RVs" = 24)
)
```

Shape should continue to represent RV count. Fill is neutral by default. To
colour external groups, supply one label per carrier and a named colour vector:

```r
p <- plot_tiger_rv_carrier_points(
  # required probability arguments omitted here
  carrier_group = study_group,
  group_colours = c("Group A" = "#3182BD", "Group B" = "#E31A1C")
)
```

The plotted points use `rv_point_alpha` so group colours remain mild and do not
obscure the probability curves. The carrier-group legend deliberately uses the
same colours at full opacity, making group membership easy to identify. With a
single ungrouped cohort, TIGER uses one neutral fill and suppresses the
unnecessary carrier-group legend. For `plot_tiger_before_after_rv()`, provide
one `carrier_group` value per input individual; the carrier-only plotting
functions require one value per plotted carrier.

Do not use group colour to encode damaging versus protective variants. Their
direction is represented by displacement above or below the baseline curve.

## Standard ggplot2 editing

Titles, fonts, axes, legends and panel styling can be changed after creation:

```r
p +
  ggplot2::labs(title = "Study-specific title") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12),
    legend.position = "right"
  )
```

The same point arguments are available in `plot_tiger_before_after_rv()` and
`plot_tiger_apoe_rv_carrier_points()`. Advanced users can inspect layers with
`p$layers` and replace scales or themes using ordinary `ggplot2` syntax.
