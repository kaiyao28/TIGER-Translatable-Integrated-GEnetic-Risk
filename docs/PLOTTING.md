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
  show_rv_labels = TRUE,
  rv_labels = rv_name,
  rv_point_size = 2.6,
  rv_point_alpha = 0.40,
  rv_point_border = "grey15",
  rv_point_stroke = 0.45,
  rv_shapes = c("1 RV" = 21, "2+ RVs" = 24)
)
```

The principal settings are:

| Argument | Default | Purpose |
| --- | --- | --- |
| `rv_point_size` | `1.25` or `2.1`, depending on the plot | Point size |
| `rv_point_alpha` | `1` | Point opacity from 0 to 1 |
| `rv_point_border` | `"grey25"` | Point-outline colour |
| `rv_point_stroke` | plot-specific | Point-outline width |
| `rv_shapes` | circle and triangle | Shapes for one and multiple RVs |
| `carrier_group` | `NULL` | Optional external group for point colour |
| `group_colours` | automatic | Named colours for external groups |
| `show_rv_labels` | `TRUE` | Display supplied RV labels |
| `rv_labels` | `NULL` | One RV name or description per plotted carrier |

With `rv_labels = NULL`, no labels are drawn even though
`show_rv_labels = TRUE`. Labels therefore become automatic when names are
supplied without affecting calls that do not include names.

RV names are shown automatically when `rv_labels` is supplied. Set
`show_rv_labels = FALSE` to hide them. Supply one value in `rv_labels` for every
plotted carrier. Labels can contain a single variant name or a combined description
such as `"SORL1; TREM2"` for carriers with multiple RVs. Overlapping labels are
separated with short guide lines when the optional `ggrepel` package is
available. Otherwise, overlapping labels are suppressed automatically.

```r
# Names are shown and moved apart when necessary.
p <- plot_tiger_rv_carrier_points(
  # required probability arguments omitted here
  rv_count = rv_count,
  rv_labels = rv_name
)

# Hide names in a dense figure.
p <- plot_tiger_rv_carrier_points(
  # required probability arguments omitted here
  rv_count = rv_count,
  rv_labels = rv_name,
  show_rv_labels = FALSE
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

The plotted points use `rv_point_alpha` (default `1`) so group colours remain
mild without becoming difficult to see and do not
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
