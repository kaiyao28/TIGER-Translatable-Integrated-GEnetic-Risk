# Customising TIGER figures

TIGER plotting functions return ordinary `ggplot` objects. Figures can be
modified through function arguments or by adding standard `ggplot2` layers.
The complete implementation is available in [`R/plotting.R`](../R/plotting.R).

Calculation and plotting are deliberately separated. Probability functions
produce the values; plotting helpers only reshape and display supplied results.
`plot_tiger_prs_methods()` compares BPC, GenoPred, and PAIR (summary) on the
same liability-scale PRS values. `plot_tiger_rv_reference()` displays the
OR-to-intrinsic-probability step used before individual RV updates.

## RV-carrier points

The RV distribution functions expose the main point controls directly:

```r
p <- plot_tiger_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_curve = p_sequence,
  prs_group = sample_group,
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
| `rv_point_size` | `1.25` or `2.6`, depending on the plot | Point size |
| `rv_point_alpha` | `1` | Point opacity from 0 to 1 |
| `rv_point_border` | `"grey15"` or `"grey25"`, depending on the plot | Point-outline colour |
| `rv_point_stroke` | plot-specific | Point-outline width |
| `rv_shapes` | circle and triangle | Shapes for one and multiple RVs |
| `prs_group` | `NULL` | Optional group label for each PRS observation |
| `show_prs_points` | `TRUE` | Show unconnected hollow PRS observations by group |
| `prs_point_size` | `1.5` | Size of hollow PRS-only group observations |
| `prs_point_alpha` | `0.80` | Opacity of hollow PRS-only group observations |
| `carrier_group` | `NULL` | Optional external group for point colour |
| `group_colours` | automatic | Named colours for external groups |
| `show_rv_labels` | `TRUE` | Display supplied RV labels |
| `rv_labels` | `NULL` | One RV name or description per plotted carrier |

With `rv_labels = NULL`, no labels are drawn even though
`show_rv_labels = TRUE`. Labels therefore become automatic when names are
supplied without affecting calls that do not include names.

RV names are shown automatically when `rv_labels` is supplied. Set
`show_rv_labels = FALSE` to hide them. Supply one value in `rv_labels` for every
plotted carrier. Labels can contain a single variant name or a combined
description such as `"SORL1; TREM2"` for carriers with multiple RVs. When the
optional `ggrepel` package is installed, TIGER uses deterministic force-based
repulsion and leader lines to separate nearby labels. Without `ggrepel`, TIGER
uses a dependency-free staggered fallback. For very dense figures, hide the
labels and provide RV names in a table, or label only a prespecified subset.

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

To show separate group-specific PRS distributions, supply one `prs_group` per
value in `prs_curve` and the matching `carrier_group` for each RV carrier. The
same named colours are used for the hollow PRS observations and filled
adjusted carrier points:

```r
p <- plot_tiger_rv_carrier_points(
  # required probability arguments omitted here
  prs_group = sample_group,
  carrier_group = carrier_sample_group,
  group_colours = c("Control" = "#3182BD", "Case" = "#E31A1C")
)
```

Unconnected hollow circles show each group PRS distribution. Larger filled
points are RV-adjusted carriers. Shape continues to distinguish one from
multiple RVs. If `prs_group` is omitted, TIGER draws one PRS curve and can still
colour carrier points using `carrier_group`.

RV-adjusted points use `rv_point_alpha = 1` by default. PRS-only observations
are lighter through `prs_point_alpha`. With a single ungrouped cohort, TIGER
uses one neutral carrier fill and suppresses the unnecessary group legend.

Do not use group colour to encode damaging versus protective variants. Their
direction is represented by displacement above or below the baseline curve.

## APOE genotype colours

The APOE plots use coordinated shades by default. Genotypes containing e4 use
warm red shades, genotypes containing e2 use cool blue shades, and e3/e3 is
neutral grey. Both APOE plotting functions accept the same named override:

```r
apoe_colours <- c(
  "e4/e4" = "#B2182B", "e3/e4" = "#EF8A62", "e2/e4" = "#F4A582",
  "e3/e3" = "#737373", "e2/e3" = "#67A9CF", "e2/e2" = "#2166AC"
)

plot_tiger_apoe_curves(
  # required probability arguments omitted here
  apoe_colours = apoe_colours
)
```

When `plot_tiger_apoe_rv_carrier_points()` is called without `prs_group`, the
same `apoe_colours` vector preserves identical genotype colours across the two
line-based APOE figures. When `prs_group` is supplied, the combined plot instead
uses unconnected points in one panel. Supply each individual's observed APOE
genotype through `prs_apoe`. Colour identifies the group, while the six
genotype-specific distributions are separated by their probability trajectories
and labelled directly. Filled RV-carrier points retain the group colour.

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
