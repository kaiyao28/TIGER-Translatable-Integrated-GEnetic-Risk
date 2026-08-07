\name{TIGER-plots}
\alias{tiger_plot_theme}
\alias{plot_tiger_prs_methods}
\alias{plot_tiger_rv_reference}
\alias{plot_tiger_before_after_rv}
\alias{plot_tiger_before_after_apoe}
\alias{plot_tiger_apoe_curves}
\alias{plot_tiger_high_impact_curves}
\alias{plot_tiger_high_impact_rv_carrier_points}
\alias{plot_tiger_rv_carrier_points}
\alias{plot_tiger_apoe_rv_carrier_points}
\title{TIGER plotting functions}
\description{
Optional, consistently themed plots for PRS probability methods, RV reference
effects, common high-impact genotypes, and RV-adjusted carriers. Distributions
use hollow points. RV carriers use larger filled points, with shape indicating
one versus multiple RVs and colour retaining an optional sample group.
Damaging or protective direction is shown by the probability displacement,
not by colour.

RV labels use black text on a white background and may be hidden with
\code{show_rv_labels = FALSE}. APOE and generic high-impact genotype labels may
also be hidden or resized. Plot functions expose the main point, alpha, colour,
shape, and label controls and return ordinary \code{ggplot} objects for further
editing. See \code{docs/PLOTTING.md} for inputs and examples. \code{ggplot2} is
required only when a plot function is called; \code{ggrepel} is used when
available for label placement.
}
