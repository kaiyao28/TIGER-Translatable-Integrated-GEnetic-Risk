\name{TIGER-plots}
\alias{tiger_plot_theme}
\alias{plot_tiger_before_after_rv}
\alias{plot_tiger_before_after_apoe}
\alias{plot_tiger_apoe_curves}
\alias{plot_tiger_rv_carrier_points}
\alias{plot_tiger_apoe_rv_carrier_points}
\title{TIGER plotting functions}
\description{
Consistently themed plotting functions for PRS probability curves, APOE or
other common high-impact variants, and adjusted points for RV carriers.
RV-carrier shape distinguishes one from multiple carried variants. Point fill
is reserved for optional user-defined carrier groups and does not encode
damaging or protective effect direction.
The RV plotting functions expose point size, alpha, border, stroke and shape
arguments directly and return standard \code{ggplot} objects for further
editing. See \code{docs/PLOTTING.md} and \code{R/plotting.R}.
\code{ggplot2} is required only when these functions are called.
}
