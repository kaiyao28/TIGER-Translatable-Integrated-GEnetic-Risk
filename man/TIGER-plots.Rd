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
Consistently themed plotting functions for PRS probability curves, APOE or
other common high-impact variants, and adjusted points for RV carriers.
Optional PRS group labels draw separate coloured distributions as unconnected
hollow circles, with filled RV-adjusted carrier points retaining the
corresponding group colour.
Hollow PRS and genotype-distribution points consistently default to size
\code{1.25} and stroke \code{0.65}. Their opacity is reduced in denser
genotype panels. Solid
RV-carrier points are intentionally larger and fully opaque by default.
RV-carrier shape distinguishes one from multiple carried variants. Point fill
is reserved for optional user-defined carrier groups and does not encode
damaging or protective effect direction.
The RV plotting functions expose point size, alpha, border, stroke, shape,
group-colour and optional RV-label arguments directly. Labels are drawn when
\code{rv_labels} is supplied and can be disabled with
\code{show_rv_labels = FALSE}. The functions return standard \code{ggplot}
objects for further editing. See \code{docs/PLOTTING.md} and
\code{R/plotting.R}.
APOE reference distributions and direct labels use a shared warm-to-cool
genotype palette and accept a named \code{apoe_colours} vector for user-defined
shades.
For \code{plot_tiger_apoe_curves()}, supplying the four \code{sample_*}
arguments displays individual genotype-specific probabilities as hollow points.
Their colours represent the sample group, such as blue controls and red cases,
while direct labels identify the six APOE distributions.
\code{plot_tiger_high_impact_curves()} applies the same dense hollow-point
distribution display to a generic mutually exclusive genotype reference,
including a single 0/1/2 biallelic variant.
\code{plot_tiger_high_impact_rv_carrier_points()} adds the same labelled RV
carrier overlay used by the APOE plot.
For the combined APOE and RV plot, supplying \code{prs_group} and
\code{prs_apoe} shows observed genotype-specific distributions in one panel
using unconnected hollow points. Group determines colour. APOE distributions
are labelled directly, and filled RV-carrier points retain the group colour.
Genotype labels can be hidden with \code{show_genotype_labels = FALSE} or
resized with \code{genotype_label_size}.
\code{ggplot2} is required only when these functions are called.
}
