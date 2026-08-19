\name{TIGER-functions}
\alias{effective_sample_size}
\alias{liability_prs_scale_factor}
\alias{prepare_liability_prs_inputs}
\alias{prepare_sbayesr_summary_statistics}
\alias{observed_to_liability_r2}
\alias{liability_to_observed_r2}
\alias{bpc_probability}
\alias{genopred_probability}
\alias{pair_probability_summary}
\alias{pair_probability_sample}
\alias{intrinsic_rv_probability}
\alias{combine_independent_rv_probabilities}
\alias{apply_rv_probability}
\alias{prepare_rv_reference}
\alias{prepare_rv_carrier_matrix}
\alias{apply_rv_carriers}
\alias{harmonise_scz_reference}
\alias{harmonise_ad_reference}
\alias{read_tiger_reference}
\title{TIGER probability and input functions}
\description{
Public functions for preparing liability-scale PRS inputs, estimating absolute
probabilities, preparing input references, and applying rare-variant carrier
effects. Each implementation is available as plain R
source in the package's \code{R/} directory.
Included example inputs are synthetic or hypothetical. TIGER does not
distribute study-specific or literature-derived reference datasets.
}
\details{
Use \code{pair_probability_summary()} as the default worked-example method.
Read the repository's docs/GETTING_STARTED.md, docs/USAGE.md and docs/METHODS.md before applying
the functions to study data. In particular, exclude any separately modelled
region from the PRS and use compatible liability-scale inputs.
Use \code{genopred_probability()} for the GenoPred conversion.

For \code{prepare_rv_reference()}, each row requires a unique \code{RV_IDs} and
either a positive \code{OR} or frequencies from which to calculate it. The
calculation is \code{[Case_freq/(1-Case_freq)] /
[Control_freq/(1-Control_freq)]}. A supplied \code{Population_freq} replaces a
missing \code{Control_freq}. The returned source columns record these choices.
\code{RV_IDs} is matched to the individual carrier-status field. \code{ID},
\code{RV_ID}, and \code{Variant_ID} are accepted aliases. \code{Symbol} and
\code{Class} are optional annotations. Column matching is case-insensitive.
}
