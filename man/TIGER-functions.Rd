\name{TIGER-functions}
\alias{effective_sample_size}
\alias{bpc_liability_scale_factor}
\alias{prepare_bpc_inputs}
\alias{prepare_sbayesr_summary_statistics}
\alias{observed_to_liability_r2}
\alias{liability_to_observed_r2}
\alias{bpc_probability}
\alias{pain_probability}
\alias{pair_probability_summary}
\alias{pair_probability_sample}
\alias{pair_probability}
\alias{intrinsic_rv_probability}
\alias{combine_independent_rv_probabilities}
\alias{apply_rv_probability}
\alias{prepare_rv_reference}
\alias{prepare_rv_carrier_matrix}
\alias{apply_rv_carriers}
\alias{prepare_high_impact_reference}
\alias{apoe_genotype_reference}
\alias{apply_high_impact_probability}
\alias{harmonise_scz_reference}
\alias{harmonise_ad_reference}
\alias{read_tiger_reference}
\alias{tiger_example_file}
\alias{tiger_reference_file}
\title{TIGER probability and input functions}
\description{
Public functions for preparing liability-scale PRS inputs, estimating absolute
probabilities, validating APOE or other high-impact references, and applying
rare-variant carrier effects. Each implementation is available as plain R
source in the package's \code{R/} directory.
}
\details{
Use \code{pair_probability_summary()} as the default worked-example method.
Read the repository's docs/GETTING_STARTED.md, docs/USAGE.md and docs/METHODS.md before applying
the functions to study data. In particular, exclude any separately modelled
region from the PRS and use compatible liability-scale inputs.
}
