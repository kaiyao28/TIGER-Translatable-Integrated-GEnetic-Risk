\name{TIGER-workflow}
\alias{tiger_probabilities}
\alias{tiger_default_apoe_reference}
\title{Calculate TIGER probabilities from one merged individual table}
\description{
Calculate PRS probabilities and optionally append rare-variant, APOE or other
common high-impact-variant, and combined probabilities while preserving the
input rows and columns.
}
\usage{
tiger_probabilities(
  data, K, SP = 0.5,
  method = c("PAIR (summary)", "BPC", "GenoPred", "PAIR (sample)"),
  id_col = "ID", prs_col = "PRS_liability", group_col = NULL,
  r2_liability = NULL, reference_prs_liability = NULL,
  r2_observed = NULL, case_mean = NULL, case_sd = NULL,
  control_mean = NULL, control_sd = NULL, n_quantiles = 100,
  include_rv = FALSE, rv_reference = NULL,
  rv_columns = NULL, rv_status_col = NULL,
  rv_delimiter = "[;,|]", rv_prevalence = 0.5,
  include_apoe = FALSE, apoe_col = "APOE", apoe_reference = NULL
)

tiger_default_apoe_reference()
}
\arguments{
  \item{data}{Non-empty data frame with one row per individual.}
  \item{K}{Population disorder prevalence.}
  \item{SP}{Target-sample prevalence.}
  \item{method}{PRS-to-probability method.}
  \item{id_col, prs_col, group_col, apoe_col}{Configured individual-table column names.}
  \item{r2_liability, reference_prs_liability, r2_observed}{Method inputs.}
  \item{case_mean, case_sd, control_mean, control_sd}{PAIR (sample) inputs.}
  \item{n_quantiles}{Number of GenoPred probability quantiles.}
  \item{include_rv, include_apoe}{Logical flags enabling optional layers.}
  \item{rv_reference}{RV effect reference accepted by \code{prepare_rv_reference()}.}
  \item{rv_columns}{Optional named character vector mapping RV-reference IDs to 0/1 or logical columns.}
  \item{rv_status_col}{Optional column containing delimited RV-reference IDs.}
  \item{rv_delimiter}{Regular expression separating RV IDs.}
  \item{rv_prevalence}{Background probability used to translate RV odds ratios.}
  \item{apoe_reference}{High-impact genotype reference. If omitted, the bundled illustrative APOE reference is used with a message.}
}
\value{
The input data frame with aligned probability and optional RV-component columns
appended. Reference tables used by the calculation are retained as attributes.
}
\details{
The PRS must already be on the liability scale. The RV and APOE layers default
to off. RV IDs must match the RV reference. APOE labels must match the genotype
reference. The combined condition recalculates the selected PRS method using
genotype-specific prevalence parameters before applying the same RV operator.
The bundled APOE frequencies are illustrative and should be replaced when an
appropriate population- and phenotype-matched reference is available.
}
