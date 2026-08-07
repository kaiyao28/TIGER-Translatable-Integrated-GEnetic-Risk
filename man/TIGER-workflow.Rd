\name{TIGER-workflow}
\alias{tiger_probabilities}
\alias{tiger_default_apoe_reference}
\title{Calculate TIGER probabilities from one merged individual table}
\description{
Calculate PRS probabilities and optionally append rare-variant (RV), common
high-impact-genotype, and combined probabilities while preserving the input
rows and columns. APOE is supported as a six-genotype high-impact component.
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
  include_high_impact = FALSE, high_impact_col = NULL,
  high_impact_reference = NULL,
  include_apoe = FALSE, apoe_col = "APOE", apoe_reference = NULL
)

tiger_default_apoe_reference()
}
\arguments{
  \item{data}{Non-empty data frame with one row per individual.}
  \item{K}{Population disorder prevalence.}
  \item{SP}{Target-sample prevalence.}
  \item{method}{PRS-to-probability method.}
  \item{id_col, prs_col, group_col}{Configured ID, liability-scale PRS and optional grouping column names. The group column is retained for plotting and does not alter probability calculations.}
  \item{r2_liability, reference_prs_liability, r2_observed}{Method inputs.}
  \item{case_mean, case_sd, control_mean, control_sd}{PAIR (sample) inputs.}
  \item{n_quantiles}{Number of GenoPred probability quantiles.}
  \item{include_rv}{Logical flag enabling the optional RV layer.}
  \item{include_high_impact}{Logical flag enabling a generic separately modelled common high-impact genotype layer. Do not also enable \code{include_apoe}.}
  \item{high_impact_col}{Individual-table column containing mutually exclusive genotype labels.}
  \item{high_impact_reference}{Case/control genotype-frequency reference accepted by \code{prepare_high_impact_reference()}.}
  \item{include_apoe}{Logical flag enabling the six-genotype APOE component. Do not also enable \code{include_high_impact}.}
  \item{apoe_col}{Individual-table column containing APOE genotype labels.}
  \item{rv_reference}{RV effect reference accepted by \code{prepare_rv_reference()}.}
  \item{rv_columns}{Optional named character vector mapping RV-reference IDs to 0/1 or logical columns.}
  \item{rv_status_col}{Optional column containing delimited RV-reference IDs.}
  \item{rv_delimiter}{Regular expression separating RV IDs.}
  \item{rv_prevalence}{Disease-probability prior used to translate RV odds ratios. A value of 0.50 gives a balanced-sample update. This is not RV allele frequency.}
  \item{apoe_reference}{High-impact genotype reference. If omitted, the bundled illustrative APOE reference is used with a message.}
}
\value{
The input data frame with aligned probability and optional RV-component columns
appended. Reference tables used by the calculation are retained as attributes.
}
\details{
The PRS must already be on the liability scale. RV and high-impact layers
default to off. Any separately modelled common high-impact variant, or an
appropriate LD-aware region around it, must be excluded from both PRS
construction and the PRS R-squared estimate to avoid double counting. RV IDs
must match the RV reference. High-impact genotype labels must match their
reference. The combined condition recalculates the selected PRS method using
genotype-specific population prevalence and sample prevalence before applying
the same RV operator.
The generic high-impact route returns \code{Probability_PRS_HIGH_IMPACT} and,
when RVs are enabled, \code{Probability_PRS_HIGH_IMPACT_RV}. The APOE route
returns \code{Probability_PRS_APOE} and, with RVs,
\code{Probability_PRS_APOE_RV}.
The bundled APOE frequencies are illustrative and should be replaced when an
appropriate population- and phenotype-matched reference is available.
}
\examples{
# Begin with PRS only
individuals <- data.frame(
  ID = c("P001", "P002", "P003"),
  PRS_liability = c(-0.5, 0, 0.5),
  RV_status = c("", "RISK_A", "PROTECT_C"),
  APOE = c("e2/e3", "e3/e3", "e3/e4")
)

prs_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1
)

# Add RV and APOE as the high-impact component after PRS-only works
rv_reference <- data.frame(
  ID = c("RISK_A", "PROTECT_C"),
  Symbol = c("GENE_A", "GENE_C"),
  Class = c("PTV", "Protective"),
  Case_freq = c(0.010, 0.003),
  Control_freq = c(0.002, 0.008),
  OR = c(8.0, 0.4)
)

full_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_status", rv_prevalence = 0.5,
  include_apoe = TRUE, apoe_col = "APOE",
  apoe_reference = tiger_default_apoe_reference()
)

full_results[c(
  "ID", "Probability_PRS", "Probability_PRS_RV",
  "Probability_PRS_APOE", "Probability_PRS_APOE_RV"
)]

# The generic route also accepts a single 0/1/2 high-impact variant
single_variant_reference <- biallelic_genotype_reference(0.20, 0.10)
individuals$Variant_genotype <- c("0", "1", "2")
single_variant_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1,
  include_high_impact = TRUE,
  high_impact_col = "Variant_genotype",
  high_impact_reference = single_variant_reference
)
}
