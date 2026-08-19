\name{TIGER-workflow}
\alias{tiger_probabilities}
\alias{tiger_default_apoe_reference}
\alias{check_tiger_inputs}
\title{Validate inputs and calculate TIGER probabilities}
\description{
Calculate PRS probabilities and optionally append rare-variant (RV), common
high-impact-genotype, and combined probabilities while preserving the input
rows and columns. APOE is supported as a six-genotype high-impact component.
\code{check_tiger_inputs()} validates the same individual and reference fields
before calculation.
}
\usage{
tiger_probabilities(
  data, K, SP = 0.5,
  method = c("PAIR (summary)", "BPC", "GenoPred", "PAIR (sample)"),
  id_col = "ID", prs_col = "PRS_liability",
  r2_liability = NULL, reference_prs_liability = NULL,
  r2_observed = NULL, case_mean = NULL, case_sd = NULL,
  control_mean = NULL, control_sd = NULL, n_quantiles = 100,
  include_rv = FALSE, rv_reference = NULL,
  rv_columns = NULL, rv_status_col = NULL,
  rv_delimiter = "[;,|]", rv_prevalence = 0.5,
  include_high_impact = FALSE, high_impact_col = NULL,
  high_impact_reference = NULL,
  include_apoe = FALSE, apoe_col = "APOE_genotype", apoe_reference = NULL
)

tiger_default_apoe_reference()

check_tiger_inputs(
  data, id_col = "ID", prs_col = "PRS_liability",
  include_rv = FALSE, rv_reference = NULL,
  rv_columns = NULL, rv_status_col = NULL, rv_delimiter = "[;,|]",
  include_high_impact = FALSE, high_impact_col = NULL,
  high_impact_reference = NULL,
  include_apoe = FALSE, apoe_col = "APOE_genotype",
  apoe_reference = NULL, verbose = TRUE
)
}
\arguments{
  \item{data}{Non-empty data frame with one row per individual.}
  \item{K}{Population disorder prevalence.}
  \item{SP}{Target-sample prevalence.}
  \item{method}{PRS-to-probability method.}
  \item{id_col, prs_col}{Configured ID and liability-scale PRS column names. Other input columns, including any plotting group, are retained unchanged.}
  \item{r2_liability, reference_prs_liability, r2_observed}{Method inputs.}
  \item{case_mean, case_sd, control_mean, control_sd}{PAIR (sample) inputs.}
  \item{n_quantiles}{Number of GenoPred probability quantiles.}
  \item{include_rv}{Logical flag enabling the optional RV layer.}
  \item{include_high_impact}{Logical flag enabling a generic separately modelled common high-impact genotype layer. Do not also enable \code{include_apoe}.}
  \item{high_impact_col}{Individual-table column containing mutually exclusive
  genotype labels matching \code{High_impact_genotype} in the high-impact reference.}
  \item{high_impact_reference}{Case/control genotype-frequency reference accepted by \code{prepare_high_impact_reference()}.}
  \item{include_apoe}{Logical flag enabling the six-genotype APOE component. Do not also enable \code{include_high_impact}.}
  \item{apoe_col}{Individual-table column containing APOE genotype labels
  matching \code{APOE_genotype} in the APOE reference.}
  \item{rv_reference}{RV effect reference accepted by
  \code{prepare_rv_reference()}. Each row requires unique \code{RV_IDs} and either
  positive \code{OR} or frequencies from which it can be calculated;
  \code{RV_IDs} is matched to carrier status. Common imported-table aliases are accepted.}
  \item{rv_columns}{Optional named character vector mapping RV-reference IDs to 0/1 or logical columns.}
  \item{rv_status_col}{Optional column containing delimited values matching
  \code{RV_IDs} in the RV reference.}
  \item{rv_delimiter}{Regular expression separating RV IDs.}
  \item{rv_prevalence}{Disease-probability prior used to translate RV odds ratios. A value of 0.50 gives a balanced-sample update. This is not RV allele frequency.}
  \item{apoe_reference}{High-impact genotype reference. If omitted, a
  hypothetical APOE example is used with a message.}
  \item{verbose}{Whether \code{check_tiger_inputs()} prints a concise validation summary.}
}
\value{
\code{tiger_probabilities()} returns the input data frame with aligned
probability and optional RV-component columns appended. Reference tables used
by the calculation are retained as attributes. \code{check_tiger_inputs()}
invisibly returns a validation-summary list.
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
The included APOE frequencies are hypothetical and must not be treated as
empirical reference data. Supply a reviewed population- and phenotype-matched
reference for an applied analysis.
}
\examples{
# Begin with PRS only
individuals <- data.frame(
  ID = c("P001", "P002", "P003"),
  PRS_liability = c(-0.5, 0, 0.5),
  RV_IDs = c("", "RISK_A", "PROTECT_C"),
  APOE_genotype = c("e2/e3", "e3/e3", "e3/e4")
)

prs_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1
)

# Add RV and APOE as the high-impact component after PRS-only works.
# All values in this example are hypothetical.
rv_reference <- data.frame(
  RV_IDs = c("RISK_A", "PROTECT_C"),
  Symbol = c("GENE_A", "GENE_C"),
  Class = c("PTV", "CNV"),
  Case_freq = c(0.010, 0.003),
  Control_freq = c(0.002, 0.008),
  OR = c(8.0, 0.4)
)

check_tiger_inputs(
  individuals,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_IDs",
  include_apoe = TRUE, apoe_col = "APOE_genotype",
  apoe_reference = tiger_default_apoe_reference()
)

full_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_IDs", rv_prevalence = 0.5,
  include_apoe = TRUE, apoe_col = "APOE_genotype",
  apoe_reference = tiger_default_apoe_reference()
)

full_results[c(
  "ID", "Probability_PRS", "Probability_PRS_RV",
  "Probability_PRS_APOE", "Probability_PRS_APOE_RV"
)]

# The generic route also accepts a single 0/1/2 high-impact variant
single_variant_reference <- biallelic_genotype_reference(0.20, 0.10)
individuals$High_impact_genotype <- c("0", "1", "2")
single_variant_results <- tiger_probabilities(
  individuals, K = 0.01, SP = 0.5,
  method = "PAIR (summary)", r2_liability = 0.1,
  include_high_impact = TRUE,
  high_impact_col = "High_impact_genotype",
  high_impact_reference = single_variant_reference
)
}
