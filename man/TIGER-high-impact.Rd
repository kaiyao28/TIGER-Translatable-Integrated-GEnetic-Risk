\name{TIGER-high-impact}
\alias{prepare_high_impact_reference}
\alias{apoe_genotype_reference}
\alias{biallelic_genotype_reference}
\alias{apply_high_impact_probability}
\alias{high_impact_prevalence_parameters}
\alias{high_impact_method_probability}
\title{Prepare and apply common high-impact genotype references}
\description{
Prepare mutually exclusive case/control genotype frequencies and use them to
update liability-scale PRS probabilities. The same interface supports APOE's
six genotypes or a single biallelic variant coded as 0, 1 and 2 effect alleles.
}
\usage{
prepare_high_impact_reference(x, tolerance = 1e-8)

apoe_genotype_reference(case_allele_frequencies,
                        control_allele_frequencies)

biallelic_genotype_reference(case_effect_allele_frequency,
                             control_effect_allele_frequency,
                             genotype_labels = c("0", "1", "2"))

apply_high_impact_probability(p_prs, genotype, reference)

high_impact_prevalence_parameters(reference, K, SP = 0.5)

high_impact_method_probability(
  prs_liability, genotype, reference, K, SP = 0.5,
  method = c("PAIR (summary)", "PAIR (sample)", "BPC", "GenoPred"),
  r2_liability = NULL, reference_prs_liability = NULL,
  r2_observed = NULL, case_mean = NULL, case_sd = NULL,
  control_mean = NULL, control_sd = NULL, n_quantiles = 100
)
}
\arguments{
  \item{x}{Data frame with unique \code{Genotype}, \code{Case_freq} and \code{Control_freq} columns. Case and control frequencies must each sum to one.}
  \item{tolerance}{Numerical tolerance used when checking frequency sums.}
  \item{case_allele_frequencies, control_allele_frequencies}{Named e2, e3 and e4 allele-frequency vectors for cases and controls.}
  \item{case_effect_allele_frequency, control_effect_allele_frequency}{Effect-allele frequencies in cases and controls for one biallelic variant.}
  \item{genotype_labels}{Three unique labels corresponding to 0, 1 and 2 effect alleles.}
  \item{p_prs}{PRS-derived probabilities in the required target scale.}
  \item{prs_liability}{Individual PRSs already transformed to the liability scale.}
  \item{genotype}{Observed genotype label for each individual. Labels must match the reference.}
  \item{reference}{A high-impact genotype reference accepted by \code{prepare_high_impact_reference()}.}
  \item{K}{Population disorder prevalence.}
  \item{SP}{Target-sample prevalence.}
  \item{method}{PRS-to-probability conversion method.}
  \item{r2_liability, reference_prs_liability, r2_observed}{Inputs required by the selected summary-based method.}
  \item{case_mean, case_sd, control_mean, control_sd}{Case and control PRS moments required by PAIR (sample).}
  \item{n_quantiles}{Number of probability quantiles used by GenoPred.}
}
\value{
\code{prepare_high_impact_reference()}, \code{apoe_genotype_reference()} and
\code{biallelic_genotype_reference()} return validated genotype-reference data
frames. \code{high_impact_prevalence_parameters()} appends genotype-specific
\code{K_genotype} and \code{SP_genotype}. The probability functions return
numeric probability vectors.
}
\details{
Allele-frequency helpers derive genotype frequencies under Hardy-Weinberg
equilibrium within cases and controls. Prefer observed case/control genotype
frequencies when suitable estimates are available.

\code{high_impact_method_probability()} is the canonical TIGER update. It
derives genotype-specific population and sample prevalences, then recalculates
the selected PRS conversion for each genotype. The lower-level
\code{apply_high_impact_probability()} instead applies an equivalent
case/control genotype likelihood-ratio update to existing probabilities.

The separately modelled variant, or an appropriate LD-aware region around it,
must be excluded from PRS construction and from the PRS R-squared estimate to
avoid double counting. References should match the phenotype, ancestry and
target population as closely as possible.
}
\examples{
# One biallelic high-impact variant under HWE in cases and controls
reference <- biallelic_genotype_reference(
  case_effect_allele_frequency = 0.20,
  control_effect_allele_frequency = 0.10
)

high_impact_method_probability(
  prs_liability = c(-0.5, 0, 0.5),
  genotype = c("0", "1", "2"),
  reference = reference,
  K = 0.01, SP = 0.5,
  method = "PAIR (summary)",
  r2_liability = 0.1
)

# APOE: named e2/e3/e4 allele frequencies are converted to six genotypes
apoe_reference <- apoe_genotype_reference(
  case_allele_frequencies = c(e2 = 0.04, e3 = 0.60, e4 = 0.36),
  control_allele_frequencies = c(e2 = 0.08, e3 = 0.78, e4 = 0.14)
)
}
