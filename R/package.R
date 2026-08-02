#' TIGER: Translatable Integrated Genetic Risk framework
#'
#' TIGER provides transparent R functions for combining liability-scale
#' polygenic risk, separately modelled common high-impact variants, and rare
#' variant carrier effects. Start with the repository README and
#' `docs/GETTING_STARTED.md` guide.
#'
#' @name TIGER-package
#' @aliases TIGER
#' @docType package
NULL

utils::globalVariables(c(
  "APOE", "Carrier_group", "Condition", "Genotype", "Probability", "Probability_after",
  "Probability_after_RV", "Probability_before", "PRS", "RV_count_group",
  "RV_group"
))
