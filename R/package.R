#' TIGER: Translatable Integrated Genetic Risk
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
  "APOE", "Carrier_group", "Condition", "Display", "Effect", "Genotype",
  "Group", "Intrinsic_probability", "Label_hjust", "Label_y", "Method", "Probability",
  "Probability_after",
  "Probability_after_RV", "Probability_before", "PRS", "PRS_group",
  "RV_count_group", "RV_group"
))
