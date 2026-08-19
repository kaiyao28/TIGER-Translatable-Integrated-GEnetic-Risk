# PRS and rare-variant probability toolkit.
# Copyright (C) 2026 TIGER study authors
# Licensed under GNU GPL v3 or later; distributed WITHOUT ANY WARRANTY.
# See LICENSE.

# Translate an odds ratio into an intrinsic damaging or protective probability.
# For OR >= 1, p = prevalence*(OR-1)/(1 + prevalence*(OR-1)). For OR < 1,
# TIGER applies the inverted effect to the complementary, disease-free
# probability, which simplifies to
# (1-prevalence)*(1-OR)/(1-prevalence*(1-OR)).
intrinsic_rv_probability <- function(odds_ratio, prevalence = 0.5) {
  if (length(prevalence) != 1L || !is.finite(prevalence) ||
      prevalence <= 0 || prevalence >= 1) {
    stop("prevalence must be one finite value strictly between 0 and 1")
  }
  if (!is.numeric(odds_ratio) || any(!is.finite(odds_ratio)) ||
      any(odds_ratio <= 0)) stop("odds_ratio must contain positive finite values")
  damaging <- odds_ratio >= 1
  result <- numeric(length(odds_ratio))
  result[damaging] <- prevalence * (odds_ratio[damaging] - 1) /
    (1 + prevalence * (odds_ratio[damaging] - 1))
  result[!damaging] <- (1 - prevalence) * (1 - odds_ratio[!damaging]) /
    (1 - prevalence * (1 - odds_ratio[!damaging]))
  pmin(pmax(result, 0), 1)
}

# Combine independent intrinsic probabilities of the same direction.
combine_independent_rv_probabilities <- function(probabilities) {
  if (!length(probabilities)) return(0)
  if (any(!is.finite(probabilities)) || any(probabilities < 0) ||
      any(probabilities > 1)) stop("probabilities must be in [0, 1]")
  1 - prod(1 - probabilities)
}

# Add independent damaging and protective RV effects to any PRS probability:
# (1-p_protective) * {1-(1-p_prs)(1-p_damaging)}.
apply_rv_probability <- function(p_prs, p_damaging = 0,
                                 p_protective = 0) {
  inputs <- c(p_prs, p_damaging, p_protective)
  if (any(!is.finite(inputs)) || any(inputs < 0) || any(inputs > 1)) {
    stop("all probability inputs must be finite and in [0, 1]")
  }
  result <- (1 - p_protective) *
    (1 - (1 - p_prs) * (1 - p_damaging))
  pmin(pmax(result, 0), 1)
}

# Validate the reusable reference schema. Each row represents one RV. ID plus
# either OR or frequencies are required. Missing control frequency may be
# replaced by a supplied population frequency before OR calculation.
prepare_rv_reference <- function(x, source = "RV") {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  normalise <- function(z) tolower(gsub("[^[:alnum:]]", "", z))
  keys <- normalise(names(x))
  locate <- function(aliases, required = TRUE) {
    hit <- match(normalise(aliases), keys, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) return(names(x)[hit[1]])
    if (required) stop("Missing column: ", aliases[1])
    NULL
  }
  id <- locate(c("RV_IDs", "RV_ID", "ID", "Variant_ID"))
  symbol <- locate(c("Symbol", "Label", "Variant"), FALSE)
  case_frequency <- locate(c("Case_freq", "Case_frequency"), FALSE)
  control_frequency <- locate(c("Control_freq", "Control_frequency"), FALSE)
  population_frequency <- locate(
    c("Population_freq", "Population_frequency", "Pop_freq"), FALSE
  )
  odds_ratio <- locate(c("OR", "Odds_ratio"), FALSE)
  if (is.null(odds_ratio) &&
      (is.null(case_frequency) ||
       (is.null(control_frequency) && is.null(population_frequency)))) {
    stop(
      "Supply OR, or supply Case_freq with Control_freq or Population_freq"
    )
  }
  class <- locate(c("Class", "Variant_class", "Type"), FALSE)
  ids <- x[[id]]
  symbols <- if (is.null(symbol)) ids else x[[symbol]]
  numeric_column <- function(column) {
    if (is.null(column)) rep(NA_real_, nrow(x)) else as.numeric(x[[column]])
  }
  case_values <- numeric_column(case_frequency)
  control_values <- numeric_column(control_frequency)
  population_values <- numeric_column(population_frequency)
  control_replaced <- !is.finite(control_values) & is.finite(population_values)
  control_values[control_replaced] <- population_values[control_replaced]
  or_values <- numeric_column(odds_ratio)
  or_source <- ifelse(is.finite(or_values) & or_values > 0, "Supplied", NA_character_)
  derive_or <- is.na(or_source)
  derivable <- derive_or & is.finite(case_values) & is.finite(control_values) &
    case_values > 0 & case_values < 1 &
    control_values > 0 & control_values < 1
  or_values[derivable] <-
    (case_values[derivable] / (1 - case_values[derivable])) /
    (control_values[derivable] / (1 - control_values[derivable]))
  or_source[derivable] <- "Calculated from frequencies"
  output <- data.frame(
    ID = as.character(ids),
    Symbol = as.character(symbols),
    Class = if (is.null(class)) source else as.character(x[[class]]),
    Case_freq = case_values,
    Control_freq = control_values,
    Population_freq = population_values,
    OR = or_values,
    OR_source = or_source,
    Control_freq_source = ifelse(
      control_replaced, "Population frequency",
      ifelse(is.finite(control_values), "Control frequency", NA_character_)
    ),
    stringsAsFactors = FALSE
  )
  valid <- !is.na(output$ID) & nzchar(output$ID) &
    !is.na(output$Symbol) & nzchar(output$Symbol) &
    !is.na(output$Class) & nzchar(output$Class) &
    is.finite(output$OR) & output$OR > 0
  supplied_frequency_values <- cbind(
    output$Case_freq, output$Control_freq, output$Population_freq
  )
  invalid_frequency <- apply(supplied_frequency_values, 1, function(z) {
    z <- z[!is.na(z)]
    any(!is.finite(z) | z < 0 | z > 1)
  })
  valid <- valid & !invalid_frequency
  if (any(!valid)) warning(sum(!valid), " invalid RV rows removed")
  output <- output[valid, , drop = FALSE]
  if (!nrow(output)) stop("No valid RV rows remain")
  if (any(!nzchar(output$ID)) || anyDuplicated(output$ID)) {
    stop("RV reference IDs must be non-empty and unique")
  }
  output$Direction <- ifelse(output$OR >= 1, "Damaging", "Protective")
  rownames(output) <- NULL
  output
}

# Convert a long individual-by-variant carrier table into a logical matrix.
# Required columns are ID, Variant_ID and Carrier. Missing ID/variant pairs in
# the table are treated as non-carriers only when a complete individual_ids and
# variant_ids set is supplied explicitly.
prepare_rv_carrier_matrix <- function(x, individual_ids = NULL,
                                      variant_ids = NULL) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  required <- c("ID", "Variant_ID", "Carrier")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing column(s): ", paste(missing, collapse = ", "))
  ids <- as.character(x$ID)
  variants <- as.character(x$Variant_ID)
  carrier <- x$Carrier
  if (is.numeric(carrier) && all(carrier %in% c(0, 1))) carrier <- as.logical(carrier)
  if (!is.logical(carrier) || anyNA(carrier) || any(!nzchar(ids)) ||
      any(!nzchar(variants))) {
    stop("ID/Variant_ID must be non-empty and Carrier must be logical or 0/1")
  }
  key <- paste(ids, variants, sep = "\r")
  if (anyDuplicated(key)) stop("carrier table contains duplicate ID/Variant_ID rows")
  sparse_table <- nrow(x) < length(unique(ids)) * length(unique(variants))
  if (sparse_table && (is.null(individual_ids) || is.null(variant_ids))) {
    stop("for a presence-only carrier table, supply complete individual_ids and variant_ids")
  }
  if (is.null(individual_ids)) individual_ids <- unique(ids)
  if (is.null(variant_ids)) variant_ids <- unique(variants)
  individual_ids <- as.character(individual_ids)
  variant_ids <- as.character(variant_ids)
  if (!length(individual_ids) || !length(variant_ids) ||
      anyNA(individual_ids) || anyNA(variant_ids) ||
      any(!nzchar(individual_ids)) || any(!nzchar(variant_ids)) ||
      anyDuplicated(individual_ids) || anyDuplicated(variant_ids) ||
      any(!ids %in% individual_ids) || any(!variants %in% variant_ids)) {
    stop(
      "individual_ids and variant_ids must be non-empty, unique, and contain ",
      "all table values"
    )
  }
  out <- matrix(
    FALSE, nrow = length(individual_ids), ncol = length(variant_ids),
    dimnames = list(individual_ids, variant_ids)
  )
  out[cbind(match(ids, individual_ids), match(variants, variant_ids))] <- carrier
  out
}

# Apply zero, one or multiple independently acting RVs per person. carrier_matrix
# has individuals in rows and variants in columns. odds_ratios must be named to
# match the columns (preferred), or supplied in the same order.
apply_rv_carriers <- function(p_prs, carrier_matrix, odds_ratios,
                              prevalence = 0.5) {
  .check_rv_probability <- function(x, name) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(x < 0) || any(x > 1)) stop(name, " must contain values in [0, 1]")
  }
  .check_rv_probability(p_prs, "p_prs")
  if (!is.matrix(carrier_matrix) || !is.logical(carrier_matrix) ||
      anyNA(carrier_matrix) || nrow(carrier_matrix) != length(p_prs) ||
      is.null(colnames(carrier_matrix)) ||
      anyNA(colnames(carrier_matrix)) || any(!nzchar(colnames(carrier_matrix))) ||
      anyDuplicated(colnames(carrier_matrix))) {
    stop("carrier_matrix must be a named-column logical matrix with one row per person")
  }
  probability_ids <- names(p_prs)
  carrier_ids <- rownames(carrier_matrix)
  if (!is.null(probability_ids) && !is.null(carrier_ids) &&
      !identical(probability_ids, carrier_ids)) {
    stop("named p_prs IDs must exactly match carrier_matrix row names and order")
  }
  if (!is.numeric(odds_ratios) || length(odds_ratios) != ncol(carrier_matrix) ||
      any(!is.finite(odds_ratios)) || any(odds_ratios <= 0)) {
    stop("odds_ratios must contain one positive finite value per variant")
  }
  if (!is.null(names(odds_ratios))) {
    if (anyNA(names(odds_ratios)) || any(!nzchar(names(odds_ratios))) ||
        anyDuplicated(names(odds_ratios))) {
      stop("named odds_ratios must have non-empty unique names")
    }
    if (!all(colnames(carrier_matrix) %in% names(odds_ratios))) {
      stop("named odds_ratios must include every carrier-matrix variant")
    }
    odds_ratios <- odds_ratios[colnames(carrier_matrix)]
  }
  intrinsic <- intrinsic_rv_probability(odds_ratios, prevalence)
  damaging <- odds_ratios >= 1
  combine_rows <- function(columns) {
    if (!any(columns)) return(rep(0, nrow(carrier_matrix)))
    selected <- carrier_matrix[, columns, drop = FALSE]
    selected_probabilities <- intrinsic[columns]
    apply(selected, 1, function(z) {
      combine_independent_rv_probabilities(selected_probabilities[z])
    })
  }
  p_damaging <- combine_rows(damaging)
  p_protective <- combine_rows(!damaging)
  data.frame(
    RV_count = rowSums(carrier_matrix),
    Damaging_RV_count = rowSums(carrier_matrix[, damaging, drop = FALSE]),
    Protective_RV_count = rowSums(carrier_matrix[, !damaging, drop = FALSE]),
    p_damaging = p_damaging,
    p_protective = p_protective,
    probability_before = p_prs,
    probability_after = apply_rv_probability(
      p_prs, p_damaging = p_damaging, p_protective = p_protective
    ),
    row.names = rownames(carrier_matrix)
  )
}
