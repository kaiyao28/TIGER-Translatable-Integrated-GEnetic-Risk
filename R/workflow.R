# High-level, merged-table TIGER workflow.
# Copyright (C) 2026 TIGER study authors
# Licensed under GNU GPL v3 or later; distributed WITHOUT ANY WARRANTY.

# Bundled six-genotype APOE reference used by the worked example. Users should
# replace it when population- and phenotype-matched frequencies are available.
tiger_default_apoe_reference <- function() {
  prepare_high_impact_reference(data.frame(
    Genotype = c("e2/e2", "e2/e3", "e2/e4", "e3/e3", "e3/e4", "e4/e4"),
    Case_freq = c(0.0016, 0.0528, 0.0240, 0.4356, 0.3960, 0.0900),
    Control_freq = c(0.0064, 0.1264, 0.0208, 0.6241, 0.2054, 0.0169)
  ))
}

.tiger_column <- function(data, column, argument, required = TRUE) {
  if (is.null(column)) {
    if (required) stop(argument, " must name one column")
    return(NULL)
  }
  if (length(column) != 1L || is.na(column) || !nzchar(column) ||
      !column %in% names(data)) {
    stop(argument, " must name one column in data")
  }
  column
}

.tiger_prs_probability <- function(
    prs, method, K, SP, r2_liability, reference_prs_liability,
    r2_observed, case_mean, case_sd, control_mean, control_sd, n_quantiles) {
  switch(
    method,
    "BPC" = bpc_probability(prs, K, SP, r2_liability),
    "GenoPred" = {
      if (is.null(r2_observed) && !is.null(r2_liability)) {
        r2_observed <- liability_to_observed_r2(r2_liability, K, SP)
      }
      genopred_probability(
        prs, reference_prs_liability, r2_observed, K, SP,
        n_quantiles = n_quantiles
      )
    },
    "PAIR (summary)" = pair_probability_summary(
      prs, K, SP, r2_liability
    ),
    "PAIR (sample)" = pair_probability_sample(
      prs, case_mean, case_sd, control_mean, control_sd, K, SP
    )
  )
}

.tiger_merged_carrier_matrix <- function(
    data, ids, reference, rv_columns = NULL, rv_status_col = NULL,
    rv_delimiter = "[;,|]") {
  if (!is.null(rv_columns) && !is.null(rv_status_col)) {
    stop("supply rv_columns or rv_status_col, not both")
  }
  variant_ids <- reference$ID
  out <- matrix(
    FALSE, nrow = nrow(data), ncol = length(variant_ids),
    dimnames = list(ids, variant_ids)
  )
  if (!is.null(rv_columns)) {
    if (!is.character(rv_columns) || !length(rv_columns)) {
      stop("rv_columns must contain merged-table column names")
    }
    if (is.null(names(rv_columns)) || any(!nzchar(names(rv_columns)))) {
      names(rv_columns) <- rv_columns
    }
    if (any(!rv_columns %in% names(data))) {
      stop("unknown rv_columns: ",
           paste(setdiff(rv_columns, names(data)), collapse = ", "))
    }
    if (any(!names(rv_columns) %in% variant_ids)) {
      stop("rv_columns names must match RV-reference IDs")
    }
    for (variant in names(rv_columns)) {
      value <- data[[rv_columns[[variant]]]]
      if (is.numeric(value) && all(value %in% c(0, 1))) {
        value <- as.logical(value)
      }
      if (!is.logical(value) || anyNA(value)) {
        stop("RV columns must contain logical or 0/1 values")
      }
      out[, variant] <- value
    }
    return(out)
  }
  rv_status_col <- .tiger_column(
    data, rv_status_col, "rv_status_col", required = TRUE
  )
  status <- data[[rv_status_col]]
  if ((is.logical(status) || is.numeric(status)) && length(variant_ids) == 1L) {
    if (is.numeric(status) && all(status %in% c(0, 1), na.rm = TRUE)) {
      status <- as.logical(status)
    }
    if (!is.logical(status)) stop("single-RV status must be logical or 0/1")
    status[is.na(status)] <- FALSE
    out[, 1] <- status
    return(out)
  }
  status <- as.character(status)
  status[is.na(status)] <- ""
  for (i in seq_along(status)) {
    if (!nzchar(trimws(status[i]))) next
    carried <- trimws(strsplit(status[i], rv_delimiter)[[1]])
    carried <- unique(carried[nzchar(carried)])
    unknown <- setdiff(carried, variant_ids)
    if (length(unknown)) {
      stop("unknown RV ID(s) in row ", i, ": ",
           paste(unknown, collapse = ", "))
    }
    out[i, carried] <- TRUE
  }
  out
}

# Calculate aligned TIGER probabilities from one merged individual table.
tiger_probabilities <- function(
    data, K, SP = 0.5,
    method = c("PAIR (summary)", "BPC", "GenoPred", "PAIR (sample)"),
    id_col = "ID", prs_col = "PRS_liability", group_col = NULL,
    r2_liability = NULL, reference_prs_liability = NULL,
    r2_observed = NULL, case_mean = NULL, case_sd = NULL,
    control_mean = NULL, control_sd = NULL, n_quantiles = 100,
    include_rv = FALSE, rv_reference = NULL,
    rv_columns = NULL, rv_status_col = NULL,
    rv_delimiter = "[;,|]", rv_prevalence = 0.5,
    include_apoe = FALSE, apoe_col = "APOE", apoe_reference = NULL) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a non-empty data.frame")
  }
  method <- match.arg(method)
  if (!is.logical(include_rv) || length(include_rv) != 1L || is.na(include_rv)) {
    stop("include_rv must be TRUE or FALSE")
  }
  if (!is.logical(include_apoe) || length(include_apoe) != 1L ||
      is.na(include_apoe)) {
    stop("include_apoe must be TRUE or FALSE")
  }
  id_col <- .tiger_column(data, id_col, "id_col")
  prs_col <- .tiger_column(data, prs_col, "prs_col")
  if (!is.null(group_col)) .tiger_column(data, group_col, "group_col")
  ids <- as.character(data[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("ID values must be non-missing and unique")
  }
  prs <- data[[prs_col]]
  if (!is.numeric(prs) || any(!is.finite(prs))) {
    stop("the PRS column must contain finite numeric liability-scale values")
  }
  .check_probability(K, "K")
  .check_probability(SP, "SP")
  p_prs <- .tiger_prs_probability(
    prs, method, K, SP, r2_liability, reference_prs_liability,
    r2_observed, case_mean, case_sd, control_mean, control_sd, n_quantiles
  )
  output <- data
  output$Probability_PRS <- p_prs

  carrier_matrix <- NULL
  odds_ratios <- NULL
  if (isTRUE(include_rv)) {
    if (is.null(rv_reference)) {
      stop("rv_reference is required when include_rv = TRUE")
    }
    rv_reference <- prepare_rv_reference(rv_reference)
    carrier_matrix <- .tiger_merged_carrier_matrix(
      data, ids, rv_reference, rv_columns, rv_status_col, rv_delimiter
    )
    odds_ratios <- stats::setNames(rv_reference$OR, rv_reference$ID)
    rv <- apply_rv_carriers(
      stats::setNames(p_prs, ids), carrier_matrix, odds_ratios,
      prevalence = rv_prevalence
    )
    output$RV_count <- rv$RV_count
    output$RV_damaging_probability <- rv$p_damaging
    output$RV_protective_probability <- rv$p_protective
    output$Probability_PRS_RV <- rv$probability_after
  }

  if (isTRUE(include_apoe)) {
    apoe_col <- .tiger_column(data, apoe_col, "apoe_col")
    if (is.null(apoe_reference)) {
      apoe_reference <- tiger_default_apoe_reference()
      message(
        "Using TIGER's bundled APOE reference. Supply apoe_reference to ",
        "use population- and phenotype-matched frequencies."
      )
    } else {
      apoe_reference <- prepare_high_impact_reference(apoe_reference)
    }
    p_apoe <- high_impact_method_probability(
      prs, data[[apoe_col]], apoe_reference, K = K, SP = SP,
      method = method, r2_liability = r2_liability,
      reference_prs_liability = reference_prs_liability,
      r2_observed = r2_observed, case_mean = case_mean, case_sd = case_sd,
      control_mean = control_mean, control_sd = control_sd,
      n_quantiles = n_quantiles
    )
    output$Probability_PRS_APOE <- p_apoe
    if (isTRUE(include_rv)) {
      combined <- apply_rv_carriers(
        stats::setNames(p_apoe, ids), carrier_matrix, odds_ratios,
        prevalence = rv_prevalence
      )
      output$Probability_PRS_APOE_RV <- combined$probability_after
    }
  }
  attr(output, "TIGER_method") <- method
  attr(output, "TIGER_RV_reference") <- if (isTRUE(include_rv)) {
    rv_reference
  } else NULL
  attr(output, "TIGER_APOE_reference") <- if (isTRUE(include_apoe)) {
    apoe_reference
  } else NULL
  output
}
