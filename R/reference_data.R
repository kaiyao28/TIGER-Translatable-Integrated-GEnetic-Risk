# TIGER disorder-specific RV reference harmonisation.
# Copyright (C) 2026 TIGER study authors
# Licensed under GNU GPL v3 or later; distributed WITHOUT ANY WARRANTY.

# Harmonise one raw SCZ reference table to the TIGER RV schema.
harmonise_scz_reference <- function(x, class = c("PTV", "MPC2", "CNV")) {
  class <- match.arg(class)
  if (!is.data.frame(x)) stop("x must be a data.frame")
  required <- switch(
    class,
    PTV = c("ID", "Symbol", "CAP_freq", "COP_freq", "CP_OR"),
    MPC2 = c("ID", "Symbol", "CAM2_freq", "COM2_freq", "CM2_OR"),
    CNV = c("Syndrome", "Case_freq", "Control_freq")
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing ", class, " column(s): ",
                            paste(missing, collapse = ", "))
  if (class == "PTV") {
    out <- data.frame(
      ID = x$ID, Symbol = x$Symbol, Class = class,
      Case_freq = x$CAP_freq, Control_freq = x$COP_freq, OR = x$CP_OR
    )
  } else if (class == "MPC2") {
    out <- data.frame(
      ID = x$ID, Symbol = x$Symbol, Class = class,
      Case_freq = x$CAM2_freq, Control_freq = x$COM2_freq, OR = x$CM2_OR
    )
  } else {
    odds_ratio <- (x$Case_freq / (1 - x$Case_freq)) /
      (x$Control_freq / (1 - x$Control_freq))
    out <- data.frame(
      ID = x$Syndrome, Symbol = x$Syndrome, Class = class,
      Case_freq = x$Case_freq, Control_freq = x$Control_freq, OR = odds_ratio
    )
  }
  prepare_rv_reference(out, source = paste0("SCZ_", class))
}

# Harmonise one raw AD reference table. LOF and REVEL classes use gene-level
# Case_OR; SRV uses published variant IDs and OR values directly.
harmonise_ad_reference <- function(
    x, class = c("LOF", "REVEL25", "REVEL50", "REVEL75", "SRV")) {
  class <- match.arg(class)
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (class == "SRV") {
    required <- c("rsID", "Gene", "Case_freq", "Control_freq", "OR")
    missing <- setdiff(required, names(x))
    if (length(missing)) stop("Missing SRV column(s): ",
                              paste(missing, collapse = ", "))
    out <- data.frame(
      ID = x$rsID,
      Symbol = paste(x$Gene, x$rsID, sep = ":"),
      Class = class,
      Case_freq = x$Case_freq,
      Control_freq = x$Control_freq,
      OR = x$OR
    )
  } else {
    required <- c("Gene", "Case_freq", "Control_freq", "Case_OR")
    missing <- setdiff(required, names(x))
    if (length(missing)) stop("Missing ", class, " column(s): ",
                              paste(missing, collapse = ", "))
    out <- data.frame(
      ID = paste(class, x$Gene, sep = ":"),
      Symbol = x$Gene,
      Class = class,
      Case_freq = x$Case_freq,
      Control_freq = x$Control_freq,
      OR = x$Case_OR
    )
  }
  prepare_rv_reference(out, source = paste0("AD_", class))
}

# Read and harmonise one supplied or user-created Excel reference workbook.
read_tiger_reference <- function(path, disorder = c("SCZ", "AD"), class) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("read_tiger_reference() requires the readxl package")
  }
  if (length(path) != 1L || !file.exists(path)) stop("reference file not found")
  disorder <- match.arg(disorder)
  raw <- readxl::read_excel(path)
  if (disorder == "SCZ") {
    harmonise_scz_reference(raw, class)
  } else {
    harmonise_ad_reference(raw, class)
  }
}
