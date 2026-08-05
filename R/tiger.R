# Load all TIGER functions from a repository checkout.

load_tiger <- function(root = ".", envir = parent.frame()) {
  files <- file.path(root, "R", c(
    "files.R",
    "probability_methods.R",
    "high_impact_variants.R",
    "rare_variant_probability.R",
    "reference_data.R",
    "workflow.R",
    "plotting.R"
  ))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop(
      "TIGER files were not found. Set root to the TIGER repository folder. ",
      "Missing: ", paste(missing, collapse = ", ")
    )
  }
  for (file in files) sys.source(file, envir = envir)
  invisible(files)
}
