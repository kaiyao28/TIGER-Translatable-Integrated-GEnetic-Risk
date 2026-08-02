# Locate a supplied synthetic example file in either an installed package or a
# repository checkout.
tiger_example_file <- function(filename, must_work = TRUE) {
  .tiger_extdata_file("example", filename, must_work = must_work)
}

# Locate a supplied disorder reference file. These are literature-derived
# inputs requiring user review, not universal defaults.
tiger_reference_file <- function(disorder, filename, must_work = TRUE) {
  disorder <- match.arg(disorder, c("SCZ", "AD"))
  .tiger_extdata_file("reference", disorder, filename,
                      must_work = must_work)
}

.tiger_extdata_file <- function(..., must_work = TRUE) {
  components <- c(...)
  installed <- do.call(
    system.file,
    c(as.list(c("extdata", components)), list(package = "TIGER"))
  )
  path <- if (nzchar(installed)) installed else
    do.call(file.path, as.list(c("inst", "extdata", components)))
  if (isTRUE(must_work) && !file.exists(path)) {
    stop("TIGER supplied file was not found: ", path)
  }
  path
}
