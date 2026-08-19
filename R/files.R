# Locate a supplied synthetic example file in either an installed package or a
# repository checkout.
tiger_example_file <- function(filename, must_work = TRUE) {
  .tiger_extdata_file("example", filename, must_work = must_work)
}

# Locate a blank or hypothetical CSV input template.
tiger_template_file <- function(filename, must_work = TRUE) {
  .tiger_extdata_file("templates", filename, must_work = must_work)
}

.tiger_extdata_file <- function(..., must_work = TRUE) {
  components <- c(...)
  local <- do.call(file.path, as.list(c("inst", "extdata", components)))
  installed <- do.call(
    system.file,
    c(as.list(c("extdata", components)), list(package = "TIGER"))
  )
  # A repository checkout should use its own files, even when an older TIGER
  # version is installed in the active R library.
  path <- if (file.exists(local)) local else installed
  if (isTRUE(must_work) && !file.exists(path)) {
    stop("TIGER data file was not found: ", path)
  }
  path
}
