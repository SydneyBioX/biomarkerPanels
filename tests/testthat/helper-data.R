fixture_path <- function(name) {
  file.path(testthat::test_path("..", "data"), name)
}

# Skip helpers for slow tests
# These allow fine-grained control over which tests run in different contexts

#' Skip tests that are too slow for CRAN/Bioconductor builds
#'
#' Use for tests that run NSGA optimization with small populations
#' (typically 5-30 seconds per test)
skip_slow_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_bioc()
}

#' Skip very intensive tests unless explicitly enabled
#'
#' Use for tests that run NSGA with large populations or many generations
#' (typically > 60 seconds per test). Only run when BIOMARKERPANELS_EXTENDED_TESTS=true
skip_extended_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_bioc()
  if (!identical(Sys.getenv("BIOMARKERPANELS_EXTENDED_TESTS"), "true")) {
    testthat::skip("Extended tests disabled. Set BIOMARKERPANELS_EXTENDED_TESTS=true to run.")
  }
}

