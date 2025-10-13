fixture_path <- function(name) {
  file.path(testthat::test_path("..", "data"), name)
}
