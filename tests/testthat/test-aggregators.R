# Tests for cohort aggregator registry and implementations

# Access internal functions via :::
.get_aggregator <- biomarkerPanels:::.get_aggregator
.apply_cohort_aggregator <- biomarkerPanels:::.apply_cohort_aggregator
.aggregator_registry <- biomarkerPanels:::.aggregator_registry
aggregator_none <- biomarkerPanels:::aggregator_none
aggregator_pairwise_ratios <- biomarkerPanels:::aggregator_pairwise_ratios
aggregator_pairwise_log_ratios <- biomarkerPanels:::aggregator_pairwise_log_ratios
aggregator_reference_norm <- biomarkerPanels:::aggregator_reference_norm

test_that("aggregator registry stores default aggregators", {
  reg <- aggregator_registry()
  expect_type(reg, "list")
  expect_true("none" %in% names(reg))
  expect_true("pairwise_ratios" %in% names(reg))
  expect_true("pairwise_log_ratios" %in% names(reg))
  expect_true("reference_norm" %in% names(reg))
})

test_that("register_aggregator adds new aggregators", {
  custom_agg <- function(x) x * 2
  register_aggregator("test_doubler", custom_agg, "Doubles all values", overwrite = TRUE)

  reg <- aggregator_registry()
  expect_true("test_doubler" %in% names(reg))
  expect_identical(reg$test_doubler$fun, custom_agg)
  expect_equal(reg$test_doubler$description, "Doubles all values")

  # Clean up
  rm("test_doubler", envir = .aggregator_registry)
})

test_that("register_aggregator prevents overwriting without overwrite=TRUE", {
  expect_error(
    register_aggregator("none", function(x) x, overwrite = FALSE),
    "already registered"
  )
})

test_that("register_aggregator validates inputs", {
  expect_error(register_aggregator(123, function(x) x), "is.character")
  expect_error(register_aggregator("", function(x) x), "nzchar")
  expect_error(register_aggregator("test", "not a function"), "is.function")
})

test_that(".get_aggregator retrieves registered aggregators", {
  agg <- .get_aggregator("none")
  expect_type(agg, "list")
  expect_true("fun" %in% names(agg))
  expect_true("description" %in% names(agg))
})

test_that(".get_aggregator errors on unknown aggregator", {
  expect_error(
    .get_aggregator("nonexistent_aggregator"),
    "is not registered"
  )
})

test_that("aggregator_none returns input unchanged", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- paste0("F", 1:4)
  result <- aggregator_none(x)
  expect_identical(result, x)
})

test_that("aggregator_pairwise_ratios calls pairwise_col_diff", {
  x <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, ncol = 2)
  colnames(x) <- c("A", "B")
  result <- aggregator_pairwise_ratios(x)

  expect_equal(ncol(result), 1)
  expect_true(grepl("--", colnames(result)[1]))
  expect_equal(result[, 1], x[, "A"] - x[, "B"])
})

test_that("aggregator_pairwise_ratios generates all pairs", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  result <- aggregator_pairwise_ratios(x)

  expected_pairs <- 4 * 3 / 2
  expect_equal(ncol(result), expected_pairs)
  expect_true(all(grepl("--", colnames(result))))
})

test_that("aggregator_pairwise_ratios errors on single feature", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "A"

  expect_error(
    aggregator_pairwise_ratios(x),
    "requires at least two features"
  )
})

test_that("aggregator_pairwise_log_ratios computes log ratios", {
  x <- matrix(c(10, 20, 5, 10), nrow = 2, ncol = 2)
  colnames(x) <- c("A", "B")
  result <- aggregator_pairwise_log_ratios(x)

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result)[1], "A--B")
  expect_equal(result[, 1], log(x[, "A"] / x[, "B"]))
})

test_that("aggregator_pairwise_log_ratios generates all pairs", {
  x <- matrix(abs(rnorm(12)) + 1, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  result <- aggregator_pairwise_log_ratios(x)

  expected_pairs <- 4 * 3 / 2
  expect_equal(ncol(result), expected_pairs)
  expect_true(all(grepl("--", colnames(result))))
})

test_that("aggregator_pairwise_log_ratios errors on single feature", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "A"

  expect_error(
    aggregator_pairwise_log_ratios(x),
    "requires at least two features"
  )
})

test_that("aggregator_pairwise_log_ratios requires column names", {
  x <- matrix(1:6, nrow = 3, ncol = 2)

  expect_error(
    aggregator_pairwise_log_ratios(x),
    "must have column names"
  )
})

test_that("aggregator_reference_norm requires reference_feature attribute", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")

  expect_error(
    aggregator_reference_norm(x),
    "requires 'reference_feature' attribute"
  )
})

test_that("aggregator_reference_norm normalizes to reference", {
  x <- matrix(c(10, 20, 30, 5, 10, 15, 2, 4, 6), nrow = 3, ncol = 3)
  colnames(x) <- c("A", "B", "Ref")
  attr(x, "reference_feature") <- "Ref"

  result <- aggregator_reference_norm(x)

  expect_equal(ncol(result), 2)
  expect_true(all(c("A--Ref", "B--Ref") %in% colnames(result)))
  expect_equal(result[, "A--Ref"], x[, "A"] - x[, "Ref"])
  expect_equal(result[, "B--Ref"], x[, "B"] - x[, "Ref"])
})

test_that("aggregator_reference_norm validates reference feature exists", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  attr(x, "reference_feature") <- "NonExistent"

  expect_error(
    aggregator_reference_norm(x),
    "not found in data"
  )
})

test_that("aggregator_reference_norm warns when only reference feature exists", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "Ref"
  attr(x, "reference_feature") <- "Ref"

  expect_warning(
    result <- aggregator_reference_norm(x),
    "requires at least one non-reference feature"
  )
  expect_identical(result, x)
})

test_that(".apply_cohort_aggregator applies aggregator to all matrices", {
  mat1 <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(mat1) <- c("A", "B")
  mat2 <- matrix(7:12, nrow = 3, ncol = 2)
  colnames(mat2) <- c("A", "B")

  result <- .apply_cohort_aggregator(list(mat1, mat2), "pairwise_ratios")

  expect_length(result, 2)
  expect_equal(ncol(result[[1]]), 1)
  expect_equal(ncol(result[[2]]), 1)
})

test_that(".apply_cohort_aggregator handles empty list", {
  result <- .apply_cohort_aggregator(list(), "none")
  expect_length(result, 0)
})

test_that(".apply_cohort_aggregator works with none aggregator", {
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(mat) <- c("A", "B")

  result <- .apply_cohort_aggregator(list(mat), "none")

  expect_length(result, 1)
  expect_identical(result[[1]], mat)
})
