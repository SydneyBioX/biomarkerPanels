# Tests for feature transform registry and implementations

# Access internal functions via :::
.get_feature_transform <- biomarkerPanels:::.get_feature_transform
.apply_feature_transform <- biomarkerPanels:::.apply_feature_transform
.apply_feature_transform_single <- biomarkerPanels:::.apply_feature_transform_single
.transform_registry <- biomarkerPanels:::.transform_registry
transform_none <- biomarkerPanels:::transform_none
transform_pairwise_ratios <- biomarkerPanels:::transform_pairwise_ratios
transform_pairwise_log_ratios <- biomarkerPanels:::transform_pairwise_log_ratios
transform_reference_norm <- biomarkerPanels:::transform_reference_norm

test_that("feature transform registry stores default transforms", {
  reg <- feature_transform_registry()
  expect_type(reg, "list")
  expect_true("none" %in% names(reg))
  expect_true("pairwise_ratios" %in% names(reg))
  expect_true("pairwise_log_ratios" %in% names(reg))
  expect_true("reference_norm" %in% names(reg))
})

test_that("register_feature_transform adds new transforms", {
  custom_transform <- function(x) x * 2
  register_feature_transform("test_doubler", custom_transform, "Doubles all values", overwrite = TRUE)

  reg <- feature_transform_registry()
  expect_true("test_doubler" %in% names(reg))
  expect_identical(reg$test_doubler$fun, custom_transform)
  expect_equal(reg$test_doubler$description, "Doubles all values")

  # Clean up
  rm("test_doubler", envir = .transform_registry)
})

test_that("register_feature_transform prevents overwriting without overwrite=TRUE", {
  expect_error(
    register_feature_transform("none", function(x) x, overwrite = FALSE),
    "already registered"
  )
})

test_that("register_feature_transform validates inputs", {
  expect_error(register_feature_transform(123, function(x) x), "is.character")
  expect_error(register_feature_transform("", function(x) x), "nzchar")
  expect_error(register_feature_transform("test", "not a function"), "is.function")
})

test_that(".get_feature_transform retrieves registered transforms", {
  transform <- .get_feature_transform("none")
  expect_type(transform, "list")
  expect_true("fun" %in% names(transform))
  expect_true("description" %in% names(transform))
})

test_that(".get_feature_transform errors on unknown transform", {
  expect_error(
    .get_feature_transform("nonexistent_transform"),
    "is not registered"
  )
})

test_that("transform_none returns input unchanged", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- paste0("F", 1:4)
  result <- transform_none(x)
  expect_identical(result, x)
})

test_that("transform_pairwise_ratios calls pairwise_col_diff", {
  x <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, ncol = 2)
  colnames(x) <- c("A", "B")
  result <- transform_pairwise_ratios(x)

  expect_equal(ncol(result), 1)
  expect_true(grepl("--", colnames(result)[1]))
  expect_equal(result[, 1], x[, "A"] - x[, "B"])
})

test_that("transform_pairwise_ratios generates all pairs", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  result <- transform_pairwise_ratios(x)

  expected_pairs <- 4 * 3 / 2
  expect_equal(ncol(result), expected_pairs)
  expect_true(all(grepl("--", colnames(result))))
})

test_that("transform_pairwise_ratios errors on single feature", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "A"

  expect_error(
    transform_pairwise_ratios(x),
    "requires at least two features"
  )
})

test_that("transform_pairwise_log_ratios computes log ratios", {
  x <- matrix(c(10, 20, 5, 10), nrow = 2, ncol = 2)
  colnames(x) <- c("A", "B")
  result <- transform_pairwise_log_ratios(x)

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result)[1], "A--B")
  expect_equal(result[, 1], log(x[, "A"] / x[, "B"]))
})

test_that("transform_pairwise_log_ratios generates all pairs", {
  x <- matrix(abs(rnorm(12)) + 1, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  result <- transform_pairwise_log_ratios(x)

  expected_pairs <- 4 * 3 / 2
  expect_equal(ncol(result), expected_pairs)
  expect_true(all(grepl("--", colnames(result))))
})

test_that("transform_pairwise_log_ratios errors on single feature", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "A"

  expect_error(
    transform_pairwise_log_ratios(x),
    "requires at least two features"
  )
})

test_that("transform_pairwise_log_ratios requires column names", {
  x <- matrix(1:6, nrow = 3, ncol = 2)

  expect_error(
    transform_pairwise_log_ratios(x),
    "must have column names"
  )
})

test_that("transform_reference_norm requires reference_feature attribute", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")

  expect_error(
    transform_reference_norm(x),
    "requires 'reference_feature' attribute"
  )
})

test_that("transform_reference_norm normalizes to reference", {
  x <- matrix(c(10, 20, 30, 5, 10, 15, 2, 4, 6), nrow = 3, ncol = 3)
  colnames(x) <- c("A", "B", "Ref")
  attr(x, "reference_feature") <- "Ref"

  result <- transform_reference_norm(x)

  expect_equal(ncol(result), 2)
  expect_true(all(c("A--Ref", "B--Ref") %in% colnames(result)))
  expect_equal(result[, "A--Ref"], x[, "A"] - x[, "Ref"])
  expect_equal(result[, "B--Ref"], x[, "B"] - x[, "Ref"])
})

test_that("transform_reference_norm validates reference feature exists", {
  x <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(x) <- c("A", "B", "C", "D")
  attr(x, "reference_feature") <- "NonExistent"

  expect_error(
    transform_reference_norm(x),
    "not found in data"
  )
})

test_that("transform_reference_norm warns when only reference feature exists", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "Ref"
  attr(x, "reference_feature") <- "Ref"

  expect_warning(
    result <- transform_reference_norm(x),
    "requires at least one non-reference feature"
  )
  expect_identical(result, x)
})

test_that(".apply_feature_transform applies transform to all matrices", {
  mat1 <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(mat1) <- c("A", "B")
  mat2 <- matrix(7:12, nrow = 3, ncol = 2)
  colnames(mat2) <- c("A", "B")

  result <- .apply_feature_transform(list(mat1, mat2), "pairwise_ratios")

  expect_length(result, 2)
  expect_equal(ncol(result[[1]]), 1)
  expect_equal(ncol(result[[2]]), 1)
})

test_that(".apply_feature_transform handles empty list", {
  result <- .apply_feature_transform(list(), "none")
  expect_length(result, 0)
})

test_that(".apply_feature_transform works with none transform", {
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(mat) <- c("A", "B")

  result <- .apply_feature_transform(list(mat), "none")

  expect_length(result, 1)
  expect_identical(result[[1]], mat)
})

test_that(".apply_feature_transform_single applies transform to single matrix", {
  x <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(x) <- c("A", "B")

  result <- .apply_feature_transform_single(x, "pairwise_ratios")

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result)[1], "A--B")
  expect_equal(result[, 1], x[, "A"] - x[, "B"])
})

test_that(".apply_feature_transform_single returns unchanged for none", {
  x <- matrix(1:6, nrow = 3, ncol = 2)
  colnames(x) <- c("A", "B")

  result <- .apply_feature_transform_single(x, "none")

  expect_identical(result, x)
})

test_that(".apply_feature_transform_single returns unchanged for single column", {
  x <- matrix(1:3, nrow = 3, ncol = 1)
  colnames(x) <- "A"

  # Should return unchanged since pairwise transform needs >= 2 features
  result <- .apply_feature_transform_single(x, "pairwise_ratios")

  expect_identical(result, x)
})
