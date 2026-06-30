# Tests for computing empirical decision thresholds to achieve targeted classification sensitivity/specificity.
test_that("find_threshold_for_sensitivity returns correct structure", {
  skip_if_not_installed("pROC")

  panel_features <- c("g1", "g2")

  # Create well-separated data
  set.seed(123)
  x <- matrix(c(rnorm(30, -1), rnorm(30, 1)), nrow = 30, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(c(rep("No", 15), rep("Yes", 15)), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  result <- find_threshold_for_sensitivity(panel, x, y, target_sensitivity = 0.90)

  expect_type(result, "list")
  expect_named(result, c("threshold", "sensitivity", "specificity",
                          "target_sensitivity", "auc"))
  expect_true(is.numeric(result$threshold))
  expect_true(result$sensitivity >= 0 && result$sensitivity <= 1)
  expect_true(result$specificity >= 0 && result$specificity <= 1)
  expect_equal(result$target_sensitivity, 0.90)
  expect_true(result$auc >= 0 && result$auc <= 1)
})

test_that("find_threshold_for_sensitivity achieves target sensitivity", {
  skip_if_not_installed("pROC")

  panel_features <- c("g1", "g2")

  # Create well-separated data for reliable ROC
  set.seed(456)
  n <- 200
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(rep(c("No", "Yes"), each = n/2), levels = c("No", "Yes"))
  x[(n/2 + 1):n, "g1"] <- x[(n/2 + 1):n, "g1"] + 2.0

  panel <- .create_test_panel_with_model(panel_features, x, y)

  result <- find_threshold_for_sensitivity(panel, x, y, target_sensitivity = 0.90)

  # Sensitivity should be at or near target (within tolerance due to discrete points)
  expect_gte(result$sensitivity, 0.85)
})

test_that("find_threshold_for_sensitivity validates inputs", {
  skip_if_not_installed("pROC")

  panel_features <- c("g1", "g2")

  set.seed(123)
  x <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  # Invalid target_sensitivity
  expect_error(
    find_threshold_for_sensitivity(panel, x, y, target_sensitivity = 1.5),
    "between 0 and 1"
  )

  expect_error(
    find_threshold_for_sensitivity(panel, x, y, target_sensitivity = -0.1),
    "between 0 and 1"
  )
})

# ============================================================================
# evaluate_panel_at_sensitivity tests
# ============================================================================

test_that("evaluate_panel_at_sensitivity combines threshold finding and evaluation", {
  skip_if_not_installed("pROC")

  panel_features <- c("g1", "g2")

  set.seed(789)
  n <- 200
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(rep(c("No", "Yes"), each = n/2), levels = c("No", "Yes"))
  x[(n/2 + 1):n, "g1"] <- x[(n/2 + 1):n, "g1"] + 2.0

  panel <- .create_test_panel_with_model(panel_features, x, y)

  result <- evaluate_panel_at_sensitivity(panel, x, y, target_sensitivity = 0.90)

  # Should have standard evaluate_panel structure
  expect_type(result, "list")
  expect_true("metrics" %in% names(result))
  expect_true("confusion" %in% names(result))
  expect_true("roc" %in% names(result))

  # Should also have target_sensitivity_info
  expect_true("target_sensitivity_info" %in% names(result))
  expect_equal(result$target_sensitivity_info$target_sensitivity, 0.90)

  # The cutoff used should match the threshold found
  expect_equal(result$cutoff, result$target_sensitivity_info$threshold)

  # The sensitivity in the result should be near target
  expect_gte(result$metrics["sensitivity"], 0.85)
})

test_that("evaluate_panel_at_sensitivity uses correct threshold vs default", {
  skip_if_not_installed("pROC")

  panel_features <- c("g1", "g2")

  # Create well-separated data for reliable threshold adjustment
  set.seed(999)
  n <- 200
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(c(rep("No", 100), rep("Yes", 100)), levels = c("No", "Yes"))
  # Add strong signal: Yes samples have substantially higher g1
  x[101:200, "g1"] <- x[101:200, "g1"] + 2.5

  panel <- .create_test_panel_with_model(panel_features, x, y)

  # Default evaluation at cutoff=0.5
  eval_default <- evaluate_panel(panel, x, y, cutoff_prob = 0.5)

  # Evaluation at target sensitivity
  eval_target <- evaluate_panel_at_sensitivity(panel, x, y, target_sensitivity = 0.90)

  # Target-based evaluation should use the computed threshold, not 0.5
  expect_equal(eval_target$cutoff, eval_target$target_sensitivity_info$threshold)

  # With well-separated data, should achieve or get close to target
  # (discrete ROC points may prevent exact match)
  expect_gte(eval_target$metrics["sensitivity"], 0.85)
})
