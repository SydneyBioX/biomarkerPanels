# Tests for evaluating a panel's performance on a single dataset, testing scoring functions and confusion matrices.
test_that("evaluate_panel accepts custom cutoff for confusion and ROC highlight", {
  panel_features <- c("g1", "g2")

  set.seed(123)
  x <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("g1", "g2")))
  linear <- 1.2 * x[, "g1"] - 0.4 * x[, "g2"]
  prob <- stats::plogis(linear)
  y <- factor(ifelse(runif(20) < prob, "Yes", "No"), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  eval_res <- evaluate_panel(
    panel = panel,
    x = x,
    y = y,
    cutoff_prob = 0.3
  )

  highlight <- eval_res$roc$highlight
  expect_equal(highlight$sensitivity, metric_sensitivity(y, eval_res$scores, cutoff_prob = 0.3))
  expect_equal(highlight$specificity, metric_specificity(y, eval_res$scores, cutoff_prob = 0.3))
})

# ============================================================================
# Input validation tests (silent failure audit fixes)
# ============================================================================

test_that("evaluate_panel stops on feature name mismatch (no positional fallback)", {
  panel_features <- c("g1", "g2")

  set.seed(42)
  x <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  # Same column count but different names should stop(), not warn()
  x_misnamed <- x
  colnames(x_misnamed) <- c("a1", "a2")
  expect_error(
    evaluate_panel(panel, x_misnamed, y),
    "not found in validation data|Feature name mismatch"
  )
})

test_that("evaluate_panel treats a data.frame as a single cohort, not a cohort list", {
  panel_features <- c("g1", "g2")

  set.seed(321)
  x <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, panel_features))
  prob <- stats::plogis(1.2 * x[, "g1"] - 0.4 * x[, "g2"])
  y <- factor(ifelse(runif(20) < prob, "Yes", "No"), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  # Character response + data.frame predictors (the analysis-notebook pattern)
  y_chr <- as.character(y)
  eval_df <- evaluate_panel(panel, as.data.frame(x), y_chr)
  eval_mat <- evaluate_panel(panel, x, y_chr)

  expect_equal(eval_df$scores, eval_mat$scores)
  expect_equal(eval_df$metrics, eval_mat$metrics)
})

test_that(".is_cohort_list distinguishes cohort lists from data.frames", {
  expect_true(.is_cohort_list(list(matrix(1:4, 2))))
  expect_false(.is_cohort_list(data.frame(a = 1:2)))
  expect_false(.is_cohort_list(matrix(1:4, 2)))
  expect_false(.is_cohort_list(NULL))
})
