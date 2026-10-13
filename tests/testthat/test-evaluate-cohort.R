# Tests for evaluating trained biomarker panels against multiple held-out validation cohorts and cross-cohort comparisons.
# ============================================================================
# evaluate_panel_by_cohort tests
# ============================================================================

test_that("evaluate_panel_by_cohort preserves cohort names", {
  panel_features <- c("g1", "g2")

  make_cohort <- function(seed) {
    set.seed(seed)
    x <- matrix(rnorm(40), nrow = 20, ncol = 2)
    colnames(x) <- c("g1", "g2")
    y <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  train <- make_cohort(0)
  v1 <- make_cohort(1)
  v2 <- make_cohort(2)
  v3 <- make_cohort(3)

  panel <- .create_test_panel_with_model(panel_features, train$x, train$y)

  result <- evaluate_panel_by_cohort(
    panel = panel,
    x = list(train = v1$x, valid = v2$x, test = v3$x),
    y = list(train = v1$y, valid = v2$y, test = v3$y)
  )

  expect_equal(result$cohort, c("train", "valid", "test"))
})

# ============================================================================
# plot_cohort_comparison tests
# ============================================================================

test_that("plot_cohort_comparison creates ggplot object", {
  skip_if_not_installed("ggplot2")

  cohort_metrics <- data.frame(
    cohort = c("A", "B", "C"),
    n_samples = c(50, 48, 52),
    n_positive = c(25, 22, 28),
    n_negative = c(25, 26, 24),
    sensitivity = c(0.92, 0.86, 0.89),
    specificity = c(0.76, 0.81, 0.79),
    auc = c(0.89, 0.87, 0.88)
  )

  p <- plot_cohort_comparison(cohort_metrics)

  expect_s3_class(p, "ggplot")
})
