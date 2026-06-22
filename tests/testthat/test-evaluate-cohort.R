# Tests for evaluating trained biomarker panels against multiple held-out validation cohorts and cross-cohort comparisons.
test_that("evaluate_panel scores held-out multi-cohort data", {
  panel_features <- c("gene_common1", "gene_common2")

  make_data <- function(seed, cols) {
    set.seed(seed)
    n <- 30L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 1.1 * x[, "gene_common1"] - 0.6 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  # Create training and validation data
  train <- make_data(100, panel_features)
  v1 <- make_data(101, c("gene_common1", "gene_common2", "gene_extra1"))
  v2 <- make_data(202, c("gene_common2", "gene_common1", "gene_extra2"))

  # Create panel with fitted model
  panel <- .create_test_panel_with_model(panel_features, train$x, train$y)

  eval_res <- evaluate_panel(
    panel = panel,
    x = list(v1$x, v2$x),
    y = list(v1$y, v2$y)
  )

  expect_named(eval_res$metrics, c("sensitivity", "specificity", "auc"))
  expect_true(all(is.finite(eval_res$metrics)))
  expect_true(is.matrix(eval_res$confusion))
  expect_equal(dim(eval_res$confusion), c(2, 2))
  expect_true(all(c("No", "Yes") %in% rownames(eval_res$confusion)))
  expect_true(all(c("No", "Yes") %in% colnames(eval_res$confusion)))
  expect_true(is.list(eval_res$roc))
  expect_true(all(c("curve", "highlight", "auc", "plot", "evalm") %in% names(eval_res$roc)))
  expect_s3_class(eval_res$roc$curve, "data.frame")
  expect_true(all(c("threshold", "tpr", "fpr", "sensitivity", "specificity") %in% colnames(eval_res$roc$curve)))
  expect_s3_class(eval_res$roc$highlight, "data.frame")
  expect_equal(nrow(eval_res$roc$highlight), 1L)
  expect_true(is.null(eval_res$roc$plot) || inherits(eval_res$roc$plot, "ggplot"))
  expect_s3_class(eval_res$roc$evalm, "data.frame")
  expect_true(all(c("obs", "prob") %in% colnames(eval_res$roc$evalm)))
})

# ============================================================================
# evaluate_panel_by_cohort tests
# ============================================================================

test_that("evaluate_panel_by_cohort returns per-cohort metrics", {
  panel_features <- c("gene_common1", "gene_common2")

  make_cohort <- function(seed, n = 30) {
    set.seed(seed)
    x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(x) <- c("gene_common1", "gene_common2", "gene_extra")
    linear <- 1.1 * x[, "gene_common1"] - 0.6 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  train <- make_cohort(100)
  v1 <- make_cohort(101)
  v2 <- make_cohort(202)
  v3 <- make_cohort(303)

  panel <- .create_test_panel_with_model(panel_features, train$x, train$y)

  result <- evaluate_panel_by_cohort(
    panel = panel,
    x = list(v1$x, v2$x, v3$x),
    y = list(v1$y, v2$y, v3$y)
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_true(all(c("cohort", "n_samples", "n_positive", "n_negative",
                    "sensitivity", "specificity", "auc") %in% names(result)))

  expect_true(all(result$sensitivity >= 0 & result$sensitivity <= 1))
  expect_true(all(result$specificity >= 0 & result$specificity <= 1))
  expect_true(all(result$auc >= 0 & result$auc <= 1))
  expect_equal(result$n_samples, c(30L, 30L, 30L))

  expect_equal(attr(result, "cutoff_prob"), 0.5)
  expect_equal(attr(result, "positive"), "Yes")
})

test_that("evaluate_panel_by_cohort handles single-cohort input", {
  panel_features <- c("g1", "g2")

  set.seed(123)
  x <- matrix(rnorm(60), nrow = 30, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(rep(c("No", "Yes"), each = 15), levels = c("No", "Yes"))

  panel <- .create_test_panel_with_model(panel_features, x, y)

  result <- evaluate_panel_by_cohort(
    panel = panel,
    x = x,
    y = y
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$cohort, "cohort_01")
})

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

test_that("plot_cohort_comparison handles method column", {
  skip_if_not_installed("ggplot2")

  cohort_metrics <- data.frame(
    cohort = rep(c("A", "B"), each = 2),
    method = rep(c("MOO", "LASSO"), 2),
    sensitivity = c(0.92, 0.85, 0.86, 0.82),
    specificity = c(0.76, 0.80, 0.81, 0.78)
  )

  p <- plot_cohort_comparison(cohort_metrics)

  expect_s3_class(p, "ggplot")
})

test_that("plot_cohort_comparison respects facet_by argument", {
  skip_if_not_installed("ggplot2")

  cohort_metrics <- data.frame(
    cohort = c("A", "B"),
    sensitivity = c(0.9, 0.85),
    specificity = c(0.8, 0.75)
  )

  p1 <- plot_cohort_comparison(cohort_metrics, facet_by = "metric")
  p2 <- plot_cohort_comparison(cohort_metrics, facet_by = "cohort")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_cohort_comparison errors without cohort column", {
  skip_if_not_installed("ggplot2")

  cohort_metrics <- data.frame(
    sensitivity = c(0.9, 0.85),
    specificity = c(0.8, 0.75)
  )

  expect_error(
    plot_cohort_comparison(cohort_metrics),
    "cohort"
  )
})
