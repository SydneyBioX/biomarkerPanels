# Helper to create a simple panel with fitted model for testing
.create_test_panel_with_model <- function(features, x_train, y_train) {
  # Fit a simple glm model
  df <- data.frame(
    .response = as.integer(y_train) - 1L,
    as.data.frame(x_train[, features, drop = FALSE], check.names = TRUE)
  )
  model <- stats::glm(.response ~ ., data = df, family = stats::binomial())

  metrics <- c(sensitivity = 0.5, specificity = 0.5, auc = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(features), nrow(objective_df)))

  new(
    "BiomarkerPanelResult",
    features = features,
    metrics = metrics,
    objectives = objective_df,
    control = list(positive_class = "Yes"),
    training_data = list(),
    model = model
  )
}

# Simple scoring function for testing (used when no model needed)
.simple_scoring_fn <- function(x_selected, selected_features, truth, cohort = NULL, ...) {
  # Use row means as simple score
  rowMeans(x_selected)
}

test_that("evaluate_panel errors when no model and no scoring_fn", {
  panel_features <- c("gene_common1", "gene_common2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, auc = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))

  # Panel without model
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list(),
    model = NULL
  )

  set.seed(123)
  x <- matrix(rnorm(60), nrow = 30, ncol = 2)
  colnames(x) <- panel_features
  y <- factor(rep(c("No", "Yes"), each = 15), levels = c("No", "Yes"))

  expect_error(
    evaluate_panel(panel, x, y),
    "fit_panel"
  )
})

test_that("evaluate_panel works with scoring_fn when no model", {
  panel_features <- c("gene_common1", "gene_common2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, auc = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))

  # Panel without model
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list(),
    model = NULL
  )

  set.seed(123)
  x <- matrix(rnorm(60), nrow = 30, ncol = 2)
  colnames(x) <- panel_features
  y <- factor(rep(c("No", "Yes"), each = 15), levels = c("No", "Yes"))

  # Should work with custom scoring_fn
  eval_res <- evaluate_panel(panel, x, y, scoring_fn = .simple_scoring_fn)

  expect_type(eval_res, "list")
  expect_true("metrics" %in% names(eval_res))
})

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

test_that("evaluate_panel requires column names for single-cohort data", {
  panel_features <- c("gene_common1", "gene_common2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, auc = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list(),
    model = NULL  # We'll use scoring_fn
  )

  x <- matrix(rnorm(20), nrow = 10)
  y <- factor(rep(c("No", "Yes"), each = 5), levels = c("No", "Yes"))

  expect_error(
    evaluate_panel(panel, x, y, scoring_fn = .simple_scoring_fn),
    "Selected feature\\(s\\) not found"
  )
})

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
  expect_equal(highlight$threshold, 0.3)
  expect_equal(highlight$sensitivity, loss_sensitivity(y, eval_res$scores, cutoff_prob = 0.3))
  expect_equal(highlight$specificity, loss_specificity(y, eval_res$scores, cutoff_prob = 0.3))
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

# ============================================================================
# find_threshold_for_sensitivity tests
# ============================================================================

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
  n <- 100
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  # Create signal: g1 predicts outcome
  linear <- 2 * x[, "g1"]
  prob <- stats::plogis(linear)
  y <- factor(ifelse(prob > 0.5, "Yes", "No"), levels = c("No", "Yes"))

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
  n <- 100
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  linear <- 2 * x[, "g1"]
  prob <- stats::plogis(linear)
  y <- factor(ifelse(prob > 0.5, "Yes", "No"), levels = c("No", "Yes"))

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
  n <- 100
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("g1", "g2")
  y <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  # Add strong signal: Yes samples have substantially higher g1
  x[51:100, "g1"] <- x[51:100, "g1"] + 3

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
