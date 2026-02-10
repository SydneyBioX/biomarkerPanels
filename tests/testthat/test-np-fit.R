# Tests for Neyman-Pearson panel fitting functions

test_that(".prepare_np_labels with minimize_FPR = TRUE uses standard encoding", {
  truth <- factor(c("No", "No", "Yes", "Yes", "No"), levels = c("No", "Yes"))

  result <- .prepare_np_labels(truth, minimize_FPR = TRUE)

  expect_equal(result$y_nproc, c(0L, 0L, 1L, 1L, 0L))
  expect_false(result$labels_flipped)
  expect_equal(result$original_positive, "Yes")
})

test_that(".prepare_np_labels with minimize_FPR = FALSE uses flipped encoding", {
  truth <- factor(c("No", "No", "Yes", "Yes", "No"), levels = c("No", "Yes"))

  result <- .prepare_np_labels(truth, minimize_FPR = FALSE)

  # Flipped: Yes=0, No=1
  expect_equal(result$y_nproc, c(1L, 1L, 0L, 0L, 1L))
  expect_true(result$labels_flipped)
  expect_equal(result$original_positive, "Yes")
})

test_that(".invert_np_predictions inverts when labels_flipped = TRUE", {
  preds <- c(0.1, 0.3, 0.7, 0.9)

  # No inversion when not flipped
  result_no_flip <- .invert_np_predictions(preds, labels_flipped = FALSE)
  expect_equal(result_no_flip, preds)

  # Inversion when flipped
  result_flip <- .invert_np_predictions(preds, labels_flipped = TRUE)
  expect_equal(result_flip, c(0.9, 0.7, 0.3, 0.1))
})

test_that(".validate_np_params catches invalid inputs", {
  # Valid params should not error
  expect_silent(.validate_np_params(alpha = 0.05, delta = 0.05, method = "penlog"))
  expect_silent(.validate_np_params(alpha = 0.1, delta = 0.1, method = "svm"))

  # Invalid alpha
  expect_error(
    .validate_np_params(alpha = 0, delta = 0.05, method = "penlog"),
    "alpha"
  )
  expect_error(
    .validate_np_params(alpha = 1, delta = 0.05, method = "penlog"),
    "alpha"
  )
  expect_error(
    .validate_np_params(alpha = -0.1, delta = 0.05, method = "penlog"),
    "alpha"
  )
  expect_error(
    .validate_np_params(alpha = c(0.05, 0.1), delta = 0.05, method = "penlog"),
    "alpha"
  )

  # Invalid delta
  expect_error(
    .validate_np_params(alpha = 0.05, delta = 0, method = "penlog"),
    "delta"
  )
  expect_error(
    .validate_np_params(alpha = 0.05, delta = 1, method = "penlog"),
    "delta"
  )

  # Invalid method
  expect_error(
    .validate_np_params(alpha = 0.05, delta = 0.05, method = "invalid"),
    "method"
  )
})

test_that(".check_nproc_available works", {
  # This test only verifies the function exists and handles the check
  # If nproc is installed, it should return invisibly
  # If not installed, it will error (skip in that case)
  skip_if_not_installed("nproc")
  expect_invisible(.check_nproc_available())
})

# Integration tests that require nproc package
test_that("fit_np_panel with minimize_FPR = TRUE returns TransferablePanelResult", {
  skip_if_not_installed("nproc")
  skip_if_not_installed("rmoo")
  skip_slow_tests()

  # Create small test data
  set.seed(42)
  sim <- simulate_expression_data(p = 50, n = 40, k = 2, seed = 42)

  # Get a small feature pool
  feature_pool <- colnames(sim$x_list[[1]])[1:20]

  # Use list format for multi-cohort (optimize_panel infers cohort from list structure)
  opt <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    feature_pool = feature_pool,
    feature_transform = "none",
    objectives = define_objectives(
      metrics = c("sensitivity", "specificity")
    ),
    nsga_control = list(popSize = 20, maxiter = 5)
  )

  # Fit NP panel
  panel <- fit_np_panel(
    opt,
    minimize_FPR = TRUE,
    alpha = 0.15,
    delta = 0.1,
    method = "penlog"
  )

  expect_s4_class(panel, "TransferablePanelResult")
  expect_true(length(panel@features) > 0)
  expect_true(is.numeric(panel@np_threshold))
  expect_equal(panel@np_alpha, 0.15)
  expect_equal(panel@np_delta, 0.1)
  expect_false(panel@control$labels_flipped)
  expect_true(panel@control$minimize_FPR)
  expect_s3_class(panel@model, "npc")
})

test_that("fit_np_panel with minimize_FPR = FALSE flips labels correctly", {
  skip_if_not_installed("nproc")
  skip_if_not_installed("rmoo")
  skip_slow_tests()

  set.seed(123)
  sim <- simulate_expression_data(p = 50, n = 40, k = 2, seed = 123)

  feature_pool <- colnames(sim$x_list[[1]])[1:20]

  opt <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    feature_pool = feature_pool,
    feature_transform = "none",
    objectives = define_objectives(
      metrics = c("sensitivity", "specificity")
    ),
    nsga_control = list(popSize = 20, maxiter = 5)
  )

  # Fit NP panel for rule-out (control FNR)
  panel <- fit_np_panel(
    opt,
    minimize_FPR = FALSE,
    alpha = 0.10,  # Target 90% sensitivity
    delta = 0.1,
    method = "penlog"
  )

  expect_s4_class(panel, "TransferablePanelResult")
  expect_true(panel@control$labels_flipped)
  expect_false(panel@control$minimize_FPR)

  # Per-cohort metrics should be computed
  expect_true(nrow(panel@per_cohort_metrics) >= 1)
  expect_true("sensitivity" %in% names(panel@per_cohort_metrics))
  expect_true("specificity" %in% names(panel@per_cohort_metrics))
})

test_that("fit_np_panel with explicit features works", {
  skip_if_not_installed("nproc")
  skip_if_not_installed("rmoo")
  skip_slow_tests()

  set.seed(456)
  sim <- simulate_expression_data(p = 50, n = 40, k = 2, seed = 456)

  feature_pool <- colnames(sim$x_list[[1]])[1:20]

  opt <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    feature_pool = feature_pool,
    feature_transform = "none",
    objectives = define_objectives(
      metrics = c("sensitivity", "specificity")
    ),
    nsga_control = list(popSize = 20, maxiter = 5)
  )

  # Fit with explicit features
  explicit_features <- c("gene_0001", "gene_0003", "gene_0005")
  panel <- fit_np_panel(
    opt,
    features = explicit_features,
    minimize_FPR = TRUE,
    alpha = 0.1,
    delta = 0.1,
    method = "penlog"
  )

  expect_equal(panel@features, explicit_features)
  expect_true(is.na(panel@control$fitted_solution_id))
})

test_that("evaluate_panel works with npc model", {
  skip_if_not_installed("nproc")
  skip_if_not_installed("rmoo")
  skip_slow_tests()

  set.seed(789)
  sim <- simulate_expression_data(p = 50, n = 40, k = 3, seed = 789)

  # Use first 2 cohorts for training
  x_train_list <- sim$x_list[1:2]
  y_train_list <- sim$y_list[1:2]

  # Use third cohort for testing
  x_test <- sim$x_list[[3]]
  y_test <- sim$y_list[[3]]

  feature_pool <- colnames(sim$x_list[[1]])[1:20]

  opt <- optimize_panel(
    x = x_train_list,
    y = y_train_list,
    feature_pool = feature_pool,
    feature_transform = "none",
    objectives = define_objectives(
      metrics = c("sensitivity", "specificity")
    ),
    nsga_control = list(popSize = 20, maxiter = 5)
  )

  panel <- fit_np_panel(
    opt,
    minimize_FPR = TRUE,
    alpha = 0.15,
    delta = 0.1,
    method = "penlog"
  )

  # Evaluate on held-out data
  eval_result <- evaluate_panel(panel, x_test, y_test)

  expect_type(eval_result, "list")
  expect_true("metrics" %in% names(eval_result))
  expect_true("confusion" %in% names(eval_result))
  expect_true("roc" %in% names(eval_result))
  expect_true("scores" %in% names(eval_result))
  expect_equal(length(eval_result$scores), nrow(x_test))
})

test_that("fit_np_panel errors without nproc package", {
  # This test mocks the nproc check to verify error handling
  # Skip if nproc IS installed (we want to test the error path)
  skip_if(requireNamespace("nproc", quietly = TRUE),
          "nproc is installed, cannot test missing package error")

  # Create minimal mock optimization result
  mock_opt <- new("OptimizationResult",
    solutions = data.frame(
      solution_id = 1L,
      sensitivity = 0.9,
      features = I(list(c("gene_0001", "gene_0002")))
    ),
    feature_pool = c("gene_0001", "gene_0002", "gene_0003"),
    control = list(),
    training_signature = list(),
    aggregated_x = matrix(1:6, nrow = 3, dimnames = list(NULL, c("gene_0001", "gene_0002"))),
    aggregated_y = factor(c("No", "Yes", "Yes"), levels = c("No", "Yes")),
    aggregated_cohort = factor(c("A", "A", "A"))
  )

  expect_error(
    fit_np_panel(mock_opt),
    "nproc.*not installed"
  )
})

test_that("fit_np_panel validates optimization_result type", {
  skip_if_not_installed("nproc")

  expect_error(
    fit_np_panel("not an optimization result"),
    "OptimizationResult"
  )

  expect_error(
    fit_np_panel(list(solutions = data.frame())),
    "OptimizationResult"
  )
})
