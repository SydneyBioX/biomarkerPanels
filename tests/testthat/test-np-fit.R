# Tests for Neyman-Pearson panel fitting functions

test_that(".prepare_np_labels with minimize_FPR = TRUE uses standard encoding", {
  truth <- factor(c("No", "No", "Yes", "Yes", "No"), levels = c("No", "Yes"))

  result <- .prepare_np_labels(truth, minimize_FPR = TRUE)

  expect_equal(result$y_nproc, c(0L, 0L, 1L, 1L, 0L))
  expect_false(result$labels_flipped)
})

test_that(".prepare_np_labels with minimize_FPR = FALSE uses flipped encoding", {
  truth <- factor(c("No", "No", "Yes", "Yes", "No"), levels = c("No", "Yes"))

  result <- .prepare_np_labels(truth, minimize_FPR = FALSE)

  # Flipped: Yes=0, No=1
  expect_equal(result$y_nproc, c(1L, 1L, 0L, 0L, 1L))
  expect_true(result$labels_flipped)
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

# Integration tests that require nproc package
test_that("fit_np_panel with explicit features works", {
  skip_if_not_installed("nproc")
  skip_if_not_installed("rmoo")
  skip_slow_tests()

  set.seed(456)
  sim <- simulate_expression_data(p = 50, n = 200, k = 2, seed = 456)

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

  # Fit with explicit features (suppress warnings from small simulated n)
  explicit_features <- c("gene_0001", "gene_0003", "gene_0005")
  panel <- suppressWarnings(fit_np_panel(
    opt,
    features = explicit_features,
    minimize_FPR = TRUE,
    alpha = 0.1,
    delta = 0.1,
    method = "penlog"
  ))

  expect_equal(panel@features, explicit_features)
})
