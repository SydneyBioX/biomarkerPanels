# Tests covering how objective and constraint functions influence feature selection, such as varying panel sizes and directionality.
test_that("num_features objective allows variable panel sizes", {
  skip_slow_tests()
  set.seed(999)
  # Create data with varying signal strength across features
  n <- 200L
  p <- 10L
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  # Strong signal on first 2 genes
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 0.5
  # Moderate signal on genes 3-5
  x[y == "Yes", 3:5] <- x[y == "Yes", 3:5] + 0.3
  # Weak signal on genes 6-8
  x[y == "Yes", 6:8] <- x[y == "Yes", 6:8] + 0.1

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 8,
    feature_pool = colnames(x),
    feature_transform = "none",
    nsga_control = list(popSize = 16, maxiter = 15)
  )

  expect_s4_class(res, "OptimizationResult")

  # Extract num_features values from all Pareto solutions
  sols <- solutions(res)
  nf_values <- sols$num_features

  # All sizes should be between 2 (min for regularized) and max_features
  expect_true(all(nf_values >= 2),
              info = paste("Some panel sizes < 2:", paste(nf_values[nf_values < 2], collapse = ", ")))
  expect_true(all(nf_values <= 8),
              info = paste("Some panel sizes > 8:", paste(nf_values[nf_values > 8], collapse = ", ")))

  # With threshold-based selection and varying signal, we expect some variation
  unique_sizes <- unique(nf_values)
  message("Panel sizes in Pareto front: ", paste(sort(unique_sizes), collapse = ", "))

  # Primary check: mechanism allows variable sizes (not all forced to max_features)
  expect_true(
    min(nf_values) < 8,
    info = "Threshold selection should allow panels smaller than max_features"
  )
})

test_that("all Pareto solutions have consistent feature counts", {
  set.seed(123)
  n <- 200L
  p <- 6L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("g", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 0.5

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 4,
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  sols <- solutions(res)
  expect_true(nrow(sols) > 0)

  # Check each solution's num_features matches actual base feature count
  for (i in seq_len(nrow(sols))) {
    stored_count <- sols$num_features[i]
    base_feats <- sols$base_features[[i]]
    expect_equal(
      stored_count, length(base_feats),
      info = paste("Solution", i, "base feature count mismatch:",
                   "stored =", stored_count, "actual =", length(base_feats))
    )
  }
})

test_that("min_metric_constraint builds feasible constraint", {
  constraint <- min_metric_constraint("sensitivity", threshold = 0.9)
  expect_type(constraint$label, "character")
  expect_true(is.function(constraint$fun))
  expect_equal(constraint$metric, "sensitivity")
  expect_equal(constraint$direction, "maximize")
})

test_that(".generate_sparse_suggestions creates diverse panel sizes", {
  set.seed(123)
  suggestions <- biomarkerPanels:::.generate_sparse_suggestions(
    n_features = 20,
    n_suggestions = 10,
    min_features = 2,
    max_features = 15
  )

  expect_true(is.matrix(suggestions))
  expect_equal(ncol(suggestions), 20)
  expect_true(nrow(suggestions) >= 2)

  # All values should be in [0, 1]
  expect_true(all(suggestions >= 0 & suggestions <= 1))

  # Each row should have a mix of high and low values
  for (i in seq_len(nrow(suggestions))) {
    high_weights <- sum(suggestions[i, ] > 0.5)
    low_weights <- sum(suggestions[i, ] < 0.5)
    expect_true(high_weights >= 2)
    expect_true(low_weights >= 1)
  }
})

test_that("adaptive threshold produces variable panel sizes", {
  skip_slow_tests()
  set.seed(123)
  n <- 200L
  p <- 30L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  # Strong signal on first 5 genes
  x[y == "Yes", 1:5] <- x[y == "Yes", 1:5] + 0.5

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("balanced_accuracy", "num_features")),
    max_features = 15,
    feature_transform = "none",
    regularized = FALSE,
    selection_threshold = "adaptive",
    nsga_control = list(popSize = 40, maxiter = 30)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$selection_threshold, "adaptive")

  sols <- solutions(res)
  nf_values <- sols$num_features
  message("Panel sizes in Pareto front (adaptive): ", paste(sort(unique(nf_values)), collapse = ", "))

  expect_true(
    min(nf_values) < 15,
    info = "Adaptive threshold should allow panels smaller than max_features"
  )
})

test_that("fixed threshold 0.5 is backward compatible", {
  skip_slow_tests()
  set.seed(456)
  n <- 200L
  p <- 15L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:3] <- x[y == "Yes", 1:3] + 0.5

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 8,
    feature_transform = "none",
    regularized = FALSE,
    selection_threshold = 0.5,
    nsga_control = list(popSize = 20, maxiter = 15)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$selection_threshold, 0.5)

  sols <- solutions(res)
  for (i in seq_len(nrow(sols))) {
    expect_true(length(sols$features[[i]]) >= 1)
    expect_true(length(sols$features[[i]]) <= 8)
  }
})

test_that("selection_threshold stored in control slot", {
  set.seed(42)
  n <- 200L
  p <- 6L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("g", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 0.5

  # Test with adaptive
  res_adaptive <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 4,
    feature_transform = "none",
    selection_threshold = "adaptive",
    nsga_control = list(popSize = 12, maxiter = 5)
  )
  expect_equal(res_adaptive@control$selection_threshold, "adaptive")

  # Test with numeric
  res_fixed <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 4,
    feature_transform = "none",
    selection_threshold = 0.6,
    nsga_control = list(popSize = 12, maxiter = 5)
  )
  expect_equal(res_fixed@control$selection_threshold, 0.6)
})
