test_that("optimize_panel returns an OptimizationResult", {
  skip_slow_tests()
  set.seed(123)

  sim <- simulate_expression_data(p = 30, n = 25, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_true(n_solutions(res) >= 1)
  expect_true(length(res@feature_pool) > 0)

  # Check solutions data frame structure
  sols <- solutions(res)
  expect_true("solution_id" %in% names(sols))
  expect_true("features" %in% names(sols))
  expect_true("sensitivity" %in% names(sols))
  expect_true("specificity" %in% names(sols))
  expect_true("num_features" %in% names(sols))

  # Check training data is stored
  expect_true(!is.null(res@aggregated_x))
  expect_true(!is.null(res@aggregated_y))
})

test_that("num_features objective allows variable panel sizes", {
  skip_slow_tests()
  set.seed(999)
  # Create data with varying signal strength across features
  n <- 40L
  p <- 10L
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  # Strong signal on first 2 genes
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 3
  # Moderate signal on genes 3-5
  x[y == "Yes", 3:5] <- x[y == "Yes", 3:5] + 1.5
  # Weak signal on genes 6-8
  x[y == "Yes", 6:8] <- x[y == "Yes", 6:8] + 0.5

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
  n <- 20L
  p <- 6L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("g", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 2

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 4,
    feature_transform = "none",
    nsga_control = list(popSize = 8, maxiter = 5)
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

test_that("optimize_panel handles multiple cohorts", {
  skip_slow_tests()
  skip_if_not(file.exists(fixture_path("fake_gene_expression.Rds")))
  fixture <- readRDS(fixture_path("fake_gene_expression.Rds"))

  res <- optimize_panel(
    x = fixture$x_list,
    y = fixture$y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 4,
    feature_pool = colnames(fixture$x_list[[1]])[seq_len(12)],
    feature_transform = "none",
    nsga_control = list(popSize = 16, maxiter = 15)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_true(n_solutions(res) >= 1)

  # Check feature count constraint
  sols <- solutions(res)
  for (i in seq_len(nrow(sols))) {
    expect_true(length(sols$features[[i]]) <= 4)
  }
})

test_that("optimize_panel intersects feature sets across cohorts", {
  skip_slow_tests()
  make_cohort <- function(seed, cols) {
    set.seed(seed)
    n <- 40L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 0.9 * x[, "gene_common1"] - 0.7 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  c1 <- make_cohort(1, c("gene_common1", "gene_common2", "gene_unique1"))
  c2 <- make_cohort(2, c("gene_common2", "gene_common1", "gene_unique2"))

  res <- optimize_panel(
    x = list(c1$x, c2$x),
    y = list(c1$y, c2$y),
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")

  # All features in solutions should be from the intersection
  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(all_features %in% c("gene_common1", "gene_common2")))
})

test_that("min_metric_constraint builds feasible constraint", {
  constraint <- min_metric_constraint("sensitivity", threshold = 0.9)
  expect_type(constraint$label, "character")
  expect_true(is.function(constraint$fun))
  expect_equal(constraint$metric, "sensitivity")
  expect_equal(constraint$direction, "maximize")
})

test_that("optimize_panel enforces minimum metric constraints", {
  skip_slow_tests()
  set.seed(5001)
  n <- 60L
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("gene_common1", "gene_common2")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", "gene_common1"] <- x[y == "Yes", "gene_common1"] + 3
  x[y == "No", "gene_common1"] <- x[y == "No", "gene_common1"] - 3

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("specificity")),
    max_features = 1,
    feature_pool = colnames(x),
    constraints = list(min_metric_constraint("sensitivity", threshold = 0.9)),
    feature_transform = "none",
    nsga_control = list(popSize = 16, maxiter = 15)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_true("min_sensitivity_0.9" %in% res@control$constraints)

  # Fit and evaluate to verify constraint is satisfied
  panel <- fit_panel(res)
  eval_res <- evaluate_panel(panel, x, y)
  expect_gte(eval_res$metrics["sensitivity"], 0.9)
})

test_that("optimize_panel errors when constraints infeasible", {
  skip_slow_tests()
  set.seed(5002)
  n <- 40L
  x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(x) <- c("gene_common1", "gene_common2", "gene_common3")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

  expect_error(
    optimize_panel(
      x = x,
      y = y,
      objectives = define_objectives(metrics = c("specificity")),
      max_features = 2,
      feature_pool = colnames(x),
      constraints = list(min_metric_constraint("sensitivity", threshold = 1.01)),
      feature_transform = "none",
      nsga_control = list(popSize = 12, maxiter = 10)
    ),
    "No solutions satisfied the supplied constraints"
  )
})

test_that("pairwise cohort aggregator produces contrast features", {
  skip_slow_tests()
  set.seed(1234)
  make_mat <- function(seed) {
    set.seed(seed)
    m <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(m) <- c("A", "B", "C")
    m
  }
  x_list <- list(make_mat(1), make_mat(2))
  y_list <- lapply(1:2, function(i) factor(rep(c("No", "Yes"), length.out = 10), levels = c("No", "Yes")))

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "pairwise_ratios",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "pairwise_ratios")
  expect_true(all(grepl("--", res@control$feature_pool)))
})

test_that("feature_pool accepts base features with pairwise aggregator", {
  skip_slow_tests()
  set.seed(777)
  x <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(x) <- c("GeneA", "GeneB", "GeneC")
  y <- factor(sample(c("No", "Yes"), 20, replace = TRUE), levels = c("No", "Yes"))

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 1,
    feature_pool = c("GeneA", "GeneC"),
    feature_transform = "pairwise_ratios",
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(sort(res@control$base_feature_pool), sort(c("GeneA", "GeneC")))

  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(grepl("--", all_features, fixed = TRUE)))
  components <- unique(unlist(strsplit(all_features, "--", fixed = TRUE)))
  expect_true(all(components %in% c("GeneA", "GeneC")))
})

test_that("optimize_panel works with pairwise_log_ratios aggregator", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 20, n = 30, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  # Ensure positive values for log-ratios
  x <- abs(x) + 1
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "pairwise_log_ratios",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "pairwise_log_ratios")
  expect_true(all(grepl("--", res@control$feature_pool)))
})

test_that("optimize_panel works with reference_norm aggregator", {
  skip_slow_tests()
  set.seed(789)
  sim <- simulate_expression_data(p = 15, n = 25, k = 1, seed = 44)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Set reference feature attribute
  attr(x, "reference_feature") <- colnames(x)[1]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "reference_norm",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "reference_norm")
})

test_that("optimize_panel rejects unregistered aggregator", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 45)

  expect_error(
    optimize_panel(
      x = sim$x_list[[1]],
      y = sim$y_list[[1]],
      objectives = define_objectives(metrics = c("sensitivity", "specificity")),
      max_features = 2,
      feature_transform = "nonexistent_transform",
      nsga_control = list(popSize = 8, maxiter = 5)
    ),
    "Unknown aggregator"
  )
})

test_that("custom aggregator can be registered and used", {
  skip_slow_tests()
  # Register custom aggregator
  custom_agg <- function(x) {
    # Simple centering aggregator
    result <- scale(x, center = TRUE, scale = FALSE)
    colnames(result) <- colnames(x)
    result
  }
  register_feature_transform("center_features", custom_agg, "Center each feature", overwrite = TRUE)

  set.seed(321)
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 46)

  res <- optimize_panel(
    x = sim$x_list[[1]],
    y = sim$y_list[[1]],
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "center_features",
    nsga_control = list(popSize = 12, maxiter = 6)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "center_features")

  # Clean up
  rm("center_features", envir = .transform_registry)
})

test_that(".compute_nsga3_partitions returns appropriate values", {
  # 2-3 objectives: 12 partitions
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(2), 12L)
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(3), 12L)

  # 4-5 objectives: 6 partitions
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(4), 6L)
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(5), 6L)

  # 6+ objectives: 4 partitions
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(6), 4L)
  expect_equal(biomarkerPanels:::.compute_nsga3_partitions(10), 4L)
})

test_that("optimize_panel uses adaptive defaults without explicit nsga_control", {
  skip_extended_tests()  # Very slow: popSize=156, maxiter=180
  set.seed(999)
  sim <- simulate_expression_data(p = 40, n = 30, k = 1, seed = 99)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Use 35 features which falls into the 31-100 range
  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(35)],
    feature_transform = "none"
  )

  expect_s4_class(res, "OptimizationResult")
  # Default algorithm is NSGA-III
  expect_equal(res@control$algorithm, "NSGA-III")
  # Verify that the stored nsga2 settings reflect adaptive NSGA-III defaults
  expect_equal(res@control$nsga2$popSize, 156)
  expect_equal(res@control$nsga2$maxiter, 180)
})

# Issue 4: feature_alignment tests
test_that("feature_alignment = 'majority' keeps features in >= 50% cohorts", {
  skip_slow_tests()
  # Create 4 cohorts with varying feature overlap
  make_cohort <- function(seed, features) {
    set.seed(seed)
    n <- 40L
    x <- matrix(rnorm(n * length(features)), nrow = n, ncol = length(features))
    colnames(x) <- features
    linear <- 0.8 * x[, "common1"] + 0.5 * x[, "common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohorts <- list(
    make_cohort(1, c("common1", "common2", "partial1", "partial2", "unique1")),
    make_cohort(2, c("common1", "common2", "partial1", "partial2", "unique2")),
    make_cohort(3, c("common1", "common2", "partial1", "unique3")),
    make_cohort(4, c("common1", "common2", "unique4"))
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "none",
    feature_alignment = "majority",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_alignment, "majority")

  pool_features <- res@control$feature_pool
  expect_true(all(c("common1", "common2") %in% pool_features))
  expect_true("partial1" %in% pool_features)  # 3/4 = 75%
  expect_true("partial2" %in% pool_features)  # 2/4 = 50%
  expect_false(any(c("unique1", "unique2", "unique3", "unique4") %in% pool_features))
})

test_that("feature_alignment = 'intersection' is default and drops partial features", {
  skip_slow_tests()
  make_cohort <- function(seed, features) {
    set.seed(seed)
    n <- 30L
    x <- matrix(rnorm(n * length(features)), nrow = n, ncol = length(features))
    colnames(x) <- features
    y <- factor(sample(c("No", "Yes"), n, replace = TRUE), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohorts <- list(
    make_cohort(1, c("geneA", "geneB", "geneC")),
    make_cohort(2, c("geneA", "geneB", "geneD"))
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "none",
    feature_alignment = "intersection",
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  expect_equal(res@control$feature_alignment, "intersection")
  expect_equal(sort(res@control$feature_pool), c("geneA", "geneB"))
})

# Regularized scoring tests
test_that("regularized = TRUE uses regularized scoring during optimization", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 20, n = 40, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    feature_transform = "none",
    regularized = TRUE,
    regularized_alpha = 0.5,
    nsga_control = list(popSize = 16, maxiter = 10)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_true(res@control$regularized)
  expect_equal(res@control$regularized_alpha, 0.5)
})

test_that("regularized = FALSE uses unregularized scoring", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 20, n = 40, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    feature_transform = "none",
    regularized = FALSE,
    nsga_control = list(popSize = 16, maxiter = 10)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_false(res@control$regularized)
  expect_null(res@control$regularized_alpha)
})

# NSGA algorithm tests
test_that("optimize_panel uses NSGA-III by default", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 15, n = 25, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$algorithm, "NSGA-III")
})

test_that("optimize_panel respects explicit NSGA-II algorithm selection", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 15, n = 25, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    algorithm = "NSGA-II",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$algorithm, "NSGA-II")
})

# Adaptive threshold tests
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
  n <- 60L
  p <- 30L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  # Strong signal on first 5 genes
  x[y == "Yes", 1:5] <- x[y == "Yes", 1:5] + 2

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
  n <- 40L
  p <- 15L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:3] <- x[y == "Yes", 1:3] + 2

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
  n <- 20L
  p <- 6L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("g", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:2] <- x[y == "Yes", 1:2] + 2

  # Test with adaptive
  res_adaptive <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "num_features")),
    max_features = 4,
    feature_transform = "none",
    selection_threshold = "adaptive",
    nsga_control = list(popSize = 8, maxiter = 5)
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
    nsga_control = list(popSize = 8, maxiter = 5)
  )
  expect_equal(res_fixed@control$selection_threshold, 0.6)
})

test_that("Pareto front contains no dominated solutions", {
  skip_slow_tests()
  set.seed(777)
  n <- 40L
  p <- 10L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:3] <- x[y == "Yes", 1:3] + 2

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 6,
    feature_pool = colnames(x),
    feature_transform = "none",
    nsga_control = list(popSize = 20, maxiter = 20)
  )

  sols <- solutions(res)
  n_sol <- nrow(sols)
  if (n_sol < 2) skip("Only one Pareto solution; dominance check trivially passes")

  obj_cols <- setdiff(names(sols), c("solution_id", "base_features", "features"))
  directions <- res@control$objective_directions[obj_cols]

  # For each pair of solutions, verify neither dominates the other on ALL objectives
  for (i in seq_len(n_sol - 1)) {
    for (j in (i + 1):n_sol) {
      vals_i <- as.numeric(sols[i, obj_cols])
      vals_j <- as.numeric(sols[j, obj_cols])

      # Compare with direction: "maximize" means higher is better
      i_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vi > vj else vi < vj
      }, vals_i, vals_j, directions)

      i_equal_or_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vi >= vj else vi <= vj
      }, vals_i, vals_j, directions)

      j_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vj > vi else vj < vi
      }, vals_i, vals_j, directions)

      j_equal_or_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vj >= vi else vj <= vi
      }, vals_i, vals_j, directions)

      # i dominates j if i is >= on all and > on at least one
      i_dominates_j <- all(i_equal_or_better) && any(i_better)
      j_dominates_i <- all(j_equal_or_better) && any(j_better)

      expect_false(
        i_dominates_j,
        info = sprintf("Solution %d dominates solution %d", i, j)
      )
      expect_false(
        j_dominates_i,
        info = sprintf("Solution %d dominates solution %d", j, i)
      )
    }
  }
})

# ============================================================================
# Input validation tests (silent failure audit fixes)
# ============================================================================

test_that("optimize_panel validates fitness_cv_folds", {
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 1)
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], fitness_cv_folds = 1L),
    "fitness_cv_folds"
  )
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], fitness_cv_folds = "abc"),
    "fitness_cv_folds"
  )
})

test_that("optimize_panel validates regularized_alpha", {
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 1)
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], regularized_alpha = 2),
    "regularized_alpha"
  )
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], regularized_alpha = -0.1),
    "regularized_alpha"
  )
})

test_that("optimize_panel validates selection_threshold", {
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 1)
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], selection_threshold = "bad"),
    "selection_threshold"
  )
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], selection_threshold = 0),
    "selection_threshold"
  )
  expect_error(
    optimize_panel(sim$x_list[[1]], sim$y_list[[1]], selection_threshold = 1),
    "selection_threshold"
  )
})
