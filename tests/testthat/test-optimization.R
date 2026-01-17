test_that("optimize_panel returns a BiomarkerPanelResult", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 30, n = 25, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity", "num_features")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popsize = 20, generations = 20)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(length(res@features) <= 3)
  expect_named(res@metrics, c("sensitivity", "specificity", "num_features"))
  expect_s3_class(res@objectives, "data.frame")
  expect_true(all(res@objectives$objective %in% c("sensitivity", "specificity", "num_features")))
  expect_equal(res@training_data$num_cohorts, 1L)
  expect_equal(res@training_data$cohort_labels, "cohort_01")
})

test_that("optimize_panel handles multiple cohorts", {
  skip_slow_tests()
  skip_if_not(file.exists(fixture_path("fake_gene_expression.Rds")))
  fixture <- readRDS(fixture_path("fake_gene_expression.Rds"))

  res <- optimize_panel(
    x = fixture$x_list,
    y = fixture$y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity", "num_features")),
    max_features = 4,
    feature_pool = colnames(fixture$x_list[[1]])[seq_len(12)],
    cohort_aggregator = "none",
    nsga_control = list(popsize = 16, generations = 15)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(length(res@features) <= 4)
  expect_equal(res@training_data$num_cohorts, length(fixture$x_list))
  expect_true(all(names(res@training_data$cohort_counts) %in% names(fixture$x_list)))
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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "none",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(all(res@features %in% c("gene_common1", "gene_common2")))
})

test_that("min_metric_constraint builds feasible constraint", {
  constraint <- min_metric_constraint("sensitivity", threshold = 0.9)
  expect_type(constraint$label, "character")
  expect_true(is.function(constraint$fun))
  expect_equal(constraint$loss, "sensitivity")
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
    objectives = define_objectives(losses = c("specificity")),
    max_features = 1,
    feature_pool = colnames(x),
    constraints = list(min_metric_constraint("sensitivity", threshold = 0.9)),
    cohort_aggregator = "none",
    nsga_control = list(popsize = 16, generations = 15)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true("min_sensitivity_0.9" %in% res@control$constraints)
  eval_res <- evaluate_panel(res, x, y)
  expect_gte(eval_res$metrics["sensitivity"], 0.9)
})

test_that("optimize_panel errors when constraints infeasible", {
  skip_slow_tests()
  set.seed(5002)
  n <- 40L
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("gene_common1", "gene_common2")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

  expect_error(
    optimize_panel(
      x = x,
      y = y,
      objectives = define_objectives(losses = c("specificity")),
      max_features = 1,
      feature_pool = colnames(x),
      constraints = list(min_metric_constraint("sensitivity", threshold = 1.01)),
      cohort_aggregator = "none",
      nsga_control = list(popsize = 12, generations = 10)
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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "pairwise_ratios",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$cohort_aggregator, "pairwise_ratios")
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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 1,
    feature_pool = c("GeneA", "GeneC"),
    cohort_aggregator = "pairwise_ratios",
    nsga_control = list(popsize = 12, generations = 5)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(sort(res@control$base_feature_pool), sort(c("GeneA", "GeneC")))
  expect_true(all(grepl("--", res@features, fixed = TRUE)))
  components <- unique(unlist(strsplit(res@features, "--", fixed = TRUE)))
  expect_true(all(components %in% c("GeneA", "GeneC")))
})

# Integration tests for new aggregation strategies

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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "pairwise_log_ratios",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$cohort_aggregator, "pairwise_log_ratios")
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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "reference_norm",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$cohort_aggregator, "reference_norm")
})

test_that("optimize_panel rejects unregistered aggregator", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 45)

  expect_error(
    optimize_panel(
      x = sim$x_list[[1]],
      y = sim$y_list[[1]],
      objectives = define_objectives(losses = c("sensitivity", "specificity")),
      max_features = 2,
      cohort_aggregator = "nonexistent_aggregator",
      nsga_control = list(popsize = 8, generations = 5)
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
  register_aggregator("center_features", custom_agg, "Center each feature", overwrite = TRUE)

  set.seed(321)
  sim <- simulate_expression_data(p = 10, n = 20, k = 1, seed = 46)

  res <- optimize_panel(
    x = sim$x_list[[1]],
    y = sim$y_list[[1]],
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "center_features",
    nsga_control = list(popsize = 12, generations = 6)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$cohort_aggregator, "center_features")

  # Clean up
  rm("center_features", envir = .aggregator_registry)
})

# Issue 1: Adaptive NSGA-II defaults
test_that(".get_adaptive_nsga_defaults scales with feature pool size", {
  # Small feature pool: original defaults

  small <- biomarkerPanels:::.get_adaptive_nsga_defaults(20)
  expect_equal(small$popsize, 64)

  expect_equal(small$generations, 60)

  # Medium feature pool: larger population and generations
  medium <- biomarkerPanels:::.get_adaptive_nsga_defaults(50)
  expect_equal(medium$popsize, 128)
  expect_equal(medium$generations, 150)

 # Large feature pool: even larger population and generations
  large <- biomarkerPanels:::.get_adaptive_nsga_defaults(150)
  expect_equal(large$popsize, 200)
  expect_equal(large$generations, 300)

  # All should have consistent base parameters
  expect_equal(small$cprob, 0.7)
  expect_equal(medium$cprob, 0.7)
  expect_equal(large$cprob, 0.7)
})

test_that("optimize_panel uses adaptive defaults without explicit nsga_control", {
  skip_extended_tests()  # Very slow: popsize=128, generations=150
  set.seed(999)
  sim <- simulate_expression_data(p = 40, n = 30, k = 1, seed = 99)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Use 35 features which falls into the 31-100 range
  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(35)],
    cohort_aggregator = "none"
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  # Verify that the stored nsga2 settings reflect adaptive defaults
  expect_equal(res@control$nsga2$popsize, 128)
  expect_equal(res@control$nsga2$generations, 150)
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

  # common1 and common2 in all 4 cohorts (100%)
  # partial1 in 3/4 cohorts (75%)
  # partial2 in 2/4 cohorts (50%)
  # unique features in 1/4 cohorts (25%)
  cohorts <- list(
    make_cohort(1, c("common1", "common2", "partial1", "partial2", "unique1")),
    make_cohort(2, c("common1", "common2", "partial1", "partial2", "unique2")),
    make_cohort(3, c("common1", "common2", "partial1", "unique3")),
    make_cohort(4, c("common1", "common2", "unique4"))
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  # With majority alignment, should keep common1, common2, partial1, partial2
  # (all have >= 50% = 2/4 cohorts)
  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "none",
    feature_alignment = "majority",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$feature_alignment, "majority")

  # Feature pool should include partial features but not unique ones
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
    make_cohort(2, c("geneA", "geneB", "geneD"))  # geneC missing, geneD unique
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "none",
    feature_alignment = "intersection",
    nsga_control = list(popsize = 12, generations = 5)
  )

  expect_equal(res@control$feature_alignment, "intersection")
  # Only geneA and geneB should be in the pool
  expect_equal(sort(res@control$feature_pool), c("geneA", "geneB"))
})

# Issue 5: final_model_cv tests
test_that("final_model_cv = TRUE uses cross-validated coefficients", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 15, n = 50, k = 1, seed = 789)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Without CV
  res_no_cv <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    final_model_cv = FALSE,
    nsga_control = list(popsize = 12, generations = 8)
  )

  # With CV
  res_cv <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    final_model_cv = TRUE,
    cv_folds = 5L,
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res_no_cv, "BiomarkerPanelResult")
  expect_s4_class(res_cv, "BiomarkerPanelResult")

  # Check control settings are stored correctly
  expect_false(res_no_cv@control$final_model_cv)
  expect_true(res_cv@control$final_model_cv)
  expect_null(res_no_cv@control$cv_folds)
  expect_equal(res_cv@control$cv_folds, 5L)

  # The CV model should have the cv_averaged flag
  if (!is.null(res_cv@model)) {
    expect_true(isTRUE(res_cv@model$cv_averaged))
  }
})

test_that("final_model_cv falls back when data is too small", {
  skip_slow_tests()
  set.seed(789)
  # Small dataset with only 8 samples
  n <- 8L
  x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(x) <- c("geneA", "geneB", "geneC")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

  # 5-fold CV needs at least 10 samples (2 per fold), so should fall back
  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "none",
    final_model_cv = TRUE,
    cv_folds = 5L,
    nsga_control = list(popsize = 8, generations = 5)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  # Should have fallen back to standard model fitting
  if (!is.null(res@model)) {
    expect_null(res@model$cv_averaged)
  }
})

# Regularized scoring tests
test_that("regularized = TRUE produces cv.glmnet model", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 20, n = 40, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    regularized = TRUE,
    regularized_alpha = 0.5,
    nsga_control = list(popsize = 16, generations = 10)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(res@control$regularized)
  expect_equal(res@control$regularized_alpha, 0.5)

  # The stored model should be a cv.glmnet object
  expect_true(inherits(res@model, "cv.glmnet"))

  # Should have metadata for prediction
  expect_true(!is.null(res@model$biomarkerPanels_meta))
  expect_true(!is.null(res@model$biomarkerPanels_meta$feature_names))
})

test_that("regularized = FALSE produces glm model", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 20, n = 40, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    regularized = FALSE,
    nsga_control = list(popsize = 16, generations = 10)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_false(res@control$regularized)
  expect_null(res@control$regularized_alpha)

  # The stored model should be a glm object, not cv.glmnet
  expect_false(inherits(res@model, "cv.glmnet"))
  expect_true(inherits(res@model, "glm"))
})

test_that("evaluate_panel works with regularized cv.glmnet model", {
  skip_slow_tests()
  set.seed(789)
  sim <- simulate_expression_data(p = 20, n = 60, k = 1, seed = 44)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Split into train/test
  train_idx <- seq_len(40)
  test_idx <- seq(41, 60)

  x_train <- x[train_idx, ]
  y_train <- y[train_idx]
  x_test <- x[test_idx, ]
  y_test <- y[test_idx]

  res <- optimize_panel(
    x = x_train,
    y = y_train,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(12)],
    cohort_aggregator = "none",
    regularized = TRUE,
    nsga_control = list(popsize = 16, generations = 10)
  )

  # Evaluate on test set
  eval_res <- evaluate_panel(res, x_test, y_test)

  expect_type(eval_res, "list")
  expect_true("metrics" %in% names(eval_res))
  expect_true("roc" %in% names(eval_res))

  # Scores should be valid probabilities
  expect_true(all(eval_res$scores >= 0 & eval_res$scores <= 1))
  expect_equal(length(eval_res$scores), length(y_test))
})

test_that("regularized scoring with different alpha values", {
  skip_slow_tests()
  set.seed(321)
  sim <- simulate_expression_data(p = 15, n = 40, k = 1, seed = 45)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Lasso (alpha = 1)
  res_lasso <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    regularized = TRUE,
    regularized_alpha = 1.0,
    nsga_control = list(popsize = 12, generations = 8)
  )

  # Ridge (alpha = 0)
  res_ridge <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    regularized = TRUE,
    regularized_alpha = 0.0,
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res_lasso, "BiomarkerPanelResult")
  expect_s4_class(res_ridge, "BiomarkerPanelResult")

  expect_equal(res_lasso@control$regularized_alpha, 1.0)
  expect_equal(res_ridge@control$regularized_alpha, 0.0)

  # Both should produce valid cv.glmnet models
  expect_true(inherits(res_lasso@model, "cv.glmnet"))
  expect_true(inherits(res_ridge@model, "cv.glmnet"))
})

test_that("regularized scoring handles multi-cohort data", {
  skip_slow_tests()
  skip_if_not(file.exists(fixture_path("fake_gene_expression.Rds")))
  fixture <- readRDS(fixture_path("fake_gene_expression.Rds"))

  res <- optimize_panel(
    x = fixture$x_list,
    y = fixture$y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 4,
    feature_pool = colnames(fixture$x_list[[1]])[seq_len(10)],
    cohort_aggregator = "none",
    regularized = TRUE,
    regularized_alpha = 0.5,
    nsga_control = list(popsize = 12, generations = 10)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(res@control$regularized)
  expect_true(inherits(res@model, "cv.glmnet"))

  # Should have cohort info in metadata
  if (!is.null(res@model$biomarkerPanels_meta$cohort_info)) {
    expect_true(length(res@model$biomarkerPanels_meta$cohort_info$levels) > 1)
  }
})
