# Tests for feature aggregation logic, primarily focusing on contrast feature generation and pairwise cohort aggregators.
test_that("pairwise cohort aggregator produces contrast features", {
  skip_slow_tests()
  set.seed(1234)
  make_mat <- function(seed) {
    set.seed(seed)
    m <- matrix(rnorm(600), nrow = 200, ncol = 3)
    colnames(m) <- c("A", "B", "C")
    # Add signal so cv.glmnet selects a feature
    m[1:100, "A"] <- m[1:100, "A"] + 0.5
    m
  }
  x_list <- list(make_mat(1), make_mat(2))
  y_list <- lapply(1:2, function(i) factor(rep(c("No", "Yes"), length.out = 200), levels = c("No", "Yes")))

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    regularized = TRUE,
    feature_transform = "pairwise_ratios",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "pairwise_ratios")
  
  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(grepl("--", all_features, fixed = TRUE)))
})

test_that("feature_pool accepts base features with pairwise aggregator", {
  skip_slow_tests()
  set.seed(777)
  x <- matrix(rnorm(600), nrow = 200, ncol = 3)
  colnames(x) <- c("GeneA", "GeneB", "GeneC")
  y <- factor(rep(c("No", "Yes"), each = 100), levels = c("No", "Yes"))
  x[1:100, "GeneA"] <- x[1:100, "GeneA"] + 0.5

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = c("GeneA", "GeneC", "GeneB"),
    feature_transform = "pairwise_ratios",
    regularized = TRUE,
    fitness_mode = "in_sample",
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(sort(res@control$base_feature_pool), sort(c("GeneA", "GeneB", "GeneC")))

  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(grepl("--", all_features, fixed = TRUE)))
  components <- unique(unlist(strsplit(all_features, "--", fixed = TRUE)))
  expect_true(all(components %in% c("GeneA", "GeneB", "GeneC")))
})

test_that("optimize_panel works with pairwise_log_ratios aggregator", {
  skip_slow_tests()
  set.seed(456)
  sim <- simulate_expression_data(p = 20, n = 200, k = 1, seed = 43)
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
    regularized = FALSE,
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "pairwise_log_ratios")
  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(grepl("--", all_features, fixed = TRUE)))
})

test_that("optimize_panel works with reference_norm aggregator", {
  skip_slow_tests()
  set.seed(789)
  sim <- simulate_expression_data(p = 15, n = 200, k = 1, seed = 44)
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
    regularized = FALSE,
    fitness_mode = "in_sample",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "reference_norm")
})

test_that("optimize_panel rejects unregistered aggregator", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 10, n = 200, k = 1, seed = 45)

  expect_error(
    optimize_panel(
      x = sim$x_list[[1]],
      y = sim$y_list[[1]],
      objectives = define_objectives(metrics = c("sensitivity", "specificity")),
      max_features = 2,
      feature_transform = "nonexistent_transform",
      nsga_control = list(popSize = 12, maxiter = 5)
    ),
    "Unknown feature transform"
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
  sim <- simulate_expression_data(p = 10, n = 200, k = 1, seed = 46)

  res <- optimize_panel(
    x = sim$x_list[[1]],
    y = sim$y_list[[1]],
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "center_features",
    fitness_mode = "in_sample",
    nsga_control = list(popSize = 12, maxiter = 6)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_transform, "center_features")

  # Clean up
  rm("center_features", envir = biomarkerPanels:::.transform_registry)
})
