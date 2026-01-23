# Tests for OptimizationResult class, fit_panel(), and summarize_solutions()

test_that("OptimizationResult class structure is correct", {
  # Create minimal OptimizationResult
  solutions_df <- data.frame(
    solution_id = 1:2,
    sensitivity = c(0.8, 0.9),
    specificity = c(0.7, 0.6),
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list(c("g1", "g2"), c("g1", "g3")))

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = c("g1", "g2", "g3"),
    control = list(max_features = 3),
    training_signature = list(n = 100, p = 3),
    aggregated_x = matrix(rnorm(300), nrow = 100, ncol = 3,
                          dimnames = list(NULL, c("g1", "g2", "g3"))),
    aggregated_y = factor(rep(c("No", "Yes"), each = 50), levels = c("No", "Yes")),
    aggregated_cohort = factor(rep("cohort_01", 100))
  )

  expect_s4_class(opt, "OptimizationResult")
  expect_equal(nrow(solutions(opt)), 2L)
  expect_equal(n_solutions(opt), 2L)
  expect_equal(length(opt@feature_pool), 3L)
})

test_that("get_solution_features accessor works", {
  solutions_df <- data.frame(
    solution_id = 1:3,
    sensitivity = c(0.8, 0.9, 0.85),
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list(c("g1", "g2"), c("g1", "g3"), c("g2", "g3")))

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = c("g1", "g2", "g3"),
    control = list(),
    training_signature = list(),
    aggregated_x = NULL,
    aggregated_y = NULL,
    aggregated_cohort = NULL
  )

  expect_equal(get_solution_features(opt, 1), c("g1", "g2"))
  expect_equal(get_solution_features(opt, 2), c("g1", "g3"))
  expect_equal(get_solution_features(opt, 3), c("g2", "g3"))

  expect_error(get_solution_features(opt, 0), "must be between 1 and")
  expect_error(get_solution_features(opt, 4), "must be between 1 and")
})

test_that("summarize_solutions returns correct format", {
  solutions_df <- data.frame(
    solution_id = 1:3,
    sensitivity = c(0.8, 0.9, 0.85),
    specificity = c(0.7, 0.6, 0.75),
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list(c("g1", "g2"), c("g1"), c("g2", "g3", "g4")))

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = c("g1", "g2", "g3", "g4"),
    control = list(),
    training_signature = list(),
    aggregated_x = NULL,
    aggregated_y = NULL,
    aggregated_cohort = NULL
  )

  summary <- summarize_solutions(opt)

  expect_s3_class(summary, "data.frame")
  expect_equal(nrow(summary), 3L)
  expect_true("solution_id" %in% names(summary))
  expect_true("n_features" %in% names(summary))
  expect_true("sensitivity" %in% names(summary))
  expect_true("specificity" %in% names(summary))
  expect_equal(summary$n_features, c(2L, 1L, 3L))
})

test_that("summarize_solutions handles empty OptimizationResult", {
  solutions_df <- data.frame(
    solution_id = integer(),
    sensitivity = numeric(),
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list())

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = character(),
    control = list(),
    training_signature = list(),
    aggregated_x = NULL,
    aggregated_y = NULL,
    aggregated_cohort = NULL
  )

  summary <- summarize_solutions(opt)
  expect_equal(nrow(summary), 0L)
})

test_that("fit_panel selects solution automatically", {
  skip_slow_tests()
  set.seed(123)

  sim <- simulate_expression_data(p = 30, n = 40, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  expect_s4_class(opt, "OptimizationResult")

  # Auto-select best solution
  panel <- fit_panel(opt)

  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_true(length(panel@features) >= 1)
  expect_true(!is.null(panel@model))
})

test_that("fit_panel selects specified solution", {
  skip_slow_tests()
  set.seed(456)

  sim <- simulate_expression_data(p = 30, n = 40, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  n_sol <- n_solutions(opt)
  expect_true(n_sol >= 1)

  # Select specific solution
  panel <- fit_panel(opt, solution_id = 1)

  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_equal(panel@features, get_solution_features(opt, 1))
})

test_that("fit_panel accepts explicit features", {
  skip_slow_tests()
  set.seed(789)

  sim <- simulate_expression_data(p = 30, n = 40, k = 1, seed = 44)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  # Use explicit features
  explicit_features <- colnames(x)[1:2]
  panel <- fit_panel(opt, features = explicit_features)

  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_equal(panel@features, explicit_features)
})

test_that("fit_panel with regularized = TRUE produces cv.glmnet model", {
  skip_slow_tests()
  set.seed(321)

  sim <- simulate_expression_data(p = 30, n = 40, k = 1, seed = 45)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  panel <- fit_panel(opt, regularized = TRUE)

  expect_true(inherits(panel@model, "cv.glmnet"))
  expect_true(panel@control$regularized)
})

test_that("fit_panel with regularized = FALSE produces glm model", {
  skip_slow_tests()
  set.seed(654)

  sim <- simulate_expression_data(p = 30, n = 40, k = 1, seed = 46)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    regularized = FALSE,
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  panel <- fit_panel(opt, regularized = FALSE)

  expect_true(inherits(panel@model, "glm"))
  expect_false(panel@control$regularized)
})

test_that("fit_panel errors on invalid solution_id", {
  solutions_df <- data.frame(
    solution_id = 1:2,
    sensitivity = c(0.8, 0.9),
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list(c("g1", "g2"), c("g1", "g3")))

  set.seed(123)
  x_mat <- matrix(rnorm(300), nrow = 100, ncol = 3,
                  dimnames = list(NULL, c("g1", "g2", "g3")))
  y_vec <- factor(rep(c("No", "Yes"), each = 50), levels = c("No", "Yes"))

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = c("g1", "g2", "g3"),
    control = list(objective_directions = list(sensitivity = "maximize")),
    training_signature = list(),
    aggregated_x = x_mat,
    aggregated_y = y_vec,
    aggregated_cohort = factor(rep("cohort_01", 100))
  )

  expect_error(fit_panel(opt, solution_id = 5), "not found")
})

test_that("fit_panel errors on features not in pool", {
  solutions_df <- data.frame(
    solution_id = 1L,
    sensitivity = 0.8,
    stringsAsFactors = FALSE
  )
  solutions_df$features <- I(list(c("g1", "g2")))

  set.seed(123)
  x_mat <- matrix(rnorm(300), nrow = 100, ncol = 3,
                  dimnames = list(NULL, c("g1", "g2", "g3")))
  y_vec <- factor(rep(c("No", "Yes"), each = 50), levels = c("No", "Yes"))

  opt <- new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = c("g1", "g2", "g3"),
    control = list(objective_directions = list(sensitivity = "maximize")),
    training_signature = list(),
    aggregated_x = x_mat,
    aggregated_y = y_vec,
    aggregated_cohort = factor(rep("cohort_01", 100))
  )

  expect_error(fit_panel(opt, features = c("g1", "not_in_pool")), "not in feature pool")
})

test_that("new workflow: optimize_panel -> fit_panel -> evaluate_panel", {
  skip_slow_tests()
  set.seed(999)

  sim <- simulate_expression_data(p = 30, n = 60, k = 1, seed = 47)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  # Split into train/test
  train_idx <- seq_len(40)
  test_idx <- seq(41, 60)
  x_train <- x[train_idx, ]
  y_train <- y[train_idx]
  x_test <- x[test_idx, ]
  y_test <- y[test_idx]

  # Step 1: Optimize
  opt <- optimize_panel(
    x = x_train,
    y = y_train,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(10)],
    cohort_aggregator = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )
  expect_s4_class(opt, "OptimizationResult")

  # Step 2: Summarize and fit
  summary <- summarize_solutions(opt)
  expect_s3_class(summary, "data.frame")
  expect_true(nrow(summary) >= 1)

  panel <- fit_panel(opt)
  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_true(!is.null(panel@model))

  # Step 3: Evaluate
  eval_res <- evaluate_panel(panel, x_test, y_test)
  expect_type(eval_res, "list")
  expect_true("metrics" %in% names(eval_res))
  expect_true("roc" %in% names(eval_res))
  expect_true(all(eval_res$scores >= 0 & eval_res$scores <= 1))
})
