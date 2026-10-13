# Tests for OptimizationResult class, fit_panel(), and summarize_solutions()

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

test_that("fit_panel selects specified solution", {
  skip_slow_tests()
  set.seed(456)

  sim <- simulate_expression_data(p = 30, n = 200, k = 1, seed = 43)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  # Select specific solution
  panel <- fit_panel(opt, solution_id = 1)

  expect_equal(panel@features, get_solution_features(opt, 1))
})

test_that("fit_panel accepts explicit features", {
  skip_slow_tests()
  set.seed(789)

  sim <- simulate_expression_data(p = 30, n = 200, k = 1, seed = 44)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  # Use explicit features
  explicit_features <- colnames(x)[1:2]
  panel <- fit_panel(opt, features = explicit_features)

  expect_equal(panel@features, explicit_features)
})

test_that("fit_panel with regularized = TRUE produces cv.glmnet model", {
  skip_slow_tests()
  set.seed(321)

  sim <- simulate_expression_data(p = 30, n = 200, k = 1, seed = 45)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  opt <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 10)
  )

  panel <- fit_panel(opt, regularized = TRUE)

  expect_true(inherits(panel@model, "cv.glmnet"))
})

test_that(".fit_final_model_cv warns on insufficient data for CV", {
  set.seed(42)
  x_small <- matrix(rnorm(12), nrow = 4, ncol = 3,
                    dimnames = list(NULL, c("a", "b", "c")))
  y_small <- factor(c("No", "Yes", "No", "Yes"), levels = c("No", "Yes"))

  expect_warning(
    .fit_final_model_cv(x_small, y_small, cv_folds = 5L),
    "Too few samples"
  )
})
