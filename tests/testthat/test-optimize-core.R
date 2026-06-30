# Core tests for optimize_panel, validating basic structural outputs, bounds checking, and handling of base scenarios.
test_that("optimize_panel returns an OptimizationResult", {
  skip_slow_tests()
  set.seed(123)

  sim <- simulate_expression_data(p = 30, n = 200, k = 1, seed = 42)
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
