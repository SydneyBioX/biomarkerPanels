# Tests for underlying NSGA-III algorithm internals, such as reference point generation and partition parameters.
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

test_that("optimize_panel uses NSGA-III by default", {
  skip_slow_tests()
  set.seed(123)
  sim <- simulate_expression_data(p = 15, n = 200, k = 1, seed = 42)
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
  sim <- simulate_expression_data(p = 15, n = 200, k = 1, seed = 43)
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
