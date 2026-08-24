# Tests for .make_loco_fitness() and the leave-one-cohort-out fitness path
# of optimize_panel_transferable().

test_that(".make_loco_fitness requires at least two cohorts", {
  set.seed(21)
  x <- matrix(rnorm(24L * 3L), nrow = 24L)
  colnames(x) <- paste0("g", 1:3)
  y <- factor(rep(c("No", "Yes"), each = 12L), levels = c("No", "Yes"))

  expect_error(
    .make_loco_fitness(
      x = x, y = y, cohort = factor(rep("A", 24L)),
      feature_pool = colnames(x),
      max_features = 2L,
      objectives = define_objectives(metrics = c("auc", "num_features")),
      constraints = list(),
      regularized = FALSE,
      alpha = 0.5,
      feature_transform = "none",
      min_features_required = 1L,
      selection_threshold = 0.5
    ),
    "at least 2 cohorts"
  )
})

test_that("LOCO fitness scores every candidate out of cohort", {
  set.seed(22)
  x <- matrix(rnorm(40L * 3L), nrow = 40L)
  colnames(x) <- paste0("g", 1:3)
  y <- factor(rep(c("No", "Yes"), each = 20L), levels = c("No", "Yes"))
  cohort <- factor(rep(c("A", "B"), length.out = 40L))

  # An objective that just counts the rows it is handed: with LOCO every
  # sample appears exactly once in the pooled out-of-cohort predictions.
  n_scored <- list(
    n_scored = list(
      label = "n_scored",
      direction = "maximize",
      fun = function(truth, scores, selected = NULL, ...) length(scores)
    )
  )

  fit <- .make_loco_fitness(
    x = x, y = y, cohort = cohort,
    feature_pool = colnames(x),
    max_features = 2L,
    objectives = n_scored,
    constraints = list(),
    regularized = FALSE,
    alpha = 0.5,
    feature_transform = "none",
    min_features_required = 1L,
    selection_threshold = 0.5
  )

  pop <- rbind(
    c(0.9, 0.1, 0.2),
    c(0.9, 0.1, 0.2),
    c(0.1, 0.9, 0.2)
  )
  values <- fit$wrapper(pop)
  # Objectives are negated for maximization, so all rows score -40.
  expect_equal(values[, 1], rep(-40, 3L))
})

# Integration: LOCO fitness runs end-to-end on a small fixture
# (slow — uses NSGA). Kept brief, mirroring test-rotating-validation.R.

test_that("optimize_panel_transferable runs with fitness_mode = 'loco'", {
  skip_slow_tests()
  set.seed(11)
  sim <- simulate_expression_data(p = 40L, n = 40L, k = 2L, seed = 11L,
                                  shift_scale = 0.3)
  pool <- head(sim$metadata$informative_genes, 10L)

  result <- optimize_panel_transferable(
    x = sim$x_list, y = sim$y_list,
    feature_pool = pool,
    max_features = 3L,
    train_ratio = 0.6,
    val_ratio = 0.25,
    fitness_mode = "loco",
    objectives = define_objectives(
      metrics = c("auc", "num_features")
    ),
    nsga_control = list(popSize = 12L, maxiter = 3L),
    seed = 11L,
    feature_transform = "none"
  )

  expect_s4_class(result, "OptimizationResult")
  expect_gt(n_solutions(result), 0L)
  expect_identical(result@control$fitness_mode, "loco")
  expect_null(result@control$n_val_splits)
})

test_that("optimize_panel_transferable rejects fitness_mode = 'loco' on one cohort", {
  skip_slow_tests()
  set.seed(12)
  sim <- simulate_expression_data(p = 30L, n = 60L, k = 1L, seed = 12L,
                                  shift_scale = 0.3)
  pool <- head(sim$metadata$informative_genes, 8L)

  expect_error(
    optimize_panel_transferable(
      x = sim$x_list, y = sim$y_list,
      feature_pool = pool,
      max_features = 3L,
      train_ratio = 0.6,
      val_ratio = 0.25,
      fitness_mode = "loco",
      objectives = define_objectives(metrics = c("auc", "num_features")),
      nsga_control = list(popSize = 12L, maxiter = 3L),
      seed = 12L,
      feature_transform = "none"
    ),
    "requires >=2 cohorts"
  )
})
