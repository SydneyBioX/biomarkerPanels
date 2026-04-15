# Tests for plot_pareto_front() visualization
# These tests run real optimization + evaluation so are slow -- skip on CRAN.

test_that("plot_pareto_front rejects non-OptimizationResult input", {
  skip_if_not_installed("ggplot2")
  expect_error(
    plot_pareto_front(data.frame(a = 1), matrix(1), factor("Yes")),
    "OptimizationResult"
  )
})

# --- shared fixture: small optimization result --------------------------------
# Build once, reuse across tests. Uses simulate_expression_data() from helper.

.make_pareto_fixture <- function() {
  sim <- simulate_expression_data(p = 50L, n = 200L, k = 2L, seed = 42L)
  # Use single matrix (no cohort dummies) so held-out eval works cleanly
  x_all <- sim$x_list[[1]]
  y_all <- sim$y_list[[1]]
  n <- nrow(x_all)
  train_idx <- seq_len(floor(n * 0.7))
  x_train <- x_all[train_idx, , drop = FALSE]
  y_train <- y_all[train_idx]
  x_test  <- x_all[-train_idx, , drop = FALSE]
  y_test  <- y_all[-train_idx]

  top_feats <- get_top_de_features(list(x_train), list(y_train), n_features = 30)
  objectives <- define_objectives(
    metrics = c("num_features", "sensitivity", "specificity")
  )
  opt <- optimize_panel(
    x_train, y_train,
    objectives    = objectives,
    feature_pool  = top_feats,
    max_features  = 8,
    feature_transform = "none",
    nsga_control  = list(popSize = 10, maxiter = 5),
    seed          = 42,
    fitness_cv    = FALSE
  )
  list(opt = opt, x_test = x_test, y_test = y_test)
}

test_that("plot_pareto_front returns ggplot with correct data columns", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  p <- plot_pareto_front(
    fix$opt, fix$x_test, fix$y_test,
    x_metric = "sensitivity", y_metric = "specificity",
    color_by = "n_features", verbose = FALSE
  )

  expect_s3_class(p, "ggplot")

  df <- p$data
  expect_true(is.data.frame(df))
  expect_true(all(c("solution_id", "n_features", "n_base_features",
                     "sensitivity", "specificity") %in% names(df)))
  expect_true(nrow(df) >= 1L)
  expect_true(all(df$sensitivity >= 0 & df$sensitivity <= 1))
  expect_true(all(df$specificity >= 0 & df$specificity <= 1))
})

test_that("color_by = NULL produces plot without colour mapping", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  p <- plot_pareto_front(
    fix$opt, fix$x_test, fix$y_test,
    color_by = NULL, verbose = FALSE
  )

  expect_s3_class(p, "ggplot")
  # No colour aesthetic in the default layer
  mapping_names <- names(p$mapping)
  expect_false("colour" %in% mapping_names)
})




