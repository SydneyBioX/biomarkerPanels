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
  sim <- simulate_expression_data(p = 200L, n = 50L, k = 3L, seed = 42L)
  x_train <- sim$x_list[1:2]
  y_train <- sim$y_list[1:2]
  x_test  <- sim$x_list[[3]]
  y_test  <- sim$y_list[[3]]

  top_feats <- get_top_de_features(x_train, y_train, n_features = 30)
  objectives <- define_objectives(
    losses = c("num_features", "sensitivity", "specificity")
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

test_that("connect = TRUE adds a geom_line layer", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  p_no  <- plot_pareto_front(fix$opt, fix$x_test, fix$y_test,
                             connect = FALSE, verbose = FALSE)
  p_yes <- plot_pareto_front(fix$opt, fix$x_test, fix$y_test,
                             connect = TRUE,  verbose = FALSE)

  layer_classes <- function(pl) {
    vapply(pl$layers, function(l) class(l$geom)[1], character(1))
  }
  expect_false("GeomLine" %in% layer_classes(p_no))
  expect_true("GeomLine"  %in% layer_classes(p_yes))
})

test_that("custom labels are passed through", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  p <- plot_pareto_front(
    fix$opt, fix$x_test, fix$y_test,
    title = "My Title", xlab = "Xlab", ylab = "Ylab",
    color_label = "Nfeat", verbose = FALSE
  )

  expect_equal(p$labels$title, "My Title")
  expect_equal(p$labels$x, "Xlab")
  expect_equal(p$labels$y, "Ylab")
})

test_that("on_error = 'stop' aborts on failure", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  # Provide x_test with wrong column names to force evaluation failure
  bad_x <- fix$x_test
  colnames(bad_x) <- paste0("WRONG_", seq_len(ncol(bad_x)))

  expect_error(
    plot_pareto_front(fix$opt, bad_x, fix$y_test,
                      on_error = "stop", verbose = FALSE),
    "Failed on solution"
  )
})

test_that("on_error = 'warn' skips failed solutions", {
  skip_slow_tests()
  skip_if_not_installed("ggplot2")

  fix <- .make_pareto_fixture()

  # Provide x_test with wrong column names — all solutions will fail
  bad_x <- fix$x_test
  colnames(bad_x) <- paste0("WRONG_", seq_len(ncol(bad_x)))

  # All fail => error about all solutions failing
  expect_error(
    suppressWarnings(
      plot_pareto_front(fix$opt, bad_x, fix$y_test,
                        on_error = "warn", verbose = FALSE)
    ),
    "All solutions failed"
  )
})
