# Tests for Pareto solution evaluation and plot_pareto_front() visualization
# These tests run real optimization + evaluation so are slow -- skip on CRAN.

test_that("evaluate_pareto_solutions rejects non-OptimizationResult input", {
  expect_error(
    evaluate_pareto_solutions(data.frame(a = 1), matrix(1), factor("Yes")),
    "OptimizationResult"
  )
})

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

  top_feats <- select_de_features(list(x_train), list(y_train), n_features = 30)
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
    fitness_mode  = "in_sample"
  )
  list(opt = opt, x_test = x_test, y_test = y_test)
}

test_that("evaluate_pareto_solutions returns held-out metrics", {
  skip_slow_tests()

  fix <- suppressWarnings(.make_pareto_fixture())

  df <- suppressWarnings(
    evaluate_pareto_solutions(
      fix$opt, fix$x_test, fix$y_test,
      regularized = FALSE, verbose = FALSE
    )
  )

  expect_named(df, c("solution_id", "n_features", "n_base_features",
                     "sensitivity", "specificity", "auc",
                     "obj_num_features", "obj_sensitivity", "obj_specificity"))
  expect_equal(nrow(df), n_solutions(fix$opt))
  expect_equal(df$solution_id, solutions(fix$opt)$solution_id)
  expect_true(all(df$sensitivity >= 0 & df$sensitivity <= 1))
  expect_true(all(df$specificity >= 0 & df$specificity <= 1))
  expect_true(all(df$auc >= 0 & df$auc <= 1))
})

test_that("evaluate_pareto_solutions uses stored held-out data", {
  skip_slow_tests()

  fix <- suppressWarnings(.make_pareto_fixture())
  fix$opt@control$heldout_x <- fix$x_test
  fix$opt@control$heldout_y <- fix$y_test

  expect_message(
    df <- suppressWarnings(
      evaluate_pareto_solutions(
        fix$opt, regularized = FALSE, verbose = TRUE
      )
    ),
    "Using held-out data stored in OptimizationResult"
  )
  expect_equal(nrow(df), n_solutions(fix$opt))
})

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

# --- objective re-scoring -----------------------------------------------------

test_that(".count_eval_cohorts counts list cohorts and explicit vectors", {
  m <- matrix(0, nrow = 4, ncol = 2)
  expect_equal(.count_eval_cohorts(list(m, m, m)), 3L)
  expect_equal(.count_eval_cohorts(m), 1L)
  expect_equal(.count_eval_cohorts(m, cohort = c("a", "a", "b", "b")), 2L)
})

test_that(".prepare_pareto_objectives NAs cohort metrics on one cohort", {
  objs <- define_objectives(metrics = c("auc", "cohort_auc_gap"))

  expect_warning(
    single <- .prepare_pareto_objectives(objs, n_cohorts = 1L),
    "cohort_auc_gap"
  )
  truth <- factor(rep(c("Yes", "No"), each = 5), levels = c("No", "Yes"))
  scores <- seq(0.1, 0.9, length.out = 10)
  expect_true(is.na(single$cohort_auc_gap$fun(truth, scores)))
  expect_false(is.na(single$auc$fun(truth, scores)))

  multi <- .prepare_pareto_objectives(objs, n_cohorts = 2L)
  cohort <- factor(rep(c("a", "b"), 5))
  expect_false(is.na(multi$cohort_auc_gap$fun(truth, scores, cohort = cohort)))
})

test_that(".prepare_pareto_objectives converts metric errors to NA", {
  boom <- list(boom = list(
    label = "Boom", direction = "maximize",
    fun = function(truth, estimate, selected = NULL, ...) stop("nope")
  ))
  safe <- .prepare_pareto_objectives(boom, n_cohorts = 1L)
  expect_identical(safe$boom$fun(factor("Yes"), 1), NA_real_)
})

test_that(".resolve_stored_objectives prefers stored objectives", {
  stored <- define_objectives(metrics = c("auc", "num_features"))
  opt <- new("OptimizationResult", control = list(objectives = stored))
  expect_named(.resolve_stored_objectives(opt), c("auc", "num_features"))

  # NA disables re-scoring entirely
  expect_length(.resolve_stored_objectives(opt, objectives = NA), 0L)

  # Explicit override wins
  override <- define_objectives(metrics = "pauc")
  expect_named(.resolve_stored_objectives(opt, override), "pauc")
  expect_error(.resolve_stored_objectives(opt, list(1)), "named list")
})

test_that(".resolve_stored_objectives falls back to objective_directions", {
  opt <- new("OptimizationResult", control = list(
    objective_directions = c(auc = "maximize", made_up = "minimize")
  ))
  warnings_seen <- character()
  objs <- withCallingHandlers(
    .resolve_stored_objectives(opt),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_named(objs, "auc")
  expect_true(any(grepl("made_up", warnings_seen)))
  expect_true(any(grepl("no objective definitions", warnings_seen)))

  empty <- new("OptimizationResult", control = list())
  expect_length(.resolve_stored_objectives(empty), 0L)
})

test_that("evaluate_pareto_solutions reports the optimised objectives", {
  skip_slow_tests()

  fix <- suppressWarnings(.make_pareto_fixture())

  df <- suppressWarnings(
    evaluate_pareto_solutions(
      fix$opt, fix$x_test, fix$y_test,
      regularized = FALSE, verbose = FALSE
    )
  )

  obj_cols <- paste0("obj_", names(fix$opt@control$objectives))
  expect_true(all(obj_cols %in% names(df)))
  # num_features is data-independent, so it must match the panel size exactly
  expect_equal(df$obj_num_features, df$n_features)
  expect_named(attr(df, "objective_labels"), obj_cols)
  expect_equal(
    unname(attr(df, "objective_directions")[["obj_num_features"]]),
    "minimize"
  )

  skipped <- suppressWarnings(
    evaluate_pareto_solutions(
      fix$opt, fix$x_test, fix$y_test, objectives = NA,
      regularized = FALSE, verbose = FALSE
    )
  )
  expect_false(any(grepl("^obj_", names(skipped))))
  # Descriptive columns are unaffected by which objectives are re-scored
  expect_equal(skipped$auc, df$auc)
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




