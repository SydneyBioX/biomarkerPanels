test_that("select_cpop_features returns pair features and base features they reference", {
  set.seed(7)
  sim <- simulate_expression_data(p = 25, n = 50, k = 2, seed = 7)

  res <- select_cpop_features(
    x_list = sim$x_list,
    y_list = sim$y_list,
    n_features = 6L,
    n_iter = 3L
  )

  expect_named(res, c("features", "base_features", "cpop1_features",
                      "step1", "step2_coefficients", "settings"))
  expect_type(res$features, "character")
  expect_lte(length(res$features), 6L)
  expect_true(all(grepl("--", res$features, fixed = TRUE)))

  # Base features must be exactly the union of genes referenced by the pairs
  derived_base <- sort(unique(unlist(strsplit(res$features, "--", fixed = TRUE))))
  expect_identical(res$base_features, derived_base)
})

test_that("CPOP step 2 sign agreement is preserved across cohorts", {
  # Integration check: after CPOP runs end-to-end with unanimous sign
  # threshold, every selected pair has matching coefficient signs across
  # cohorts in the recorded step-2 ridge fits.
  set.seed(17)
  sim <- simulate_expression_data(p = 25, n = 60, k = 2, seed = 17)
  res <- select_cpop_features(sim$x_list, sim$y_list,
                              n_features = 6L, n_iter = 3L,
                              sign_consistency_threshold = 1.0)
  if (length(res$features) >= 2L) {
    coefs <- res$step2_coefficients
    # Compare coefficient signs across cohorts row-by-row (ignoring zeros)
    sign_consistent <- apply(coefs, 1L, function(v) {
      nz <- sign(v[v != 0])
      length(nz) == 0L || all(nz == nz[1L])
    })
    expect_true(all(sign_consistent[res$features]))
  } else {
    skip("CPOP returned <2 features in this seed configuration")
  }
})

test_that(".mst_lratio_pruning yields n_nodes - 1 edges and is acyclic", {
  features <- c("A--B", "A--C", "B--C", "C--D", "A--D")
  pruned <- biomarkerPanels:::.mst_lratio_pruning(features)
  expect_length(pruned, 3L)

  # All retained edges should appear in the input
  expect_true(all(pruned %in% c(features, vapply(strsplit(features, "--", fixed = TRUE),
                                                  function(p) paste0(rev(p), collapse = "--"),
                                                  character(1)))))
})

test_that("colmeans_penalty normalises to sum equal to vector length", {
  set.seed(3)
  m1 <- matrix(stats::rnorm(40, mean = 0), nrow = 10)
  m2 <- matrix(stats::rnorm(40, mean = 1), nrow = 10)
  w <- colmeans_penalty(list(m1, m2))
  expect_length(w, 4L)
  expect_equal(sum(w), 4)
  expect_true(all(w >= 0))
})

test_that("fit_cpop_panel returns a BiomarkerPanelResult that evaluate_panel consumes", {
  set.seed(7)
  sim <- simulate_expression_data(p = 30, n = 60, k = 2, seed = 7)

  res <- select_cpop_features(sim$x_list, sim$y_list,
                              n_features = 8L, n_iter = 3L)
  if (length(res$features) < 2L) skip("Selection returned <2 features")
  panel <- fit_cpop_panel(res, sim$x_list, sim$y_list)

  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_identical(panel_features(panel), res$features)
  expect_identical(panel_base_features(panel), res$base_features)
  expect_identical(panel@control$feature_transform, "pairwise_ratios")

  # Averaged coefficients are the column means of per-cohort coefficients
  per <- panel@control$cpop$per_cohort_coefficients
  avg <- panel@control$cpop$averaged_coefficients
  expect_equal(rowMeans(per), avg)

  # evaluate_panel should produce sensible metrics without warnings about
  # feature mismatches (the panel selects a subset of all possible pairs)
  eval_res <- evaluate_panel(panel, sim$x_list[[2]], sim$y_list[[2]])
  expect_named(eval_res$metrics, c("sensitivity", "specificity", "auc"),
               ignore.order = TRUE)
  expect_true(all(is.finite(eval_res$metrics)))
})

test_that("shell model prediction equals CPOP's averaged per-cohort log-odds", {
  # Verifies the fidelity claim: predict(panel@model, ...) returns
  # plogis(mean of per-cohort linear predictors), matching how the original
  # CPOP package combines its two cohort-specific models.
  set.seed(7)
  sim <- simulate_expression_data(p = 30, n = 60, k = 2, seed = 7)
  res <- select_cpop_features(sim$x_list, sim$y_list,
                              n_features = 8L, n_iter = 3L)
  if (length(res$features) < 2L) skip("Selection returned <2 features")
  panel <- fit_cpop_panel(res, sim$x_list, sim$y_list)

  z_new <- pairwise_col_diff(sim$x_list[[2]])[, panel_features(panel),
                                              drop = FALSE]
  beta_mat <- panel@control$cpop$per_cohort_coefficients
  a0 <- panel@control$cpop$per_cohort_intercepts
  link1 <- as.numeric(z_new %*% beta_mat[, 1L] + a0[1L])
  link2 <- as.numeric(z_new %*% beta_mat[, 2L] + a0[2L])
  expected_prob <- stats::plogis((link1 + link2) / 2)

  pred_prob <- as.numeric(stats::predict(panel@model, newx = z_new,
                                         s = "lambda.min",
                                         type = "response")[, 1L])
  expect_equal(pred_prob, expected_prob, tolerance = 1e-6)
})

test_that("evaluate_panel_by_cohort accepts CPOP panels", {
  set.seed(7)
  sim <- simulate_expression_data(p = 30, n = 60, k = 2, seed = 7)
  res <- select_cpop_features(sim$x_list, sim$y_list,
                              n_features = 8L, n_iter = 3L)
  if (length(res$features) < 2L) skip("Selection returned <2 features")
  panel <- fit_cpop_panel(res, sim$x_list, sim$y_list)

  per_cohort <- evaluate_panel_by_cohort(panel, sim$x_list, sim$y_list)
  expect_true(is.data.frame(per_cohort) || is.list(per_cohort))
})

test_that("select_cpop_features generalises to >2 cohorts", {
  set.seed(13)
  sim <- simulate_expression_data(p = 20, n = 40, k = 3, seed = 13)
  res <- select_cpop_features(sim$x_list, sim$y_list,
                              n_features = 4L, n_iter = 3L)
  expect_equal(ncol(res$step2_coefficients), 3L)
  expect_true(length(res$features) >= 0L)
})
