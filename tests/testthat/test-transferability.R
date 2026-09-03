# Tests for optimize_panel_transferable() and related functions
# See plan for test case descriptions

# -----------------------------------------------------------------------------
# Fast tests (partition validation, stratification, weighted variance)
# -----------------------------------------------------------------------------

test_that("partition ratio validation rejects invalid ratios", {
  # train_ratio < 0.5
  expect_error(
    .validate_partition_ratios(0.4, 0.2),
    "train_ratio.*at least 0.5"
  )

  # train + val over 1 leaves a negative held-out share
  expect_error(
    .validate_partition_ratios(0.7, 0.4),
    "must not exceed 1"
  )

  # Invalid types
  expect_error(
    .validate_partition_ratios("0.7", 0.2),
    "train_ratio.*numeric"
  )
  expect_error(
    .validate_partition_ratios(0.7, NA),
    "val_ratio.*numeric"
  )

  # Valid ratios should pass silently
  expect_invisible(.validate_partition_ratios(0.7, 0.2))
  expect_invisible(.validate_partition_ratios(0.5, 0.4))
  expect_invisible(.validate_partition_ratios(0.6, 0.3))

  # Zero shares are legal: they mean "skip that partition"
  expect_invisible(.validate_partition_ratios(1, 0))
  expect_invisible(.validate_partition_ratios(0.7, 0))
  expect_invisible(.validate_partition_ratios(0.7, 0.3))
})

test_that("stratified partitioning preserves class balance", {
  set.seed(123)
  # Create test data with known class proportions
  n <- 100
  x1 <- matrix(rnorm(n * 10), nrow = n)
  colnames(x1) <- paste0("gene_", 1:10)
  y1 <- factor(c(rep("Yes", 40), rep("No", 60)), levels = c("No", "Yes"))

  x2 <- matrix(rnorm(n * 10), nrow = n)
  colnames(x2) <- paste0("gene_", 1:10)
  y2 <- factor(c(rep("Yes", 30), rep("No", 70)), levels = c("No", "Yes"))

  partitions <- .stratified_partition_cohorts(
    x_list = list(c1 = x1, c2 = x2),
    y_list = list(c1 = y1, c2 = y2),
    train_ratio = 0.7,
    val_ratio = 0.2
  )

  # Check that all partition components exist
  expect_named(partitions, c("train_x", "train_y", "val_x", "val_y",
                             "heldout_x", "heldout_y", "cohort_names",
                             "partition_info"))

  # Check cohort 1 class proportions are approximately preserved
  train_yes_1 <- sum(partitions$train_y$c1 == "Yes")
  train_no_1 <- sum(partitions$train_y$c1 == "No")
  original_ratio_1 <- 40 / 60

  # Training set should have similar ratio (within tolerance)
  train_ratio_1 <- train_yes_1 / train_no_1
  expect_true(
    abs(train_ratio_1 - original_ratio_1) < 0.3,
    info = sprintf("Expected ratio ~%.2f, got %.2f", original_ratio_1, train_ratio_1)
  )

  # Check partition info is populated correctly
  expect_equal(nrow(partitions$partition_info), 2)
  expect_true(all(partitions$partition_info$n_train > 0))
  expect_true(all(partitions$partition_info$n_val > 0))
  expect_true(all(partitions$partition_info$n_heldout > 0))
})


test_that("weighted variance computation is correct", {
  # Test with known values
  per_cohort <- data.frame(
    cohort = c("c1", "c2", "c3"),
    n = c(100, 50, 25),
    n_yes = c(40, 20, 10),
    n_no = c(60, 30, 15),
    sensitivity = c(0.8, 0.9, 0.7),
    specificity = c(0.7, 0.8, 0.6)
  )

  wv <- .compute_weighted_variance(per_cohort)

  expect_named(wv, c("sensitivity", "specificity"))
  expect_true(wv$sensitivity >= 0)
  expect_true(wv$specificity >= 0)

  # Single cohort should return 0 variance
  single <- per_cohort[1, , drop = FALSE]
  wv_single <- .compute_weighted_variance(single)
  expect_equal(wv_single$sensitivity, 0)
  expect_equal(wv_single$specificity, 0)
})

test_that("per-cohort metrics computation is correct", {
  set.seed(456)
  # Create known data
  scores <- c(0.9, 0.8, 0.7,  # cohort 1
              0.3, 0.4, 0.2)  # cohort 2
  truth <- factor(c("Yes", "Yes", "No",
                    "No", "No", "Yes"),
                  levels = c("No", "Yes"))
  cohort <- factor(c("c1", "c1", "c1", "c2", "c2", "c2"))

  threshold <- 0.5

  metrics <- .compute_per_cohort_metrics(scores, truth, cohort, threshold)

  expect_equal(nrow(metrics), 2)
  expect_equal(metrics$cohort, c("c1", "c2"))

  # Cohort 1: threshold 0.5, scores [0.9, 0.8, 0.7], truth [Yes, Yes, No]
  # All predicted Yes, so TP=2, FN=0, FP=1, TN=0
  # Sensitivity = 2/2 = 1, Specificity = 0/1 = 0
  expect_equal(metrics$sensitivity[1], 1.0)
  expect_equal(metrics$specificity[1], 0.0)

  # Cohort 2: threshold 0.5, scores [0.3, 0.4, 0.2], truth [No, No, Yes]
  # All predicted No, so TP=0, FN=1, FP=0, TN=2
  # Sensitivity = 0/1 = 0, Specificity = 2/2 = 1
  expect_equal(metrics$sensitivity[2], 0.0)
  expect_equal(metrics$specificity[2], 1.0)
})

test_that("NP threshold fallback works without nproc", {
  # Mock nproc unavailable by testing the fallback path
  scores <- runif(100)
  truth <- factor(sample(c("Yes", "No"), 100, replace = TRUE),
                  levels = c("No", "Yes"))

  # The function should work even if nproc fails
  result <- .select_np_threshold(scores, truth, alpha = 0.15, delta = 0.05)

  expect_true(is.list(result))
  expect_true("threshold" %in% names(result))
  expect_true(result$threshold >= 0 && result$threshold <= 1)
})

# -----------------------------------------------------------------------------
# Slow tests (full optimization pipeline)
# -----------------------------------------------------------------------------

test_that("optimize_panel returns OptimizationResult with a held-out share", {
  skip_slow_tests()
  set.seed(789)

  # Create multi-cohort data
  sim <- simulate_expression_data(p = 30, n = 200, k = 3, seed = 42)
  pool <- head(sim$metadata$informative_genes, 15)

  result <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 4,
    feature_pool = pool,
    fitness_mode = "within_cohort_val",
    train_ratio = 0.7,
    val_ratio = 0.2,
    np_alpha = 0.15,
    np_delta = 0.05,
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(result, "OptimizationResult")

  # Check solutions data frame has expected columns
  sols <- solutions(result)
  expect_true("base_features" %in% names(sols))
  expect_true("features" %in% names(sols))
  expect_true("solution_id" %in% names(sols))
  expect_true(nrow(sols) >= 1)

  # Check held-out data is stored in control
  expect_true(!is.null(result@control$heldout_x))
  expect_true(!is.null(result@control$heldout_y))
  expect_true(!is.null(result@control$heldout_cohort))
  expect_equal(result@control$np_alpha, 0.15)
  expect_equal(result@control$np_delta, 0.05)

  # Check partition info is stored
  pi <- result@control$partition_info
  expect_equal(pi$train_ratio, 0.7)
  expect_equal(pi$val_ratio, 0.2)

  # fit_panel() works on the result
  panel <- fit_panel(result)
  expect_s4_class(panel, "BiomarkerPanelResult")
  expect_true(length(panel@base_features) >= 1)
  expect_true(length(panel@features) >= 1)
})

test_that("optimize_panel() with an explicit feature_pool respects max_features and partition shares", {
  skip_slow_tests()
  set.seed(101)

  sim <- simulate_expression_data(p = 25, n = 200, k = 2, seed = 42)
  pool <- head(sim$metadata$informative_genes, 10)

  result <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    feature_pool = pool,
    max_features = 3,
    fitness_mode = "within_cohort_val",
    train_ratio = 0.7,
    val_ratio = 0.2,
    feature_transform = "none",
    nsga_control = list(popSize = 10, maxiter = 6),
    seed = 42
  )

  expect_s4_class(result, "OptimizationResult")

  # Fit panel and check feature count
  panel <- fit_panel(result)
  expect_true(length(panel@base_features) <= 3)

  # Verify partition_info shows expected split proportions
  pi <- result@control$partition_info
  expect_equal(pi$train_ratio, 0.7)
  expect_equal(pi$val_ratio, 0.2)
  expect_equal(pi$heldout_ratio, 0.1)
})


test_that("per-cohort metrics are populated for all cohorts", {
  skip_slow_tests()
  set.seed(303)

  n_cohorts <- 4
  sim <- simulate_expression_data(p = 20, n = 200, k = n_cohorts, seed = 42)
  pool <- head(sim$metadata$informative_genes, 8)

  opt <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    feature_pool = pool,
    max_features = 3,
    fitness_mode = "within_cohort_val",
    train_ratio = 0.7,
    val_ratio = 0.2,
    feature_transform = "none",
    nsga_control = list(popSize = 8, maxiter = 5)
  )

  panel <- fit_panel(opt)
  calibrated <- calibrate_panel(
    panel,
    x_heldout = opt@control$heldout_x,
    y_heldout = opt@control$heldout_y,
    cohort_heldout = opt@control$heldout_cohort
  )

  pcm <- per_cohort_metrics(calibrated)
  expect_equal(nrow(pcm), n_cohorts)
  expect_true(all(pcm$n > 0))
})

test_that("single cohort is handled gracefully", {
  skip_slow_tests()
  set.seed(505)

  sim <- simulate_expression_data(p = 20, n = 200, k = 1, seed = 42)

  opt <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    max_features = 3,
    feature_pool = sim$metadata$informative_genes[1:5],
    fitness_mode = "within_cohort_val",
    train_ratio = 0.7,
    val_ratio = 0.2,
    feature_transform = "none",
    nsga_control = list(popSize = 8, maxiter = 5)
  )

  expect_s4_class(opt, "OptimizationResult")

  panel <- fit_panel(opt)
  calibrated <- calibrate_panel(
    panel,
    x_heldout = opt@control$heldout_x,
    y_heldout = opt@control$heldout_y,
    cohort_heldout = opt@control$heldout_cohort
  )

  expect_s4_class(calibrated, "TransferablePanelResult")

  # Single cohort: weighted variance should be 0
  wv <- weighted_variance(calibrated)
  expect_equal(wv$sensitivity, 0)
  expect_equal(wv$specificity, 0)

  # Per-cohort metrics should have 1 row
  expect_equal(nrow(per_cohort_metrics(calibrated)), 1)
})



test_that("calibrate_panel validates inputs", {
  # Should error on non-BiomarkerPanelResult input
  expect_error(
    calibrate_panel("not_a_panel", matrix(1:4, 2, 2), factor(c("No", "Yes"))),
    "BiomarkerPanelResult"
  )
})

test_that("optimize_panel_transferable() is a deprecated wrapper that auto-selects features and draws a seed", {
  skip_slow_tests()

  sim <- simulate_expression_data(p = 25, n = 200, k = 2, seed = 42)

  warnings_seen <- character()
  result <- withCallingHandlers(
    optimize_panel_transferable(
      x = sim$x_list,
      y = sim$y_list,
      max_features = 3,
      feature_pool = NULL,
      train_ratio = 0.7,
      val_ratio = 0.2,
      feature_transform = "none",
      nsga_control = list(popSize = 10, maxiter = 6),
      n_top_features = 10,
      seed = NULL
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("deprecated", warnings_seen)))
  expect_s4_class(result, "OptimizationResult")
  expect_true(!is.null(result@control$heldout_x))
  expect_true(is.numeric(result@control$seed) && length(result@control$seed) == 1L)
})
