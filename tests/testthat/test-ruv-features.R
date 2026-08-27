# Test fixture: small synthetic data with known signal/batch/null structure
make_ruv_fixture <- function() {
  make_cohort <- function(seed, batch_shift) {
    set.seed(seed)
    n <- 30L
    p_signal <- 5L
    p_batch <- 5L
    p_null <- 10L
    p <- p_signal + p_batch + p_null

    x <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
    colnames(x) <- c(
      paste0("signal_", seq_len(p_signal)),
      paste0("batch_", seq_len(p_batch)),
      paste0("null_", seq_len(p_null))
    )
    y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

    # Signal genes: strong outcome association
    for (j in seq_len(p_signal)) {
      x[, j] <- x[, j] + rep(c(-1, 1), each = n / 2) * (1.5 + 0.2 * j)
    }

    # Batch genes: cohort-specific shift, no outcome association
    for (j in seq_len(p_batch)) {
      x[, p_signal + j] <- x[, p_signal + j] + batch_shift * (0.8 + 0.1 * j)
    }

    list(x = x, y = y)
  }

  cohort1 <- make_cohort(101, -1.5)
  cohort2 <- make_cohort(202, 0.0)
  cohort3 <- make_cohort(303, 1.5)

  list(
    x_list = list(cohort1$x, cohort2$x, cohort3$x),
    y_list = list(cohort1$y, cohort2$y, cohort3$y),
    signal_genes = paste0("signal_", 1:5),
    batch_genes = paste0("batch_", 1:5),
    null_genes = paste0("null_", 1:10)
  )
}

test_that("select_ruv_features returns expected structure", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2
  ))

  expect_type(result, "list")
  expect_named(result, c("features", "scores", "ruv_summary", "corrected",
    "settings"))
  expect_type(result$features, "character")
  expect_true(length(result$features) <= 10)

  # Scores data.frame structure
  expect_s3_class(result$scores, "data.frame")
  expect_true(all(c("feature", "discrimination", "alpha_norm", "score", "rank")
    %in% names(result$scores)))

  # RUV summary
  expect_type(result$ruv_summary, "list")
  expect_equal(result$ruv_summary$k, 2L)
  expect_true(result$ruv_summary$n_controls > 0)

  # Settings
  expect_equal(result$settings$k, 2L)
  expect_equal(result$settings$scoring, "ratio")

  # Corrected should be NULL by default
  expect_null(result$corrected)
})

test_that("select_ruv_features prioritises signal genes over batch genes", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 20, k = 2
  ))

  # At least some signal genes should appear in top features
  signal_in_top <- intersect(result$features, fixture$signal_genes)
  expect_true(length(signal_in_top) >= 1)

  # Mean score of signal genes should exceed mean score of batch genes
  scores <- result$scores
  signal_scores <- scores$score[scores$feature %in% fixture$signal_genes]
  batch_scores <- scores$score[scores$feature %in% fixture$batch_genes]
  if (length(signal_scores) > 0 && length(batch_scores) > 0) {
    expect_gt(mean(signal_scores), mean(batch_scores))
  }
})

test_that("neg_control_method = 'list' works with user-supplied genes", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2,
    neg_control_method = "list",
    neg_control_genes = fixture$null_genes
  ))

  expect_type(result$features, "character")
  expect_true(length(result$features) >= 1)
  expect_true(all(result$ruv_summary$control_genes %in% fixture$null_genes))
})

test_that("k auto-selection works", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = NULL
  ))

  expect_true(result$ruv_summary$k >= 1L)
  expect_true(result$ruv_summary$k <= 10L)
  expect_type(result$features, "character")
})

test_that("iterate_controls refines without error", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2,
    iterate_controls = TRUE, max_iterations = 2
  ))

  expect_type(result$features, "character")
  expect_true(result$settings$iterate_controls)
})

test_that("both scoring methods produce valid output", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()

  result_ratio <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2, scoring = "ratio"
  ))
  result_z <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2, scoring = "composite_z"
  ))

  expect_equal(result_ratio$settings$scoring, "ratio")
  expect_equal(result_z$settings$scoring, "composite_z")

  # Both produce valid scores
  expect_true(all(is.finite(result_ratio$scores$score)))
  expect_true(all(is.finite(result_z$scores$score)))
})

test_that("return_corrected produces matrix of correct dimensions", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list, fixture$y_list,
    n_features = 10, k = 2, return_corrected = TRUE
  ))

  expect_true(is.matrix(result$corrected))
  # 3 cohorts x 30 samples each = 90 samples, 20 genes
  expect_equal(nrow(result$corrected), 90)
  expect_equal(ncol(result$corrected), 20)
})

test_that("single cohort works without list wrapping", {
  skip_if_not_installed("ruv")
  skip_slow_tests()

  fixture <- make_ruv_fixture()
  result <- suppressWarnings(select_ruv_features(
    fixture$x_list[[1]], fixture$y_list[[1]],
    n_features = 5, k = 1
  ))

  expect_type(result$features, "character")
  expect_true(length(result$features) >= 1)
})