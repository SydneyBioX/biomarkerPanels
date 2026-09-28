make_ratio_fixture <- function() {
  make_cohort <- function(shift, seed) {
    n <- 8L
    p <- 4L
    set.seed(seed)
    x <- matrix(rnorm(n * p, sd = 0.05), nrow = n, ncol = p)
    colnames(x) <- paste0("gene", seq_len(p))
    y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
    x[, 1] <- x[, 1] + rep(c(0, shift), each = n / 2)
    x[, 2] <- x[, 2] + 0.2
    x[, 3] <- x[, 3] + rep(c(0, shift / 4), each = n / 2)
    baseline <- seq(-0.2, 0.2, length.out = n)
    x[, 4] <- x[, 4] + baseline
    list(x = x, y = y)
  }

  list(
    cohorts = list(
      make_cohort(2.5, 1),
      make_cohort(2.0, 2)
    )
  )
}

test_that("select_de_features ranks shared signals", {
  fixture <- make_ratio_fixture()

  selected <- select_de_features(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_features = 3,
    combination_method = "Stouffer"
  )

  expect_type(selected, "character")
  expect_lte(length(selected), 3)
  expect_true("gene1" %in% selected)
})

test_that("select_features_for_ratios returns stable and informative sets", {
  fixture <- make_ratio_fixture()

  ratios <- select_features_for_ratios(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_stable = 2,
    n_informative = 3,
    stability_method = "precision_weighted",
    combination_method = "Stouffer"
  )

  expect_named(ratios, c("stable", "informative"))
  expect_length(ratios$stable, 2)
  expect_length(ratios$informative, 3)
  expect_true("gene2" %in% ratios$stable)
  expect_true("gene1" %in% ratios$informative)

  alt <- select_features_for_ratios(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_stable = 2,
    n_informative = 2,
    stability_method = "inverse_t_se"
  )
  expect_length(alt$stable, 2)
  expect_length(alt$informative, 2)

  expect_error(
    select_features_for_ratios(
      lapply(fixture$cohorts, `[[`, "x"),
      fixture$cohorts[[1]]$y,
      n_stable = 2,
      n_informative = 2
    ),
    "`x_list` and `y_list` must have the same length"
  )

  expect_error(
    select_features_for_ratios(
      lapply(fixture$cohorts, `[[`, "x"),
      lapply(fixture$cohorts, `[[`, "y"),
      stability_method = "unknown-method"
    ),
    "should be one of"
  )
})

test_that("select_transferable_features prioritises consistent ridge coefficients", {
  skip_slow_tests()
  make_cohort <- function(seed, delta) {
    set.seed(seed)
    n <- 40L
    x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(x) <- c("geneA", "geneB", "geneC")
    linear <- 1.5 * x[, "geneA"] - 0.5 * x[, "geneB"] + delta
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohort1 <- make_cohort(101, -0.2)
  cohort2 <- make_cohort(202, 0.2)

  set.seed(999)
  transferable <- select_transferable_features(
    list(cohort1$x, cohort2$x),
    list(cohort1$y, cohort2$y),
    n_features = 2,
    lambda_choice = "lambda_min"
  )

  expect_type(transferable$features, "character")
  expect_true("geneA" %in% transferable$features)
  expect_true(all(transferable$scores$sign_consistent))
  expect_true(all(transferable$scores$min_abs >= 0))
  expect_equal(colnames(transferable$coefficients), c("cohort_01", "cohort_02"))

  expect_error(
    select_transferable_features(
      list(cohort1$x, cohort2$x),
      cohort1$y
    ),
    "`x_list` and `y_list` must have the same length"
  )
})

test_that("select_transferable_features intersects cohort feature sets", {
  skip_slow_tests()
  make_dataset <- function(seed, cols) {
    set.seed(seed)
    n <- 40L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 1.2 * x[, "gene_common1"] - 0.8 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  ds1 <- make_dataset(111, c("gene_common1", "gene_common2", "gene_unique1"))
  ds2 <- make_dataset(222, c("gene_unique2", "gene_common2", "gene_common1"))

  result <- select_transferable_features(
    list(ds1$x, ds2$x),
    list(ds1$y, ds2$y),
    n_features = 2,
    lambda = 0.1,
    min_coefficient = 0.01,
    standardize = FALSE
  )

  expect_true(length(result$features) <= 2)
  expect_true(all(result$features %in% c("gene_common1", "gene_common2")))
  expect_equal(colnames(result$coefficients), c("cohort_01", "cohort_02"))
  expect_true(all(rownames(result$coefficients) %in% c("gene_common1", "gene_common2")))
})

# Issue 3: sign_consistency_threshold tests
test_that("sign_consistency_threshold allows partial sign consistency", {
  skip_slow_tests()
  # Create 4 cohorts where geneB has sign inconsistency in one cohort
  make_cohort <- function(seed, b_sign) {
    set.seed(seed)
    n <- 80L
    x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
    colnames(x) <- c("geneA", "geneB")
    # geneA always has positive effect, geneB has varying sign
    linear <- 1.5 * x[, "geneA"] + b_sign * 1.5 * x[, "geneB"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  # cohort 4 has opposite sign for geneB
  cohorts <- list(
    make_cohort(101, 1),   # positive geneB
    make_cohort(102, 1),   # positive geneB
    make_cohort(103, 1),   # positive geneB
    make_cohort(104, -1)   # negative geneB
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  # With 100% sign consistency (default), geneB should be excluded
  result_strict <- select_transferable_features(
    x_list, y_list,
    n_features = 2,
    require_sign_consistency = TRUE,
    sign_consistency_threshold = 1.0,
    lambda = 0.1
  )

  # With 75% sign consistency, geneB should be included (3/4 = 75% agree)
  result_lenient <- select_transferable_features(
    x_list, y_list,
    n_features = 2,
    require_sign_consistency = TRUE,
    sign_consistency_threshold = 0.75,
    lambda = 0.1
  )

  # geneA should be in both (100% sign consistency)
  expect_true("geneA" %in% result_strict$features)
  expect_true("geneA" %in% result_lenient$features)

  # geneB should be excluded from strict, included in lenient
  expect_false("geneB" %in% result_strict$features)
  expect_true("geneB" %in% result_lenient$features)

  # Verify sign_consistency_threshold is stored in settings
  expect_equal(result_strict$settings$sign_consistency_threshold, 1.0)
  expect_equal(result_lenient$settings$sign_consistency_threshold, 0.75)

  # Verify sign_agreement column exists and has proper values
  expect_true("sign_agreement" %in% names(result_lenient$scores))
})

# ==============================================================================
# select_discriminative_features tests
# ==============================================================================

test_that("select_discriminative_features ranks signal gene above offset gene", {
  fixture <- make_ratio_fixture()
  x_list <- lapply(fixture$cohorts, `[[`, "x")
  y_list <- lapply(fixture$cohorts, `[[`, "y")

  result <- select_discriminative_features(
    x_list, y_list,
    n_features = 4,
    min_auc = 0.5
  )

  expect_type(result$features, "character")
  expect_true(length(result$features) >= 1)
  expect_true("gene1" %in% result$features)

  # gene1 (strong signal) should rank above gene2 (cohort offset, weaker signal)
  scores <- result$scores
  gene1_score <- scores$score[scores$feature == "gene1"]
  gene2_score <- scores$score[scores$feature == "gene2"]
  expect_gt(gene1_score, gene2_score)

  # gene1 has different shifts per cohort (2.5 vs 2.0), so its cohort_disc
  # should be higher than gene4 (no cohort-dependent shift)
  gene1_cd <- scores$cohort_disc[scores$feature == "gene1"]
  gene4_cd <- scores$cohort_disc[scores$feature == "gene4"]
  expect_gt(gene1_cd, gene4_cd)

  # Return structure matches expected pattern
  expect_named(result, c("features", "scores", "per_cohort_auc", "settings"))
  expect_true(all(c("feature", "mean_auc", "sd_auc", "min_auc",
    "cohort_disc", "score") %in% names(result$scores)))
  expect_equal(ncol(result$per_cohort_auc), 2)
})

test_that("select_discriminative_features respects min_auc filter", {
  fixture <- make_ratio_fixture()
  x_list <- lapply(fixture$cohorts, `[[`, "x")
  y_list <- lapply(fixture$cohorts, `[[`, "y")

  result_strict <- select_discriminative_features(
    x_list, y_list,
    n_features = 4,
    min_auc = 0.95
  )

  result_lenient <- select_discriminative_features(
    x_list, y_list,
    n_features = 4,
    min_auc = 0.5
  )

  expect_true(length(result_strict$features) <= length(result_lenient$features))
  if (nrow(result_strict$scores) > 0) {
    expect_true(all(result_strict$scores$mean_auc >= 0.95))
  }
})

test_that("select_discriminative_features works with single cohort", {
  fixture <- make_ratio_fixture()
  x <- fixture$cohorts[[1]]$x
  y <- fixture$cohorts[[1]]$y

  result <- select_discriminative_features(x, y, n_features = 3, min_auc = 0.5)

  expect_type(result$features, "character")
  # With single cohort, cohort_disc should be 0 and sd_auc should be 0
  expect_true(all(result$scores$cohort_disc == 0))
  expect_true(all(result$scores$sd_auc == 0))
})

test_that("select_discriminative_features validates parameters", {
  fixture <- make_ratio_fixture()
  x_list <- lapply(fixture$cohorts, `[[`, "x")
  y_list <- lapply(fixture$cohorts, `[[`, "y")

  expect_error(
    select_discriminative_features(x_list, y_list, lambda_cohort = "bad"),
    "`lambda_cohort` must be a numeric scalar"
  )
  expect_error(
    select_discriminative_features(x_list, y_list, lambda_sd = c(1, 2)),
    "`lambda_sd` must be a numeric scalar"
  )
  expect_error(
    select_discriminative_features(x_list, y_list, min_auc = 1.5),
    "`min_auc` must be a numeric scalar"
  )
})

# ==============================================================================
# select_ruleout_features tests
# ==============================================================================

# Builds a matrix contrasting two kinds of case/control separation. Controls
# are N(0, 1). "auc_high" has 30 cases far above controls and 10 cases only
# modestly above them (mean 1.5), giving a high AUC but a poor specificity once
# the threshold must reach into that overlapping group to capture 90% of
# cases. "tail_clean" has 36 cases well above every control and 4 cases far
# below them: a lower AUC, but the top 90% of cases separate perfectly.
# "down" mirrors "tail_clean" in the opposite direction; "moderate" is a
# uniformly modest shift (specificity ~0.4 at 90% sensitivity); "noise" has no
# signal.
make_ruleout_fixture <- function(seed = 1L) {
  set.seed(seed)
  n_case <- 40L
  n_ctrl <- 40L
  ctrl <- rnorm(n_ctrl)

  auc_high <- c(rnorm(30, mean = 6), rnorm(10, mean = 1.5))
  tail_clean <- c(rnorm(36, mean = 5, sd = 0.3), rnorm(4, mean = -6))

  x <- cbind(
    auc_high = c(auc_high, ctrl),
    tail_clean = c(tail_clean, ctrl),
    down = c(-tail_clean, -ctrl),
    moderate = c(rnorm(n_case, mean = 1), rnorm(n_ctrl)),
    noise = rnorm(n_case + n_ctrl)
  )
  y <- factor(rep(c("Yes", "No"), c(n_case, n_ctrl)), levels = c("No", "Yes"))
  list(x = x, y = y)
}

test_that("select_ruleout_features ranks by tail specificity rather than AUC", {
  fx <- make_ruleout_fixture()

  # Sanity: auc_high really has the higher AUC.
  auc <- biomarkerPanels:::.compute_rank_auc(fx$x, as.integer(fx$y) - 1L)
  expect_gt(auc[["auc_high"]], auc[["tail_clean"]])

  selected <- select_ruleout_features(fx$x, fx$y, n_features = 5)

  expect_type(selected, "character")
  expect_equal(length(selected), 5)
  expect_lt(match("tail_clean", selected), match("auc_high", selected))
  expect_equal(selected[4:5], c("moderate", "noise"))
  expect_equal(length(select_ruleout_features(fx$x, fx$y, n_features = 2)), 2)
})

test_that("select_ruleout_features is direction-agnostic", {
  fx <- make_ruleout_fixture()
  spec <- biomarkerPanels:::.spec_at_sensitivity_by_feature(
    fx$x, as.integer(fx$y) - 1L, 0.90
  )
  expect_equal(spec[["down", "down"]], spec[["tail_clean", "up"]])
  expect_equal(spec[["tail_clean", "up"]], 1)
})

test_that("select_ruleout_features requires one direction across cohorts", {
  fx1 <- make_ruleout_fixture(1L)
  fx2 <- make_ruleout_fixture(2L)
  # "flip" separates perfectly in both cohorts, but its cases are high in
  # cohort 1 and low in cohort 2, so no single direction works in both.
  fx1$x <- cbind(fx1$x, flip = fx1$x[, "tail_clean"])
  fx2$x <- cbind(fx2$x, flip = fx2$x[, "down"])

  selected <- select_ruleout_features(
    list(fx1$x, fx2$x), list(fx1$y, fx2$y), n_features = 6
  )

  # Scoring each cohort's best direction before taking the worst cohort would
  # give "flip" a perfect score; a shared direction leaves it below the
  # consistently moderate feature.
  expect_lt(match("moderate", selected), match("flip", selected))
})

test_that("select_ruleout_features threshold always captures >= target sensitivity", {
  fx <- make_ruleout_fixture()
  y <- as.integer(fx$y) - 1L
  x_pos <- fx$x[y == 1L, ]
  n <- nrow(x_pos)
  for (target in c(0.8, 0.9, 0.95, 1)) {
    thr <- apply(x_pos, 2, biomarkerPanels:::.ruleout_thresholds, target = target)
    captured_up <- colMeans(sweep(x_pos, 2, thr[1L, ], ">="))
    captured_down <- colMeans(sweep(x_pos, 2, thr[2L, ], "<="))
    expect_true(all(captured_up >= target), info = paste("up", target))
    expect_true(all(captured_down >= target), info = paste("down", target))
    # No slack: exactly ceiling(n * target) cases are kept (no ties in fixture).
    expect_equal(unname(captured_up * n), rep(ceiling(n * target), ncol(x_pos)))
    expect_equal(unname(captured_down * n), rep(ceiling(n * target), ncol(x_pos)))
  }
})

test_that("select_ruleout_features 'min' aggregation penalises cohort-specific features", {
  fx1 <- make_ruleout_fixture(1L)
  fx2 <- make_ruleout_fixture(2L)
  # In cohort 2, tail_clean carries no signal at all.
  fx2$x[, "tail_clean"] <- fx2$x[, "noise"]

  x_list <- list(fx1$x, fx2$x)
  y_list <- list(fx1$y, fx2$y)

  by_min <- select_ruleout_features(x_list, y_list, n_features = 6)
  by_mean <- select_ruleout_features(x_list, y_list,
    n_features = 6,
    cohort_aggregation = "mean"
  )

  # "down" is perfect in both cohorts and tops both rankings.
  expect_equal(by_min[1], "down")
  expect_equal(by_mean[1], "down")
  # tail_clean: spec 1 in cohort 1, ~0 in cohort 2. moderate: ~0.4 in both.
  # Worst-cohort ranking prefers the consistently moderate feature; mean
  # ranking is seduced by the single excellent cohort.
  expect_lt(match("moderate", by_min), match("tail_clean", by_min))
  expect_lt(match("tail_clean", by_mean), match("moderate", by_mean))
})

test_that("select_ruleout_features handles NA values and rejects single-class cohorts", {
  fx <- make_ruleout_fixture()
  x_na <- fx$x
  # Blank three of tail_clean's clean cases: with 37 cases left only three may
  # be missed, so one of its four far-below cases must now be captured and
  # tail_clean's specificity drops from 1 to 0.
  x_na[1:3, "tail_clean"] <- NA
  x_na[41:43, "auc_high"] <- NA
  x_na[, "noise"] <- NA

  spec <- biomarkerPanels:::.spec_at_sensitivity_by_feature(
    x_na, as.integer(fx$y) - 1L, 0.90
  )
  expect_equal(max(spec["tail_clean", ]), 0)
  expect_equal(spec[["down", "down"]], 1)
  expect_true(all(is.na(spec["noise", ])))

  selected <- select_ruleout_features(x_na, fx$y, n_features = 10)
  expect_equal(length(selected), 4)
  expect_equal(selected[1], "down")

  y_one_class <- factor(rep("Yes", nrow(fx$x)), levels = c("No", "Yes"))
  expect_error(
    select_ruleout_features(fx$x, y_one_class),
    "at least one case and one control"
  )
  expect_error(select_ruleout_features(fx$x, fx$y, target_sensitivity = 0))
})
