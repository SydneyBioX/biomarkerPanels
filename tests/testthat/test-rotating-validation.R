# Tests for .generate_stratified_splits() and the rotating-validation
# fitness path of optimize_panel_transferable().

test_that(".generate_stratified_splits produces the requested number of splits", {
  set.seed(42)
  y <- factor(rep(c("Yes", "No"), each = 30), levels = c("No", "Yes"))
  cohort <- factor(rep(c("A", "B"), times = 30))

  splits <- .generate_stratified_splits(y, cohort, n_splits = 5L,
                                        val_ratio = 0.2, seed = 1L)
  expect_length(splits, 5L)
  for (s in splits) {
    expect_named(s, c("train", "val"))
    expect_true(length(intersect(s$train, s$val)) == 0L)
    expect_setequal(c(s$train, s$val), seq_along(y))
  }
})

test_that(".generate_stratified_splits preserves (cohort, class) balance", {
  set.seed(7)
  y <- factor(c(rep("Yes", 40), rep("No", 60)), levels = c("No", "Yes"))
  cohort <- factor(c(rep("A", 50), rep("B", 50)))

  splits <- .generate_stratified_splits(y, cohort, n_splits = 4L,
                                        val_ratio = 0.25, seed = 1L)

  strata <- interaction(cohort, y, drop = TRUE)
  strata_counts <- table(strata)

  for (s in splits) {
    val_counts <- table(strata[s$val])
    # Every stratum should appear in val for reasonably sized strata
    for (lv in names(strata_counts)) {
      if (strata_counts[[lv]] >= 4L) {
        expect_true(val_counts[[lv]] >= 1L,
                    info = sprintf("Stratum %s missing from val", lv))
      }
    }
    # Val proportion within each stratum should be close to requested ratio
    for (lv in names(strata_counts)) {
      if (strata_counts[[lv]] >= 4L) {
        prop <- val_counts[[lv]] / strata_counts[[lv]]
        expect_true(abs(prop - 0.25) < 0.3,
                    info = sprintf("Stratum %s val proportion %.2f", lv, prop))
      }
    }
  }
})

test_that(".generate_stratified_splits is deterministic given a seed", {
  y <- factor(rep(c("Yes", "No"), each = 20))
  cohort <- factor(rep(c("A", "B"), times = 20))

  s1 <- .generate_stratified_splits(y, cohort, n_splits = 3L,
                                    val_ratio = 0.3, seed = 123L)
  s2 <- .generate_stratified_splits(y, cohort, n_splits = 3L,
                                    val_ratio = 0.3, seed = 123L)
  expect_identical(s1, s2)
})

test_that(".generate_stratified_splits rotates different samples to val", {
  # Across K splits, a healthy rotation should expose most samples as val
  # at some point.
  y <- factor(rep(c("Yes", "No"), each = 50))
  cohort <- factor(rep(c("A", "B"), times = 50))
  splits <- .generate_stratified_splits(y, cohort, n_splits = 10L,
                                        val_ratio = 0.3, seed = 42L)
  ever_val <- unique(unlist(lapply(splits, `[[`, "val")))
  # Expect the majority of rows to land in val across 10 splits at ratio 0.3
  expect_true(length(ever_val) / length(y) > 0.8)
})

test_that(".generate_stratified_splits rejects invalid inputs", {
  y <- factor(c("Yes", "No", "Yes", "No"))
  cohort <- factor(c("A", "A", "B", "B"))
  expect_error(.generate_stratified_splits(y, cohort, 1L, 0.2),
               "n_splits.*>= 2")
  expect_error(.generate_stratified_splits(y, cohort, 3L, 0),
               "val_ratio.*between 0 and 1")
  expect_error(.generate_stratified_splits(y, cohort[1], 3L, 0.2),
               "matching length")
})

# Integration: rotating fitness runs end-to-end on a small fixture
# (slow — uses NSGA). Kept brief.

test_that("optimize_panel_transferable runs with within_cohort_rotating", {
  skip_slow_tests()
  set.seed(11)
  sim <- simulate_expression_data(p = 40L, n = 40L, k = 2L, seed = 11L,
                                  shift_scale = 0.3)
  # Use only the most informative genes as pool to keep the test fast
  pool <- head(sim$metadata$informative_genes, 10L)

  result <- optimize_panel_transferable(
    x = sim$x_list, y = sim$y_list,
    feature_pool = pool,
    max_features = 3L,
    train_ratio = 0.6,
    val_ratio = 0.25,
    fitness_mode = "within_cohort_rotating",
    n_val_splits = 4L,
    objectives = define_objectives(
      metrics = c("auc", "num_features")
    ),
    nsga_control = list(popSize = 12L, maxiter = 3L),
    seed = 11L,
    feature_transform = "none"
  )

  expect_s4_class(result, "OptimizationResult")
  expect_gt(n_solutions(result), 0L)
  expect_identical(result@control$fitness_mode, "within_cohort_rotating")
  expect_identical(result@control$n_val_splits, 4L)
})
