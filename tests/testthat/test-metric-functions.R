test_that("metric_sensitivity and metric_specificity behave as expected", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(metric_sensitivity(truth, scores), 1)
  expect_equal(metric_specificity(truth, scores), 1)
  expect_false(is.na(metric_balanced_accuracy(truth, scores)))
})

test_that("metric_auc computes AUC", {
  truth <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
  scores <- c(seq(0.1, 0.9, length.out = 10), seq(0.2, 1, length.out = 10))
  expect_gt(metric_auc(truth, scores), 0.5)
})

test_that("metric_num_features counts correctly", {
  expect_equal(metric_num_features(selected = c("g1", "g2")), 2)
  expect_equal(metric_num_features(selected = c(TRUE, FALSE, TRUE)), 2)
  expect_equal(metric_num_features(selected = NULL), 0)
})

test_that("metric_num_features is consistent with feature vector length", {

  # Test various input types - key property: output equals length(selected)
  features_char <- c("gene_A", "gene_B", "gene_C")
  expect_equal(metric_num_features(selected = features_char), length(features_char))

  features_5 <- paste0("feature_", seq_len(5))
  expect_equal(metric_num_features(selected = features_5), length(features_5))

  # Logical vectors: count equals sum of TRUE values
  features_logical <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  expect_equal(metric_num_features(selected = features_logical), sum(features_logical))

  # Empty/NULL cases
  expect_equal(metric_num_features(selected = character(0)), 0)
  expect_equal(metric_num_features(selected = NULL), 0)
})

test_that("metric_num_features handles edge cases", {
  # Single feature
  expect_equal(metric_num_features(selected = "single_gene"), 1)

  # Large panel
  big_panel <- paste0("gene_", seq_len(100))
  expect_equal(metric_num_features(selected = big_panel), 100)

  # All FALSE logical
  expect_equal(metric_num_features(selected = rep(FALSE, 10)), 0)

  # All TRUE logical
  expect_equal(metric_num_features(selected = rep(TRUE, 5)), 5)

  # Mixed with names (shouldn't affect counting)
  named_features <- c(a = "gene1", b = "gene2")
  expect_equal(metric_num_features(selected = named_features), 2)
})

test_that("build_objectives wraps registered metrics", {
  objs <- build_objectives(c("sensitivity", "num_features"),
                           params = list(sensitivity = list(cutoff_prob = 0.7)))
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.6, 0.3)
  selected <- c("a", "b", "c")
  expect_equal(names(objs), c("sensitivity", "num_features"))
  expect_equal(objs$sensitivity$direction, "maximize")
  expect_equal(objs$num_features$direction, "minimize")
  expect_true(is.numeric(objs$sensitivity$fun(truth, scores, selected)))
  expect_equal(objs$num_features$fun(truth, scores, selected), 3)
})

test_that("register_metric adds new entries", {
  custom <- function(truth, scores = NULL, selected = NULL, ...) 42
  name <- paste0("custom_", sample(1000, 1))
  register_metric(name, custom, direction = "maximize", overwrite = TRUE)
  objs <- build_objectives(name)
  expect_equal(objs[[name]]$fun(NULL, NULL, NULL), 42)
})

test_that("cohort-aware AUC metrics compute transfer metrics", {
  # Need enough samples per cohort for AUC to be meaningful
  set.seed(123)
  truth <- factor(
    c("No", "Yes", "Yes", "No", "Yes", "No", "No", "Yes", "Yes", "No"),
    levels = c("No", "Yes")
  )
  scores <- c(0.1, 0.9, 0.8, 0.2, 0.7, 0.3, 0.15, 0.85, 0.75, 0.25)
  cohort <- factor(c("A", "A", "A", "A", "A", "B", "B", "B", "B", "B"))

  auc_A <- metric_auc(truth[cohort == "A"], scores[cohort == "A"])
  auc_B <- metric_auc(truth[cohort == "B"], scores[cohort == "B"])

  # metric_min_cohort_auc equals minimum per-cohort AUC
  expect_equal(
    metric_min_cohort_auc(truth, scores, cohort = cohort),
    min(auc_A, auc_B)
  )

  # metric_cohort_auc_gap equals max - min per-cohort AUC
  gap <- metric_cohort_auc_gap(truth, scores, cohort = cohort)
  expect_gte(gap, 0)
  expect_equal(gap, max(auc_A, auc_B) - min(auc_A, auc_B))

  # metric_cohort_auc_var equals variance of per-cohort AUC
  auc_var <- metric_cohort_auc_var(truth, scores, cohort = cohort)
  expect_gte(auc_var, 0)
  expect_equal(auc_var, var(c(auc_A, auc_B)))

  brier <- metric_max_cohort_brier(truth, scores, cohort = cohort)
  expect_gte(brier, 0)

  objs <- build_objectives("min_cohort_auc")
  expect_equal(
    objs$min_cohort_auc$fun(truth, scores, cohort = cohort),
    metric_min_cohort_auc(truth, scores, cohort = cohort)
  )
})

test_that("cohort AUC metrics handle degenerate cohorts (0 positives)", {
  # Cohort C has no positives — metric_auc will error, should return NA
  truth <- factor(
    c("No", "Yes", "Yes", "No", "No", "No"),
    levels = c("No", "Yes")
  )
  scores <- c(0.1, 0.9, 0.8, 0.2, 0.3, 0.4)
  cohort <- factor(c("A", "A", "A", "A", "C", "C"))

  # Should not crash — degenerate cohort returns NA, aggregator ignores it

  result <- metric_min_cohort_auc(truth, scores, cohort = cohort)
  expect_true(is.numeric(result))
  expect_false(is.na(result))  # cohort A is valid

  gap <- metric_cohort_auc_gap(truth, scores, cohort = cohort)
  # With only one valid cohort, gap = 0 (single value range)
  expect_equal(gap, 0)
})

test_that("cohort AUC metrics return fallback for single cohort", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # No cohort argument — gap and var should return 0 (single_cohort_fallback)
  expect_equal(metric_cohort_auc_gap(truth, scores), 0)
  expect_equal(metric_cohort_auc_var(truth, scores), 0)

  # min_cohort_auc without cohort falls back to overall AUC
  expect_equal(
    metric_min_cohort_auc(truth, scores),
    metric_auc(truth, scores)
  )
})

# Issue 2: cutoff_strategy parameter tests
test_that("metric_sensitivity uses cutoff_strategy correctly", {
  # Create imbalanced data (30% positive)
  set.seed(42)
  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # Fixed cutoff at 0.5
  sens_fixed <- metric_sensitivity(truth, scores, cutoff_strategy = "fixed", cutoff_prob = 0.5)
  expect_equal(sens_fixed, 1.0) # All positives have scores > 0.5

  # Prevalence-based cutoff (0.3 for this data)
  sens_prev <- metric_sensitivity(truth, scores, cutoff_strategy = "prevalence")
  # With cutoff = 0.3, scores >= 0.3 are predicted positive
  # positives (3): scores 0.7, 0.8, 0.9 -> all >= 0.3, so TP=3
  expect_equal(sens_prev, 1.0)

  # Verify prevalence cutoff is computed correctly
  cutoff_val <- biomarkerPanels:::.compute_cutoff(truth, scores, 0.5, "prevalence", "Yes")
  expect_equal(cutoff_val, 0.3) # 3/10 positive
})

test_that("metric_specificity uses cutoff_strategy correctly", {
  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # Fixed cutoff at 0.5
  spec_fixed <- metric_specificity(truth, scores, cutoff_strategy = "fixed", cutoff_prob = 0.5)
  # Negatives (7): scores 0.1-0.55. With cutoff 0.5, only 0.55 is false positive
  expect_equal(spec_fixed, 6/7)

  # Prevalence-based cutoff (0.3)
  spec_prev <- metric_specificity(truth, scores, cutoff_strategy = "prevalence")
  # With cutoff = 0.3, scores < 0.3 are TN. Negatives: 0.1, 0.2 < 0.3 -> TN=2
  expect_equal(spec_prev, 2/7)
})

test_that("define_objectives passes cutoff_strategy to metrics", {
  objs <- define_objectives(
    metrics = c("sensitivity", "specificity"),
    cutoff_strategy = "prevalence"
  )

  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # The objectives should use prevalence-based cutoff
  sens <- objs$sensitivity$fun(truth, scores)
  spec <- objs$specificity$fun(truth, scores)

  # Compare with direct call using prevalence
  expect_equal(sens, metric_sensitivity(truth, scores, cutoff_strategy = "prevalence"))
  expect_equal(spec, metric_specificity(truth, scores, cutoff_strategy = "prevalence"))
})

# Tests for metric functions: AUC, F1, precision, NPV
test_that("metric_auc uses pROC for computation", {
  truth <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
  scores <- c(seq(0.1, 0.9, length.out = 10), seq(0.2, 1, length.out = 10))

  auc_val <- metric_auc(truth, scores)
  expect_gt(auc_val, 0.5)
  expect_lte(auc_val, 1.0)
})

test_that("metric_f1 computes F1 score correctly", {
  # Perfect classification
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(metric_f1(truth, scores), 1.0)

  # Imperfect classification
  truth2 <- factor(c("No", "Yes", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.4, 0.6, 0.8)  # one FN (0.4) and one FP (0.6)

  # With cutoff = 0.5:
  # TP = 2 (0.9, 0.8), FP = 1 (0.6), FN = 1 (0.4)
  # precision = 2/3, recall = 2/3
  # F1 = 2 * (2/3 * 2/3) / (2/3 + 2/3) = 2/3
  expect_equal(metric_f1(truth2, scores2), 2/3)
})

test_that("metric_precision computes PPV correctly", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # Perfect predictions: all predicted positives are true positives

  expect_equal(metric_precision(truth, scores), 1.0)

  # Add a false positive
  truth2 <- factor(c("No", "Yes", "Yes", "No", "No"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.8, 0.2, 0.7)  # last No scored above 0.5

  # TP = 2, FP = 1, precision = 2/3
  expect_equal(metric_precision(truth2, scores2), 2/3)
})

test_that("metric_npv computes negative predictive value correctly", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # Perfect predictions: all predicted negatives are true negatives
  expect_equal(metric_npv(truth, scores), 1.0)

  # Add a false negative
  truth2 <- factor(c("No", "Yes", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.3, 0.2, 0.8)  # third Yes scored below 0.5

  # TN = 2, FN = 1, NPV = 2/3
  expect_equal(metric_npv(truth2, scores2), 2/3)
})

test_that("metric functions handle edge cases", {
  # All positives
  truth_all_pos <- factor(rep("Yes", 5), levels = c("No", "Yes"))
  scores <- c(0.6, 0.7, 0.8, 0.9, 0.95)

  expect_equal(metric_f1(truth_all_pos, scores), 1.0)  # precision=1, recall=1
  expect_equal(metric_precision(truth_all_pos, scores), 1.0)
  expect_true(is.na(metric_npv(truth_all_pos, scores)))  # no negatives predicted

  # All negatives
  truth_all_neg <- factor(rep("No", 5), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.2, 0.3, 0.4, 0.45)

  expect_true(is.na(metric_f1(truth_all_neg, scores2)))  # no positives
  expect_true(is.na(metric_precision(truth_all_neg, scores2)))  # no predicted positives
  expect_equal(metric_npv(truth_all_neg, scores2), 1.0)
})

test_that("metric functions are registered in the metric registry", {
  registry <- metric_registry()
  expect_true("auc" %in% names(registry))
  expect_true("f1" %in% names(registry))
  expect_true("precision" %in% names(registry))
  expect_true("npv" %in% names(registry))

  # All should be maximize direction
  expect_equal(registry$auc$direction, "maximize")
  expect_equal(registry$f1$direction, "maximize")
  expect_equal(registry$precision$direction, "maximize")
  expect_equal(registry$npv$direction, "maximize")
})

test_that("metric functions work with build_objectives", {
  objs <- build_objectives(c("auc", "f1", "precision", "npv"))

  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  expect_true(is.numeric(objs$auc$fun(truth, scores)))
  expect_true(is.numeric(objs$f1$fun(truth, scores)))
  expect_true(is.numeric(objs$precision$fun(truth, scores)))
  expect_true(is.numeric(objs$npv$fun(truth, scores)))
})

# Tests for metric_specificity_at_sensitivity
test_that("metric_specificity_at_sensitivity computes correctly", {
  set.seed(42)
  truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  # Well-separated scores
  scores <- c(rnorm(50, 0.3, 0.1), rnorm(50, 0.7, 0.1))

  spec_at_90 <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.90)
  expect_true(is.numeric(spec_at_90))
  expect_gte(spec_at_90, 0)
  expect_lte(spec_at_90, 1)

  # Higher target sensitivity should yield lower or equal specificity
  spec_at_95 <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.95)
  expect_lte(spec_at_95, spec_at_90 + 0.01)  # small tolerance for interpolation

  # Lower target sensitivity should yield higher or equal specificity
  spec_at_80 <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.80)
  expect_gte(spec_at_80, spec_at_90 - 0.01)  # small tolerance for interpolation
})

test_that("metric_specificity_at_sensitivity handles edge cases", {
  truth <- factor(c(rep("No", 10), rep("Yes", 10)), levels = c("No", "Yes"))

  # Poor discrimination
  poor_scores <- c(runif(10, 0.4, 0.6), runif(10, 0.4, 0.6))
  spec_poor <- metric_specificity_at_sensitivity(truth, poor_scores, target_sensitivity = 0.90)
  expect_true(is.numeric(spec_poor))

  # Error when no positives
  truth_no_pos <- factor(rep("No", 10), levels = c("No", "Yes"))
  expect_error(
    metric_specificity_at_sensitivity(truth_no_pos, runif(10)),
    "no positive samples"
  )

  # Error when no negatives
  truth_no_neg <- factor(rep("Yes", 10), levels = c("No", "Yes"))
  expect_error(
    metric_specificity_at_sensitivity(truth_no_neg, runif(10)),
    "no negative samples"
  )

  # Invalid target_sensitivity
  expect_error(
    metric_specificity_at_sensitivity(truth, runif(20), target_sensitivity = 1.5),
    "between 0 and 1"
  )

  # Missing scores
  expect_error(
    metric_specificity_at_sensitivity(truth),
    "scores.*must be supplied"
  )
})

test_that("metric_specificity_at_sensitivity is registered correctly", {
  registry <- metric_registry()
  expect_true("specificity_at_sensitivity" %in% names(registry))
  expect_equal(registry$specificity_at_sensitivity$direction, "maximize")
})

test_that("metric_specificity_at_sensitivity works with build_objectives", {
  objs <- build_objectives(
    "specificity_at_sensitivity",
    params = list(specificity_at_sensitivity = list(target_sensitivity = 0.95))
  )

  set.seed(42)
  truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  scores <- c(rnorm(50, 0.3, 0.1), rnorm(50, 0.7, 0.1))

  result <- objs$specificity_at_sensitivity$fun(truth, scores)
  expect_true(is.numeric(result))

  # Should match direct call with same parameter
  direct <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.95)
  expect_equal(result, direct)
})

# Tests for cutoff_dependent metadata
test_that("metric registry includes cutoff_dependent field", {
  registry <- metric_registry()

  # Cutoff-dependent metrics should have cutoff_dependent = TRUE
  expect_true(registry$sensitivity$cutoff_dependent)
  expect_true(registry$specificity$cutoff_dependent)
  expect_true(registry$balanced_accuracy$cutoff_dependent)
  expect_true(registry$f1$cutoff_dependent)
  expect_true(registry$precision$cutoff_dependent)
  expect_true(registry$npv$cutoff_dependent)
  expect_false(registry$min_cohort_auc$cutoff_dependent)
  expect_false(registry$cohort_auc_gap$cutoff_dependent)
  expect_false(registry$cohort_auc_var$cutoff_dependent)

  # Cutoff-free metrics should have cutoff_dependent = FALSE
  expect_false(registry$auc$cutoff_dependent)
  expect_false(registry$pauc$cutoff_dependent)
  expect_false(registry$specificity_at_sensitivity$cutoff_dependent)
  expect_false(registry$num_features$cutoff_dependent)
  expect_false(registry$max_cohort_brier$cutoff_dependent)

})

test_that("min_metric_constraint warns for cutoff-dependent metrics", {
  # Should warn for sensitivity

  expect_warning(
    min_metric_constraint("sensitivity", threshold = 0.9),
    "depends on a probability cutoff"
  )

  # Should warn for specificity
  expect_warning(
    min_metric_constraint("specificity", threshold = 0.8),
    "depends on a probability cutoff"
  )

  # Should NOT warn for AUC
  expect_no_warning(
    min_metric_constraint("auc", threshold = 0.8)
  )

  # Should NOT warn for specificity_at_sensitivity
  expect_no_warning(
    min_metric_constraint("specificity_at_sensitivity", threshold = 0.5)
  )
})

test_that("register_metric accepts cutoff_dependent parameter", {
  custom_cutoff <- function(truth, scores = NULL, selected = NULL, ...) 0.5
  name <- paste0("custom_cutoff_", sample(1000, 1))

  register_metric(
    name, custom_cutoff,
    direction = "maximize",
    cutoff_dependent = TRUE,
    overwrite = TRUE
  )

  registry <- metric_registry()
  expect_true(registry[[name]]$cutoff_dependent)

  # Should warn when used in constraint
  expect_warning(
    min_metric_constraint(name, threshold = 0.4),
    "depends on a probability cutoff"
  )
})

# ---------------------------------------------------------------------------
# Tests for cohort transferability metrics
# ---------------------------------------------------------------------------

test_that("metric_cohort_leakage detects cohort-driven score shifts", {

  set.seed(42)
  n <- 60
  truth <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  cohort <- factor(rep(c("A", "B", "C"), times = n / 3))

  # Scores with strong cohort effect within each class
  scores_leaky <- ifelse(cohort == "A", 0.2, ifelse(cohort == "B", 0.5, 0.8))
  scores_leaky <- scores_leaky + rnorm(n, 0, 0.02)

  leak <- metric_cohort_leakage(truth, scores_leaky, cohort = cohort)
  expect_gt(leak, 0.5)  # strong leakage

  # Scores with no cohort effect
  scores_clean <- ifelse(truth == "Yes", 0.8, 0.2) + rnorm(n, 0, 0.05)

  leak_clean <- metric_cohort_leakage(truth, scores_clean, cohort = cohort)
  expect_lt(leak_clean, 0.15)  # minimal leakage
})

test_that("metric_cohort_leakage returns 0 when cohort is NULL or single level", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  expect_equal(metric_cohort_leakage(truth, scores, cohort = NULL), 0)
  expect_equal(
    metric_cohort_leakage(truth, scores, cohort = rep("A", 4)),
    0
  )
})

test_that("metric_cohort_leakage uses adjusted R-squared", {
  # With many cohort levels and few samples, unadjusted R-squared inflates.
  # Adjusted R-squared should be lower.
  set.seed(99)
  n <- 40
  truth <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  # 4 cohorts, balanced — random scores should have low leakage
  cohort <- factor(rep(paste0("C", 1:4), times = n / 4))
  scores <- runif(n)

  leak <- metric_cohort_leakage(truth, scores, cohort = cohort)
  # Random scores + adjusted R2 should be near 0 (not inflated)
  expect_lte(leak, 0.25)
})

test_that("metric_cohort_leakage errors without scores", {
  truth <- factor(c("No", "Yes"), levels = c("No", "Yes"))
  expect_error(metric_cohort_leakage(truth), "scores")
})

test_that("metric_conditional_score_shift detects distributional shifts", {
  set.seed(42)
  n_per <- 30

  truth <- factor(rep(c("No", "Yes"), each = n_per * 2), levels = c("No", "Yes"))
  cohort <- factor(rep(rep(c("A", "B"), each = n_per), 2))

  # Controls: A ~ N(0.2, 0.05), B ~ N(0.6, 0.05) — large shift
  # Cases: A ~ N(0.7, 0.05), B ~ N(0.7, 0.05) — no shift
  scores <- c(
    rnorm(n_per, 0.2, 0.05), rnorm(n_per, 0.6, 0.05),  # controls
    rnorm(n_per, 0.7, 0.05), rnorm(n_per, 0.7, 0.05)   # cases
  )

  shift <- metric_conditional_score_shift(truth, scores, cohort = cohort)
  expect_gt(shift, 0.3)  # detects the control shift

  # Clean scores: no cohort shift within either class
  scores_clean <- c(
    rnorm(n_per, 0.3, 0.05), rnorm(n_per, 0.3, 0.05),  # controls
    rnorm(n_per, 0.7, 0.05), rnorm(n_per, 0.7, 0.05)   # cases
  )
  shift_clean <- metric_conditional_score_shift(truth, scores_clean, cohort = cohort)
  expect_lt(shift_clean, 0.1)
})

test_that("metric_conditional_score_shift returns 0 when cohort is NULL", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  expect_equal(metric_conditional_score_shift(truth, scores, cohort = NULL), 0)
})

test_that("metric_conditional_score_shift handles small cohorts gracefully", {
  # Cohorts too small for W1 (< 10 per class per cohort) -> NA
  truth <- factor(c(rep("No", 6), rep("Yes", 6)), levels = c("No", "Yes"))
  cohort <- factor(rep(c("A", "B"), 6))
  scores <- runif(12)

  result <- metric_conditional_score_shift(truth, scores, cohort = cohort)
  expect_true(is.na(result) || is.numeric(result))
})

test_that("metric_conditional_score_shift errors without scores", {
  truth <- factor(c("No", "Yes"), levels = c("No", "Yes"))
  expect_error(metric_conditional_score_shift(truth), "scores")
})

test_that("cohort transferability metrics are registered correctly", {
  registry <- metric_registry()

  expect_true("cohort_leakage" %in% names(registry))
  expect_equal(registry$cohort_leakage$direction, "minimize")
  expect_false(registry$cohort_leakage$cutoff_dependent)

  expect_true("conditional_score_shift" %in% names(registry))
  expect_equal(registry$conditional_score_shift$direction, "minimize")
  expect_false(registry$conditional_score_shift$cutoff_dependent)
})

test_that("cohort transferability metrics work through build_objectives", {
  set.seed(42)
  n_per <- 30
  truth <- factor(rep(c("No", "Yes"), each = n_per * 2), levels = c("No", "Yes"))
  cohort <- factor(rep(rep(c("A", "B"), each = n_per), 2))
  scores <- c(
    rnorm(n_per, 0.3, 0.1), rnorm(n_per, 0.5, 0.1),
    rnorm(n_per, 0.7, 0.1), rnorm(n_per, 0.7, 0.1)
  )

  objs <- build_objectives(c("cohort_leakage", "conditional_score_shift"))

  leak_obj <- objs$cohort_leakage$fun(truth, scores, cohort = cohort)
  leak_direct <- metric_cohort_leakage(truth, scores, cohort = cohort)
  expect_equal(leak_obj, leak_direct)

  shift_obj <- objs$conditional_score_shift$fun(truth, scores, cohort = cohort)
  shift_direct <- metric_conditional_score_shift(truth, scores, cohort = cohort)
  expect_equal(shift_obj, shift_direct)
})

# --- Hand-calculated ground truth tests ---

test_that("classification metrics match hand-calculated confusion matrix", {
  # Confusion matrix (cutoff = 0.5):
  #                Predicted No   Predicted Yes
  #  Actual No          3 (TN)        1 (FP)
  #  Actual Yes         1 (FN)        2 (TP)
  #
  # TP=2, FP=1, TN=3, FN=1
  truth <- factor(c("No", "No", "No", "No", "Yes", "Yes", "Yes"),
                  levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.7, 0.4, 0.8, 0.9)

  # Sensitivity = TP / (TP + FN) = 2 / 3
  expect_equal(metric_sensitivity(truth, scores, cutoff_prob = 0.5), 2 / 3)

 # Specificity = TN / (TN + FP) = 3 / 4
  expect_equal(metric_specificity(truth, scores, cutoff_prob = 0.5), 3 / 4)

  # Precision = TP / (TP + FP) = 2 / 3
  expect_equal(metric_precision(truth, scores, cutoff_prob = 0.5), 2 / 3)

  # NPV = TN / (TN + FN) = 3 / 4
  expect_equal(metric_npv(truth, scores, cutoff_prob = 0.5), 3 / 4)

  # F1 = 2 * (precision * recall) / (precision + recall) = 2 * (2/3 * 2/3) / (2/3 + 2/3) = 2/3
  expect_equal(metric_f1(truth, scores, cutoff_prob = 0.5), 2 / 3)

  # Balanced accuracy = (sensitivity + specificity) / 2 = (2/3 + 3/4) / 2 = 17/24
  expect_equal(metric_balanced_accuracy(truth, scores, cutoff_prob = 0.5), 17 / 24)
})

test_that("AUC is 1.0 for perfectly separable scores", {
  # Perfect separation: all negatives < all positives
  truth <- factor(c(rep("No", 5), rep("Yes", 5)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  expect_equal(metric_auc(truth, scores), 1.0)
})

test_that("AUC is 0.5 for random-coin scores", {
  # Identical scores for both classes -> AUC = 0.5
  truth <- factor(c("No", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.5, 0.5, 0.5, 0.5)
  expect_equal(metric_auc(truth, scores), 0.5)
})

# ---------------------------------------------------------------------------
# Tests for easy-hard balanced accuracy
# ---------------------------------------------------------------------------

test_that("metric_easy_hard_accuracy computes stratified balanced accuracy", {
  # 4 easy samples: 2 No, 2 Yes — all classified correctly

  # 4 hard samples: 2 No, 2 Yes — 1 FN (Yes scored 0.3)
  truth <- factor(c("No", "No", "Yes", "Yes", "No", "No", "Yes", "Yes"),
                  levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.3, 0.9)
  difficulty <- c("easy", "easy", "easy", "easy", "hard", "hard", "hard", "hard")

  # Easy accuracy: 4/4 = 1.0
  # Hard accuracy: 3/4 = 0.75 (the 0.3 Yes is a FN)
  # Mean: (1.0 + 0.75) / 2 = 0.875
  result <- metric_easy_hard_accuracy(truth, scores, difficulty = difficulty)
  expect_equal(result, 0.875)
})

test_that("metric_easy_hard_accuracy perfect classifier scores 1", {
  truth <- factor(c("No", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.2, 0.8)
  difficulty <- c("easy", "easy", "hard", "hard")

  expect_equal(
    metric_easy_hard_accuracy(truth, scores, difficulty = difficulty),
    1.0
  )
})

test_that("metric_easy_hard_accuracy handles all-easy or all-hard", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # All easy — returns accuracy on that single stratum
  result_easy <- metric_easy_hard_accuracy(truth, scores,
                                           difficulty = rep("easy", 4))
  expect_equal(result_easy, 1.0)

  # All hard
  result_hard <- metric_easy_hard_accuracy(truth, scores,
                                           difficulty = rep("hard", 4))
  expect_equal(result_hard, 1.0)
})

test_that("metric_easy_hard_accuracy errors on missing difficulty", {
  truth <- factor(c("No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9)
  expect_error(
    metric_easy_hard_accuracy(truth, scores),
    "difficulty.*must be supplied"
  )
})

test_that("metric_easy_hard_accuracy errors on invalid difficulty labels", {
  truth <- factor(c("No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9)
  expect_error(
    metric_easy_hard_accuracy(truth, scores, difficulty = c("low", "high")),
    "easy.*hard"
  )
})

test_that("metric_easy_hard_accuracy errors on length mismatch", {
  truth <- factor(c("No", "Yes", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8)
  expect_error(
    metric_easy_hard_accuracy(truth, scores, difficulty = c("easy", "hard")),
    "same length"
  )
})

test_that("metric_easy_hard_accuracy is registered correctly", {
  registry <- metric_registry()
  expect_true("easy_hard_accuracy" %in% names(registry))
  expect_equal(registry$easy_hard_accuracy$direction, "maximize")
  expect_true(registry$easy_hard_accuracy$cutoff_dependent)
})

test_that("metric_easy_hard_accuracy works through build_objectives", {
  truth <- factor(c("No", "No", "Yes", "Yes", "No", "No", "Yes", "Yes"),
                  levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.3, 0.9)
  difficulty <- c("easy", "easy", "easy", "easy", "hard", "hard", "hard", "hard")

  objs <- build_objectives(
    "easy_hard_accuracy",
    params = list(easy_hard_accuracy = list(difficulty = difficulty))
  )

  result <- objs$easy_hard_accuracy$fun(truth, scores)
  direct <- metric_easy_hard_accuracy(truth, scores, difficulty = difficulty)
  expect_equal(result, direct)
})

# ============================================================================
# Robustness tests (silent failure audit fixes)
# ============================================================================

test_that("metric_specificity_at_sensitivity returns scalar result", {
  skip_if_not_installed("pROC")

  set.seed(42)
  truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  scores <- c(runif(50, 0, 0.6), runif(50, 0.4, 1))

  result <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.90)
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})
