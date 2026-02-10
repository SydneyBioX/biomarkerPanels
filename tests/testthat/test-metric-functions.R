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
