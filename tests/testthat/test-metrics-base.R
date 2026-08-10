# Tests for foundational single-cutoff or continuous classification metrics (e.g., sensitivity, specificity, AUC).
test_that("metric_sensitivity and metric_specificity behave as expected", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(metric_sensitivity(truth, scores), 1)
  expect_equal(metric_specificity(truth, scores), 1)
  expect_false(is.na(metric_balanced_accuracy(truth, scores)))
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

# Tests for metric_sensitivity_at_specificity
test_that("metric_sensitivity_at_specificity computes correctly", {
  set.seed(42)
  truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  # Well-separated scores
  scores <- c(rnorm(50, 0.3, 0.1), rnorm(50, 0.7, 0.1))

  sens_at_90 <- metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.90)
  expect_true(is.numeric(sens_at_90))
  expect_length(sens_at_90, 1L)
  expect_gte(sens_at_90, 0)
  expect_lte(sens_at_90, 1)

  # Higher target specificity should yield lower or equal sensitivity
  sens_at_95 <- metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.95)
  expect_lte(sens_at_95, sens_at_90 + 0.01)  # small tolerance for interpolation

  # Lower target specificity should yield higher or equal sensitivity
  sens_at_80 <- metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.80)
  expect_gte(sens_at_80, sens_at_90 - 0.01)  # small tolerance for interpolation
})

test_that("metric_sensitivity_at_specificity is the transpose of its sibling", {
  set.seed(7)
  truth <- factor(c(rep("No", 60), rep("Yes", 60)), levels = c("No", "Yes"))
  scores <- c(rnorm(60, 0.3, 0.15), rnorm(60, 0.7, 0.15))

  # Round-trip: the specificity reported at a target sensitivity should recover
  # (at least) that sensitivity when fed back in as the specificity target.
  spec <- metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.90)
  sens <- metric_sensitivity_at_specificity(truth, scores, target_specificity = spec)
  expect_gte(sens, 0.90 - 0.02)

  # A perfectly separated classifier hits 1 on both axes.
  perfect_truth <- factor(c(rep("No", 20), rep("Yes", 20)), levels = c("No", "Yes"))
  perfect_scores <- c(runif(20, 0, 0.4), runif(20, 0.6, 1))
  expect_equal(
    metric_sensitivity_at_specificity(perfect_truth, perfect_scores,
                                      target_specificity = 0.95),
    1
  )
})

test_that("metric_sensitivity_at_specificity handles edge cases", {
  truth <- factor(c(rep("No", 10), rep("Yes", 10)), levels = c("No", "Yes"))

  # Poor discrimination
  poor_scores <- c(runif(10, 0.4, 0.6), runif(10, 0.4, 0.6))
  sens_poor <- metric_sensitivity_at_specificity(truth, poor_scores, target_specificity = 0.90)
  expect_true(is.numeric(sens_poor))
  expect_false(is.na(sens_poor))

  # Degenerate single-class inputs return the worst value (0), not NA, so the
  # NSGA fitness path never sees a missing objective value.
  truth_no_pos <- factor(rep("No", 10), levels = c("No", "Yes"))
  expect_equal(metric_sensitivity_at_specificity(truth_no_pos, runif(10)), 0)

  truth_no_neg <- factor(rep("Yes", 10), levels = c("No", "Yes"))
  expect_equal(metric_sensitivity_at_specificity(truth_no_neg, runif(10)), 0)

  # Constant scores are degenerate but still return a finite value
  const <- metric_sensitivity_at_specificity(truth, rep(0.5, 20))
  expect_true(is.finite(const))

  # Invalid target_specificity
  expect_error(
    metric_sensitivity_at_specificity(truth, runif(20), target_specificity = 1.5),
    "between 0 and 1"
  )

  # Missing scores
  expect_error(
    metric_sensitivity_at_specificity(truth),
    "scores.*must be supplied"
  )
})

test_that("sensitivity_at_specificity is registered and usable as an objective", {
  registry <- metric_registry()
  expect_true("sensitivity_at_specificity" %in% names(registry))
  expect_equal(registry$sensitivity_at_specificity$direction, "maximize")
  expect_false(isTRUE(registry$sensitivity_at_specificity$cutoff_dependent))

  set.seed(11)
  truth <- factor(c(rep("No", 30), rep("Yes", 30)), levels = c("No", "Yes"))
  scores <- c(rnorm(30, 0.3, 0.1), rnorm(30, 0.7, 0.1))

  objs <- build_objectives(
    metrics = c("sensitivity_at_specificity", "num_features"),
    params = list(sensitivity_at_specificity = list(target_specificity = 0.95))
  )
  value <- objs$sensitivity_at_specificity$fun(truth, scores, selected = c("A", "B"))
  expect_equal(
    value,
    metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.95)
  )
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
