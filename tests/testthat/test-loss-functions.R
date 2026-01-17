test_that("loss_sensitivity and loss_specificity behave as expected", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(loss_sensitivity(truth, scores), 1)
  expect_equal(loss_specificity(truth, scores), 1)
  expect_false(is.na(loss_balanced_accuracy(truth, scores)))
})

test_that("loss_auc computes AUC", {
  truth <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
  scores <- c(seq(0.1, 0.9, length.out = 10), seq(0.2, 1, length.out = 10))
  expect_gt(loss_auc(truth, scores), 0.5)
})

test_that("loss_num_features counts correctly", {
  expect_equal(loss_num_features(selected = c("g1", "g2")), 2)
  expect_equal(loss_num_features(selected = c(TRUE, FALSE, TRUE)), 2)
  expect_equal(loss_num_features(selected = NULL), 0)
})

test_that("build_objectives wraps registered losses", {
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

test_that("register_loss_function adds new entries", {
  custom <- function(truth, scores = NULL, selected = NULL, ...) 42
  name <- paste0("custom_", sample(1000, 1))
  register_loss_function(name, custom, direction = "maximize", overwrite = TRUE)
  objs <- build_objectives(name)
  expect_equal(objs[[name]]$fun(NULL, NULL, NULL), 42)
})

test_that("cohort-aware losses compute transfer metrics", {
  truth <- factor(
    c("No", "Yes", "Yes", "No", "Yes", "No"),
    levels = c("No", "Yes")
  )
  scores <- c(0.1, 0.9, 0.8, 0.2, 0.6, 0.4)
  cohort <- factor(c("A", "A", "B", "B", "B", "A"))
  expect_equal(
    loss_min_cohort_sensitivity(truth, scores, cohort = cohort),
    min(
      loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
      loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
    )
  )
  expect_equal(
    loss_min_cohort_specificity(truth, scores, cohort = cohort),
    min(
      loss_specificity(truth[cohort == "A"], scores[cohort == "A"]),
      loss_specificity(truth[cohort == "B"], scores[cohort == "B"])
    )
  )
  gap <- loss_cohort_sensitivity_gap(truth, scores, cohort = cohort)
  expect_gte(gap, 0)
  expect_equal(gap, max(
    loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
    loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
  ) - min(
    loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
    loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
  ))

  brier <- loss_max_cohort_brier(truth, scores, cohort = cohort)
  expect_gte(brier, 0)

  x <- matrix(
    c(1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 7),
    nrow = 6, ncol = 2, byrow = FALSE
  )
  shift <- loss_max_cohort_mean_shift(truth, scores, cohort = cohort, x = x)
  expect_gte(shift, 0)

  objs <- build_objectives(c("min_cohort_sensitivity", "max_cohort_mean_shift"))
  expect_equal(
    objs$min_cohort_sensitivity$fun(truth, scores, cohort = cohort),
    loss_min_cohort_sensitivity(truth, scores, cohort = cohort)
  )
  expect_equal(
    objs$max_cohort_mean_shift$fun(truth, scores, cohort = cohort, x = x),
    loss_max_cohort_mean_shift(truth, scores, cohort = cohort, x = x)
  )
})

# Issue 2: cutoff_strategy parameter tests
test_that("loss_sensitivity uses cutoff_strategy correctly", {
  # Create imbalanced data (30% positive)
  set.seed(42)
  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # Fixed cutoff at 0.5
  sens_fixed <- loss_sensitivity(truth, scores, cutoff_strategy = "fixed", cutoff_prob = 0.5)
  expect_equal(sens_fixed, 1.0) # All positives have scores > 0.5

  # Prevalence-based cutoff (0.3 for this data)
  sens_prev <- loss_sensitivity(truth, scores, cutoff_strategy = "prevalence")
  # With cutoff = 0.3, scores >= 0.3 are predicted positive
  # positives (3): scores 0.7, 0.8, 0.9 -> all >= 0.3, so TP=3
  expect_equal(sens_prev, 1.0)

  # Verify prevalence cutoff is computed correctly
  cutoff_val <- biomarkerPanels:::.compute_cutoff(truth, scores, 0.5, "prevalence", "Yes")
  expect_equal(cutoff_val, 0.3) # 3/10 positive
})

test_that("loss_specificity uses cutoff_strategy correctly", {
  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # Fixed cutoff at 0.5
  spec_fixed <- loss_specificity(truth, scores, cutoff_strategy = "fixed", cutoff_prob = 0.5)
  # Negatives (7): scores 0.1-0.55. With cutoff 0.5, only 0.55 is false positive
  expect_equal(spec_fixed, 6/7)

  # Prevalence-based cutoff (0.3)
  spec_prev <- loss_specificity(truth, scores, cutoff_strategy = "prevalence")
  # With cutoff = 0.3, scores < 0.3 are TN. Negatives: 0.1, 0.2 < 0.3 -> TN=2
  expect_equal(spec_prev, 2/7)
})

test_that("define_objectives passes cutoff_strategy to losses", {
  objs <- define_objectives(
    losses = c("sensitivity", "specificity"),
    cutoff_strategy = "prevalence"
  )

  truth <- factor(c(rep("No", 7), rep("Yes", 3)), levels = c("No", "Yes"))
  scores <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.8, 0.9)

  # The objectives should use prevalence-based cutoff
  sens <- objs$sensitivity$fun(truth, scores)
  spec <- objs$specificity$fun(truth, scores)

  # Compare with direct call using prevalence
  expect_equal(sens, loss_sensitivity(truth, scores, cutoff_strategy = "prevalence"))
  expect_equal(spec, loss_specificity(truth, scores, cutoff_strategy = "prevalence"))
})

# Tests for new loss functions: AUC, F1, precision, NPV
test_that("loss_auc uses pROC for computation", {
  truth <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
  scores <- c(seq(0.1, 0.9, length.out = 10), seq(0.2, 1, length.out = 10))

  auc_val <- loss_auc(truth, scores)
  expect_gt(auc_val, 0.5)
  expect_lte(auc_val, 1.0)
})

test_that("loss_f1 computes F1 score correctly", {
  # Perfect classification
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(loss_f1(truth, scores), 1.0)

  # Imperfect classification
  truth2 <- factor(c("No", "Yes", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.4, 0.6, 0.8)  # one FN (0.4) and one FP (0.6)

  # With cutoff = 0.5:
  # TP = 2 (0.9, 0.8), FP = 1 (0.6), FN = 1 (0.4)
  # precision = 2/3, recall = 2/3
  # F1 = 2 * (2/3 * 2/3) / (2/3 + 2/3) = 2/3
  expect_equal(loss_f1(truth2, scores2), 2/3)
})

test_that("loss_precision computes PPV correctly", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # Perfect predictions: all predicted positives are true positives

  expect_equal(loss_precision(truth, scores), 1.0)

  # Add a false positive
  truth2 <- factor(c("No", "Yes", "Yes", "No", "No"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.8, 0.2, 0.7)  # last No scored above 0.5

  # TP = 2, FP = 1, precision = 2/3
  expect_equal(loss_precision(truth2, scores2), 2/3)
})

test_that("loss_npv computes negative predictive value correctly", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  # Perfect predictions: all predicted negatives are true negatives
  expect_equal(loss_npv(truth, scores), 1.0)

  # Add a false negative
  truth2 <- factor(c("No", "Yes", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.9, 0.3, 0.2, 0.8)  # third Yes scored below 0.5

  # TN = 2, FN = 1, NPV = 2/3
  expect_equal(loss_npv(truth2, scores2), 2/3)
})

test_that("new loss functions handle edge cases", {
  # All positives
  truth_all_pos <- factor(rep("Yes", 5), levels = c("No", "Yes"))
  scores <- c(0.6, 0.7, 0.8, 0.9, 0.95)

  expect_equal(loss_f1(truth_all_pos, scores), 1.0)  # precision=1, recall=1
  expect_equal(loss_precision(truth_all_pos, scores), 1.0)
  expect_true(is.na(loss_npv(truth_all_pos, scores)))  # no negatives predicted

  # All negatives
  truth_all_neg <- factor(rep("No", 5), levels = c("No", "Yes"))
  scores2 <- c(0.1, 0.2, 0.3, 0.4, 0.45)

  expect_true(is.na(loss_f1(truth_all_neg, scores2)))  # no positives
  expect_true(is.na(loss_precision(truth_all_neg, scores2)))  # no predicted positives
  expect_equal(loss_npv(truth_all_neg, scores2), 1.0)
})

test_that("new loss functions are registered in the loss registry", {
  registry <- loss_registry()
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

test_that("new loss functions work with build_objectives", {
  objs <- build_objectives(c("auc", "f1", "precision", "npv"))

  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)

  expect_true(is.numeric(objs$auc$fun(truth, scores)))
  expect_true(is.numeric(objs$f1$fun(truth, scores)))
  expect_true(is.numeric(objs$precision$fun(truth, scores)))
  expect_true(is.numeric(objs$npv$fun(truth, scores)))
})

# Tests for loss_specificity_at_sensitivity
test_that("loss_specificity_at_sensitivity computes correctly", {
  set.seed(42)
  truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
  # Well-separated scores
  scores <- c(rnorm(50, 0.3, 0.1), rnorm(50, 0.7, 0.1))

  spec_at_90 <- loss_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.90)
  expect_true(is.numeric(spec_at_90))
  expect_gte(spec_at_90, 0)
  expect_lte(spec_at_90, 1)

  # Higher target sensitivity should yield lower or equal specificity
  spec_at_95 <- loss_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.95)
  expect_lte(spec_at_95, spec_at_90 + 0.01)  # small tolerance for interpolation

  # Lower target sensitivity should yield higher or equal specificity
  spec_at_80 <- loss_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.80)
  expect_gte(spec_at_80, spec_at_90 - 0.01)  # small tolerance for interpolation
})

test_that("loss_specificity_at_sensitivity handles edge cases", {
  truth <- factor(c(rep("No", 10), rep("Yes", 10)), levels = c("No", "Yes"))

  # Poor discrimination
  poor_scores <- c(runif(10, 0.4, 0.6), runif(10, 0.4, 0.6))
  spec_poor <- loss_specificity_at_sensitivity(truth, poor_scores, target_sensitivity = 0.90)
  expect_true(is.numeric(spec_poor))

  # Error when no positives
  truth_no_pos <- factor(rep("No", 10), levels = c("No", "Yes"))
  expect_error(
    loss_specificity_at_sensitivity(truth_no_pos, runif(10)),
    "no positive samples"
  )

  # Error when no negatives
  truth_no_neg <- factor(rep("Yes", 10), levels = c("No", "Yes"))
  expect_error(
    loss_specificity_at_sensitivity(truth_no_neg, runif(10)),
    "no negative samples"
  )

  # Invalid target_sensitivity
  expect_error(
    loss_specificity_at_sensitivity(truth, runif(20), target_sensitivity = 1.5),
    "between 0 and 1"
  )

  # Missing scores
  expect_error(
    loss_specificity_at_sensitivity(truth),
    "scores.*must be supplied"
  )
})

test_that("loss_specificity_at_sensitivity is registered correctly", {
  registry <- loss_registry()
  expect_true("specificity_at_sensitivity" %in% names(registry))
  expect_equal(registry$specificity_at_sensitivity$direction, "maximize")
})

test_that("loss_specificity_at_sensitivity works with build_objectives", {
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
  direct <- loss_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.95)
  expect_equal(result, direct)
})

# Tests for cutoff_dependent metadata
test_that("loss registry includes cutoff_dependent field", {
  registry <- loss_registry()

  # Cutoff-dependent metrics should have cutoff_dependent = TRUE
  expect_true(registry$sensitivity$cutoff_dependent)
  expect_true(registry$specificity$cutoff_dependent)
  expect_true(registry$balanced_accuracy$cutoff_dependent)
  expect_true(registry$f1$cutoff_dependent)
  expect_true(registry$precision$cutoff_dependent)
  expect_true(registry$npv$cutoff_dependent)
  expect_true(registry$min_cohort_sensitivity$cutoff_dependent)
  expect_true(registry$min_cohort_specificity$cutoff_dependent)
  expect_true(registry$cohort_sensitivity_gap$cutoff_dependent)

  # Cutoff-free metrics should have cutoff_dependent = FALSE
  expect_false(registry$auc$cutoff_dependent)
  expect_false(registry$pauc$cutoff_dependent)
  expect_false(registry$specificity_at_sensitivity$cutoff_dependent)
  expect_false(registry$num_features$cutoff_dependent)
  expect_false(registry$max_cohort_brier$cutoff_dependent)
  expect_false(registry$max_cohort_mean_shift$cutoff_dependent)
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

test_that("register_loss_function accepts cutoff_dependent parameter", {
  custom_cutoff <- function(truth, scores = NULL, selected = NULL, ...) 0.5
  name <- paste0("custom_cutoff_", sample(1000, 1))

  register_loss_function(
    name, custom_cutoff,
    direction = "maximize",
    cutoff_dependent = TRUE,
    overwrite = TRUE
  )

  registry <- loss_registry()
  expect_true(registry[[name]]$cutoff_dependent)

  # Should warn when used in constraint
  expect_warning(
    min_metric_constraint(name, threshold = 0.4),
    "depends on a probability cutoff"
  )
})
