# Tests for multi-cohort performance metrics like minimum cohort AUC and transferability scores.
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

