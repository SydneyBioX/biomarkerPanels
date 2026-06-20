# Tests for the Rcpp ROC implementation - production behavior and pROC agreement.
# (Pure-R reference implementations and their equivalence/benchmark tests were
# removed once the Rcpp versions were validated.)

test_that(".compute_roc_curve handles edge cases correctly", {
  # All positive - FPR should be NA (no negatives to compute FPR)
  # TPR should range from 0 to 1 as threshold decreases
  truth_all_pos <- factor(rep("Yes", 50), levels = c("No", "Yes"))
  scores <- runif(50)
  roc <- biomarkerPanels:::.compute_roc_curve(truth_all_pos, scores, "Yes")
  expect_true(all(is.na(roc$fpr)))
  expect_true(all(roc$tpr >= 0 & roc$tpr <= 1))
  # Last row (lowest threshold) should have TPR = 1
  expect_equal(max(roc$tpr), 1)

  # All negative - TPR should be NA (no positives to compute TPR)
  truth_all_neg <- factor(rep("No", 50), levels = c("No", "Yes"))
  roc <- biomarkerPanels:::.compute_roc_curve(truth_all_neg, scores, "Yes")
  expect_true(all(is.na(roc$tpr)))
  expect_true(all(roc$fpr >= 0 & roc$fpr <= 1 | is.na(roc$fpr)))

  # Tied scores - should still produce valid ROC
  set.seed(123)
  truth <- factor(sample(c("No", "Yes"), 100, replace = TRUE), levels = c("No", "Yes"))
  scores_tied <- sample(1:10 / 10, 100, replace = TRUE)
  roc <- biomarkerPanels:::.compute_roc_curve(truth, scores_tied, "Yes")
  expect_s3_class(roc, "data.frame")
  expect_true(nrow(roc) > 0)

  # Single sample
  truth_single <- factor("Yes", levels = c("No", "Yes"))
  scores_single <- 0.5
  roc <- biomarkerPanels:::.compute_roc_curve(truth_single, scores_single, "Yes")
  expect_s3_class(roc, "data.frame")
})

test_that(".compute_roc_curve_cpp matches pROC", {
  skip_if_not_installed("pROC")

  set.seed(42)
  n <- 500
  truth <- factor(sample(c("No", "Yes"), n, replace = TRUE), levels = c("No", "Yes"))
  scores <- runif(n)

  # pROC reference
  roc_proc <- pROC::roc(truth, scores, levels = c("No", "Yes"), direction = "<", quiet = TRUE)

  # Our implementation
  roc_ours <- biomarkerPanels:::.compute_roc_curve(truth, scores, "Yes")

  # Compare key operating points
  # pROC sensitivities should be present (or very close) in our curve
  sens_distances <- sapply(roc_proc$sensitivities, function(sens) {
    min(abs(roc_ours$sensitivity - sens), na.rm = TRUE)
  })
  expect_true(all(sens_distances < 1e-6), 
              info = "Some pROC sensitivities were not found in our curve")

  # Same for specificities
  spec_distances <- sapply(roc_proc$specificities, function(spec) {
    min(abs(roc_ours$specificity - spec), na.rm = TRUE)
  })
  expect_true(all(spec_distances < 1e-6),
              info = "Some pROC specificities were not found in our curve")
})
