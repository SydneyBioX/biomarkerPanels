# Tests for Rcpp implementations - verify equivalence with pure R and speedup

test_that(".compute_roc_curve_cpp matches pure R implementation", {
  set.seed(42)
  n <- 1000
  truth <- factor(sample(c("No", "Yes"), n, replace = TRUE), levels = c("No", "Yes"))
  scores <- runif(n)

  # Both implementations should produce identical results

  roc_r <- biomarkerPanels:::.compute_roc_curve_pure_r(truth, scores, "Yes")
  roc_cpp <- biomarkerPanels:::.compute_roc_curve(truth, scores, "Yes")

  expect_equal(nrow(roc_cpp), nrow(roc_r))
  expect_equal(roc_cpp$threshold, roc_r$threshold, tolerance = 1e-10)
  expect_equal(roc_cpp$tpr, roc_r$tpr, tolerance = 1e-10)
  expect_equal(roc_cpp$fpr, roc_r$fpr, tolerance = 1e-10)
  expect_equal(roc_cpp$sensitivity, roc_r$sensitivity, tolerance = 1e-10)
  expect_equal(roc_cpp$specificity, roc_r$specificity, tolerance = 1e-10)
})

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

test_that(".compute_roc_curve_cpp is faster than pure R", {
  skip_if_not_installed("microbenchmark")
  skip_on_cran()

  set.seed(42)
  n <- 10000
  truth <- factor(sample(c("No", "Yes"), n, replace = TRUE), levels = c("No", "Yes"))
  scores <- runif(n)

  bm <- microbenchmark::microbenchmark(
    pure_r = biomarkerPanels:::.compute_roc_curve_pure_r(truth, scores, "Yes"),
    rcpp = biomarkerPanels:::.compute_roc_curve(truth, scores, "Yes"),
    times = 10
  )

  median_r <- median(bm$time[bm$expr == "pure_r"])
  median_cpp <- median(bm$time[bm$expr == "rcpp"])
  speedup <- median_r / median_cpp

  message(sprintf("ROC curve speedup: %.1fx (R: %.1fms, Rcpp: %.1fms)",
                  speedup, median_r / 1e6, median_cpp / 1e6))

  # Expect at least 2x faster
  expect_gt(speedup, 2)
})

# Feature selection function equivalence tests

test_that(".select_stable_genes_cpp matches pure R implementation", {
  set.seed(123)
  n_genes <- 500
  n_cohorts <- 4

  t_matrix <- matrix(rnorm(n_genes * n_cohorts, mean = 2, sd = 1),
                     nrow = n_genes, ncol = n_cohorts)
  rownames(t_matrix) <- paste0("gene_", seq_len(n_genes))
  se_matrix <- matrix(abs(rnorm(n_genes * n_cohorts, mean = 0.5, sd = 0.2)),
                      nrow = n_genes, ncol = n_cohorts)

  for (method in c("precision_weighted", "cv_t_stats", "inverse_t_se")) {
    result_r <- biomarkerPanels:::.select_stable_genes_pure_r(
      t_matrix, se_matrix, method, top_n = 50
    )
    result_cpp <- biomarkerPanels:::.select_stable_genes(
      t_matrix, se_matrix, method, top_n = 50
    )
    expect_equal(result_cpp, result_r,
                 info = sprintf("Method: %s", method))
  }
})

test_that(".score_transferable_features_cpp matches pure R implementation", {
  set.seed(456)
  n_features <- 100
  n_cohorts <- 5

  coefficient_matrix <- matrix(rnorm(n_features * n_cohorts),
                               nrow = n_features, ncol = n_cohorts)
  rownames(coefficient_matrix) <- paste0("feature_", seq_len(n_features))

  result_r <- biomarkerPanels:::.score_transferable_features_pure_r(
    coefficient_matrix,
    min_coefficient = 0.01,
    require_sign_consistency = TRUE,
    sign_consistency_threshold = 0.8
  )
  result_cpp <- biomarkerPanels:::.score_transferable_features(
    coefficient_matrix,
    min_coefficient = 0.01,
    require_sign_consistency = TRUE,
    sign_consistency_threshold = 0.8
  )

  expect_equal(result_cpp$feature, result_r$feature)
  expect_equal(result_cpp$mean_abs, result_r$mean_abs, tolerance = 1e-10)
  expect_equal(result_cpp$sd, result_r$sd, tolerance = 1e-10)
  expect_equal(result_cpp$sign_agreement, result_r$sign_agreement, tolerance = 1e-10)
})

test_that(".aggregate_de_pvalues_cpp matches pure R implementation", {
  set.seed(789)
  n_genes <- 1000
  n_cohorts <- 4

  t_matrix <- matrix(rnorm(n_genes * n_cohorts, mean = 0, sd = 2),
                     nrow = n_genes, ncol = n_cohorts)
  rownames(t_matrix) <- paste0("gene_", seq_len(n_genes))

  for (method in c("Stouffer", "Fisher")) {
    result_r <- biomarkerPanels:::.aggregate_de_pvalues_pure_r(t_matrix, method)
    result_cpp <- biomarkerPanels:::.aggregate_de_pvalues(t_matrix, method)

    expect_equal(names(result_cpp), names(result_r))
    expect_equal(result_cpp, result_r, tolerance = 1e-6,
                 info = sprintf("Method: %s", method))
  }
})
