#' Threshold Calibration and Per-Cohort Metrics
#'
#' Internal functions for Neyman-Pearson threshold selection and computing
#' per-cohort performance metrics for transferability evaluation.
#'
#' @name threshold_calibration
#' @keywords internal
NULL

#' Select Neyman-Pearson Classification Threshold
#'
#' Uses the nproc package to select a classification threshold that controls
#' the Type I error rate (1 - specificity) at level alpha.
#'
#' @param scores Numeric vector of predicted probabilities.
#' @param truth Binary response factor.
#' @param alpha Target Type I error rate (default 0.15).
#' @param delta Tolerance parameter (default 0.05).
#' @return List with threshold value and metadata.
#' @keywords internal
.select_np_threshold <- function(scores, truth, alpha = 0.15, delta = 0.05) {
  if (!requireNamespace("nproc", quietly = TRUE)) {
    stop(
      "Package 'nproc' is required for NP threshold selection. ",
      "Install it via install.packages('nproc').",
      call. = FALSE
    )
  }

  scores <- as.numeric(scores)
  y_binary <- as.integer(truth) - 1L # nproc expects 0/1 labels
  neg_scores <- scores[y_binary == 0L]

  empirical_fallback <- function(reason) {
    warning(
      reason,
      " Falling back to an empirical training-data quantile threshold, which ",
      "does NOT carry the (alpha, delta) Neyman-Pearson guarantee.",
      call. = FALSE
    )
    thr <- if (length(neg_scores) > 0L) {
      as.numeric(stats::quantile(neg_scores, probs = 1 - alpha, na.rm = TRUE))
    } else {
      0.5
    }
    list(threshold = thr, method = "empirical_quantile_fallback",
         alpha = alpha, delta = delta)
  }

  if (length(unique(y_binary)) < 2L) {
    return(empirical_fallback("NP threshold selection needs both classes present."))
  }

  # Genuine Neyman-Pearson threshold via nproc::npc, which controls the type-I
  # error (1 - specificity) at <= alpha with probability >= 1 - delta. We fit on
  # the 1-D score with method = "lda" (penlog/glmnet rejects a single predictor)
  # and then recover the threshold ON THE SCORE SCALE from npc's own
  # NP-controlled predictions. npc's fits[[1]]$cutoff lives on lda's internal
  # projection scale and is NOT comparable to the probability scores we threshold
  # against, so it cannot be used directly. Recovering min(score | predicted Yes)
  # reproduces npc's NP decision when applied as `score >= threshold`.
  np_err <- NULL
  np_threshold <- tryCatch({
    np_model <- nproc::npc(
      x = as.matrix(scores),
      y = y_binary,
      method = "lda",
      alpha = alpha,
      delta = delta,
      split = 1L
    )
    pred <- predict(np_model, as.matrix(scores))
    pos_scores <- scores[pred$pred.label == 1]
    if (length(pos_scores)) min(pos_scores) else NA_real_
  }, error = function(e) {
    np_err <<- conditionMessage(e)
    NA_real_
  })

  # Verify the recovered threshold actually controls the type-I error on the
  # negatives. This single check catches sign flips, boundary ties, npc errors
  # and degeneracy (no predicted positives); only if it passes is the threshold a
  # genuine NP threshold.
  realized_fpr <- if (is.finite(np_threshold) && length(neg_scores) > 0L) {
    mean(neg_scores >= np_threshold)
  } else {
    NA_real_
  }

  if (is.finite(np_threshold) && !is.na(realized_fpr) &&
      realized_fpr <= alpha + 1e-8) {
    list(threshold = as.numeric(np_threshold), method = "nproc",
         alpha = alpha, delta = delta)
  } else {
    reason <- if (!is.null(np_err)) {
      paste0("nproc::npc failed (", np_err, ").")
    } else {
      paste0("nproc did not yield a threshold controlling the type-I error at ",
             "alpha = ", alpha, " on this data.")
    }
    empirical_fallback(reason)
  }
}

#' Compute Per-Cohort Performance Metrics
#'
#' Calculates sensitivity and specificity for each cohort at the given
#' classification threshold.
#'
#' @param scores Numeric vector of predicted probabilities.
#' @param truth Binary response factor.
#' @param cohort Factor indicating cohort membership.
#' @param threshold Classification threshold.
#' @return Data frame with columns: cohort, n, n_yes, n_no, sensitivity, specificity.
#' @keywords internal
.compute_per_cohort_metrics <- function(scores, truth, cohort, threshold) {
  cohort_levels <- levels(cohort)
  results <- lapply(cohort_levels, function(coh) {
    idx <- which(cohort == coh)
    coh_scores <- scores[idx]
    coh_truth <- truth[idx]

    n <- length(idx)
    n_yes <- sum(coh_truth == "Yes")
    n_no <- sum(coh_truth == "No")

    sens <- metric_sensitivity(coh_truth, coh_scores, cutoff_prob = threshold)
    spec <- metric_specificity(coh_truth, coh_scores, cutoff_prob = threshold)

    data.frame(
      cohort = coh,
      n = n,
      n_yes = n_yes,
      n_no = n_no,
      sensitivity = sens,
      specificity = spec,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Compute Inverse Sample-Size Weighted Variance
#'
#' Calculates weighted variance of sensitivity and specificity across cohorts,
#' using inverse sample size as weights (smaller cohorts contribute less).
#'
#' @param per_cohort_metrics Data frame from .compute_per_cohort_metrics().
#' @return List with sensitivity and specificity variance values.
#' @keywords internal
.compute_weighted_variance <- function(per_cohort_metrics) {
  n_cohorts <- nrow(per_cohort_metrics)

  if (n_cohorts <= 1) {
    return(list(sensitivity = 0, specificity = 0))
  }

  # Inverse sample-size weights (normalized)
  weights <- 1 / per_cohort_metrics$n
  weights <- weights / sum(weights)

  compute_weighted_var <- function(values, w) {
    valid <- !is.na(values)
    if (sum(valid) <= 1) {
      return(0)
    }

    v <- values[valid]
    wt <- w[valid]
    wt <- wt / sum(wt)

    weighted_mean <- sum(v * wt)
    weighted_var <- sum(wt * (v - weighted_mean)^2)
    weighted_var
  }

  list(
    sensitivity = compute_weighted_var(per_cohort_metrics$sensitivity, weights),
    specificity = compute_weighted_var(per_cohort_metrics$specificity, weights)
  )
}

#' Calibrate a Panel with Neyman-Pearson Threshold
#'
#' Performs NP threshold calibration on a fitted panel using held-out data,
#' computing per-cohort metrics and weighted variance for transferability
#' assessment. Returns a `TransferablePanelResult` extending the input panel.
#'
#' @param panel A `BiomarkerPanelResult` from [fit_panel()].
#' @param x_heldout Held-out feature matrix (raw base features).
#' @param y_heldout Held-out response vector (factor with levels `"No"`, `"Yes"`).
#' @param cohort_heldout Optional cohort factor for held-out samples.
#' @param np_alpha Type I error rate for NP threshold (default 0.15).
#' @param np_delta Tolerance parameter (default 0.05).
#' @return A `TransferablePanelResult` with NP threshold, per-cohort metrics,
#'   and weighted variance.
#' @export
#' @seealso [fit_panel()], [optimize_panel_transferable()], [evaluate_panel()]
calibrate_panel <- function(panel, x_heldout, y_heldout,
                            cohort_heldout = NULL,
                            np_alpha = 0.15, np_delta = 0.05) {
  if (!inherits(panel, "BiomarkerPanelResult")) {
    stop("`panel` must be a BiomarkerPanelResult from fit_panel().", call. = FALSE)
  }
  .validate_probability(np_alpha, "np_alpha", bounds = "open")
  .validate_probability(np_delta, "np_delta", bounds = "open")

  # Get base features and feature transform from panel
  base_features <- panel@base_features
  if (is.null(base_features) || length(base_features) == 0L) {
    base_features <- panel@features
  }
  feature_transform <- panel@control$feature_transform
  if (is.null(feature_transform)) feature_transform <- "none"

  # Prepare held-out matrix: validate base features and apply the transform.
  x_ho <- as.matrix(x_heldout)
  x_selected <- .prepare_scoring_matrix(x_ho, base_features, feature_transform,
                                        context = "held-out data")$x_selected

  # Ensure response is properly formatted
  y_ho <- ensure_binary_response(y_heldout)

  # Set up cohort factor
  if (is.null(cohort_heldout)) {
    cohort_ho <- factor(rep("cohort_01", nrow(x_ho)), levels = "cohort_01")
  } else {
    cohort_ho <- factor(cohort_heldout)
  }

  # Generate predictions using stored model
  heldout_scores <- .predict_panel_model(panel@model, x_selected, cohort = cohort_ho)

  # NP threshold selection
  np_result <- .select_np_threshold(
    scores = heldout_scores,
    truth = y_ho,
    alpha = np_alpha,
    delta = np_delta
  )

  # Per-cohort metrics
  per_cohort_df <- .compute_per_cohort_metrics(
    scores = heldout_scores,
    truth = y_ho,
    cohort = cohort_ho,
    threshold = np_result$threshold
  )

  # Weighted variance
  weighted_var <- .compute_weighted_variance(per_cohort_df)

  # Construct TransferablePanelResult inheriting from panel
  new(
    "TransferablePanelResult",
    base_features = panel@base_features,
    features = panel@features,
    metrics = panel@metrics,
    objectives = panel@objectives,
    control = panel@control,
    training_data = panel@training_data,
    model = panel@model,
    np_threshold = np_result$threshold,
    np_alpha = np_alpha,
    np_delta = np_delta,
    per_cohort_metrics = per_cohort_df,
    weighted_variance = weighted_var,
    validation_metrics = list(),
    partition_info = if (!is.null(panel@control$partition_info)) {
      panel@control$partition_info
    } else {
      list()
    }
  )
}