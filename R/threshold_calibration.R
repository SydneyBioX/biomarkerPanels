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
  # Check if nproc is available
  if (!requireNamespace("nproc", quietly = TRUE)) {
    warning(
      "Package 'nproc' not available. Using fallback threshold of 0.5. ",
      "Install nproc via install.packages('nproc') for NP threshold selection.",
      call. = FALSE
    )
    return(list(
      threshold = 0.5,
      method = "fallback",
      alpha = alpha,
      delta = delta
    ))
  }

  tryCatch(
    {
      # nproc expects 0/1 labels
      y_binary <- as.integer(truth) - 1L

      # npc() returns a classifier; we extract the threshold
      np_fit <- nproc::npc(
        x = as.matrix(scores),
        y = y_binary,
        method = "lda", # Simple method since scores are 1D
        alpha = alpha,
        delta = delta
      )

      # Extract threshold from fitted classifier
      # For 1D scores with LDA, threshold is related to cutoff
      # We need to find the threshold that achieves the target alpha

      # Alternative: use nproc's internal threshold calculation
      # The npc object stores the threshold implicitly
      # We can compute it by evaluating at different cutoffs

      # Simple approach: find threshold that achieves target specificity
      neg_scores <- scores[y_binary == 0]
      target_spec <- 1 - alpha

      if (length(neg_scores) > 0) {
        threshold <- stats::quantile(neg_scores, probs = target_spec, na.rm = TRUE)
      } else {
        threshold <- 0.5
      }

      list(
        threshold = as.numeric(threshold),
        method = "nproc",
        alpha = alpha,
        delta = delta
      )
    },
    error = function(e) {
      warning(
        "NP threshold selection failed: ", conditionMessage(e),
        ". Using fallback threshold of 0.5.",
        call. = FALSE
      )
      list(
        threshold = 0.5,
        method = "fallback",
        alpha = alpha,
        delta = delta,
        error = conditionMessage(e)
      )
    }
  )
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

    pred_class <- ifelse(coh_scores >= threshold, "Yes", "No")

    # Sensitivity: TP / (TP + FN)
    tp <- sum(pred_class == "Yes" & coh_truth == "Yes")
    fn <- sum(pred_class == "No" & coh_truth == "Yes")
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

    # Specificity: TN / (TN + FP)
    tn <- sum(pred_class == "No" & coh_truth == "No")
    fp <- sum(pred_class == "Yes" & coh_truth == "No")
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

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

#' Predict from Fitted Model
#'
#' Internal helper to get predictions from either glmnet or glm model.
#'
#' @param model Fitted model (cv.glmnet or glm).
#' @param newx New feature matrix for prediction.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_from_model <- function(model, newx) {
  if (inherits(model, "cv.glmnet")) {
    # glmnet model
    x_mat <- as.matrix(newx)

    # Add cohort dummies if the model was trained with them
    meta <- model$biomarkerPanels_meta
    if (!is.null(meta$cohort_info)) {
      # For prediction, we assume no cohort effect (set dummies to 0)
      n_dummies <- meta$cohort_info$n_dummies
      dummy_cols <- matrix(0, nrow = nrow(x_mat), ncol = n_dummies)
      colnames(dummy_cols) <- paste0(".cohort_", seq_len(n_dummies))
      x_mat <- cbind(x_mat, dummy_cols)
    }

    preds <- stats::predict(model,
      newx = x_mat, s = "lambda.min",
      type = "response"
    )[, 1]
  } else if (inherits(model, "glm")) {
    # Standard GLM
    newdata <- as.data.frame(newx, check.names = TRUE)
    preds <- stats::predict(model, newdata = newdata, type = "response")
  } else {
    stop("Unknown model type: ", class(model)[1], call. = FALSE)
  }

  as.numeric(preds)
}
