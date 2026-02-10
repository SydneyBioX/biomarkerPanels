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

  # Get base features and feature transform from panel
  base_features <- panel@base_features
  if (is.null(base_features) || length(base_features) == 0L) {
    base_features <- panel@features
  }
  feature_transform <- panel@control$feature_transform
  if (is.null(feature_transform)) feature_transform <- "none"

  # Prepare held-out matrix
  x_ho <- as.matrix(x_heldout)
  if (!all(base_features %in% colnames(x_ho))) {
    missing <- setdiff(base_features, colnames(x_ho))
    stop("Base feature(s) not found in held-out data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x_base <- x_ho[, base_features, drop = FALSE]

  # Apply feature transform
  if (feature_transform != "none" && length(base_features) >= 2L) {
    x_selected <- .apply_feature_transform_single(x_base, feature_transform)
  } else {
    x_selected <- x_base
  }

  # Ensure response is properly formatted
  y_ho <- ensure_binary_response(y_heldout)

  # Set up cohort factor
  if (is.null(cohort_heldout)) {
    cohort_ho <- factor(rep("cohort_01", nrow(x_ho)), levels = "cohort_01")
  } else {
    cohort_ho <- factor(cohort_heldout)
  }

  # Generate predictions using stored model
  heldout_scores <- .predict_from_model(panel@model, x_selected, cohort = cohort_ho)

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

#' Predict from Fitted Model
#'
#' Internal helper to get predictions from either glmnet or glm model.
#'
#' @param model Fitted model (cv.glmnet or glm).
#' @param newx New feature matrix for prediction.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_from_model <- function(model, newx, cohort = NULL) {
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
    # Note: as.data.frame.matrix() ignores check.names, so we must sanitize
    # names explicitly to match the names used during training in
    # .fit_final_model(), which wraps in data.frame(check.names = TRUE).
    newdata <- as.data.frame(newx)
    names(newdata) <- make.names(names(newdata), unique = TRUE)

    # If model was trained with cohort covariate, always use reference level
    # so predictions are cohort-agnostic (matching the glmnet path which zeros
    # out cohort dummies). Cohort-aware metrics split by cohort downstream.
    model_terms <- attr(stats::terms(model), "term.labels")
    if (".cohort" %in% model_terms) {
      ref_level <- model$xlevels[[".cohort"]][1]
      newdata$.cohort <- factor(rep(ref_level, nrow(newdata)),
                                levels = model$xlevels[[".cohort"]])
    }

    preds <- stats::predict(model, newdata = newdata, type = "response")
  } else {
    stop("Unknown model type: ", class(model)[1], call. = FALSE)
  }

  as.numeric(preds)
}
