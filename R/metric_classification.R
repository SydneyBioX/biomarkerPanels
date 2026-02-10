#' Classification Metric Functions
#'
#' Cutoff-dependent metric functions for binary classification that require a
#' probability threshold. These compute metrics like sensitivity, specificity,
#' precision, NPV, F1, and balanced accuracy.
#'
#' @name metric_classification
NULL

#' Sensitivity (True Positive Rate)
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param cutoff_prob Classification probability cutoff applied to `scores`.
#'   Ignored when `cutoff_strategy` is not `"fixed"`.
#' @param cutoff_strategy Strategy for computing cutoff. One of `"fixed"`
#'   (use `cutoff_prob`), `"prevalence"` (cutoff = class prevalence), or
#'   `"youden"` (optimal Youden's J). Default is `"fixed"`.
#' @param positive Label treated as the positive ("event") class.
#' @return Sensitivity between 0 and 1, or `NA_real_` if undefined.
#' @export
metric_sensitivity <- function(truth, scores = NULL, selected = NULL,
                             cutoff_prob = 0.5,
                             cutoff_strategy = c("fixed", "prevalence", "youden"),
                             positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute sensitivity.", call. = FALSE)
  }
  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy, positive)
  positives <- truth == positive
  predicted <- scores >= cutoff
  tp <- sum(positives & predicted)
  fn <- sum(positives & !predicted)
  if ((tp + fn) == 0) {
    return(NA_real_)
  }
  tp / (tp + fn)
}

#' Specificity (True Negative Rate)
#'
#' @inheritParams metric_sensitivity
#' @return Specificity between 0 and 1, or `NA_real_` if undefined.
#' @export
metric_specificity <- function(truth, scores = NULL, selected = NULL,
                             cutoff_prob = 0.5,
                             cutoff_strategy = c("fixed", "prevalence", "youden"),
                             positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute specificity.", call. = FALSE)
  }
  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy, positive)
  negatives <- truth != positive
  predicted <- scores >= cutoff
  tn <- sum(negatives & !predicted)
  fp <- sum(negatives & predicted)
  if ((tn + fp) == 0) {
    return(NA_real_)
  }
  tn / (tn + fp)
}

#' Precision (Positive Predictive Value)
#'
#' Computes the proportion of predicted positives that are true positives.
#' Also known as Positive Predictive Value (PPV).
#'
#' @inheritParams metric_sensitivity
#' @return Precision between 0 and 1, or `NA_real_` if undefined.
#' @export
metric_precision <- function(truth, scores = NULL, selected = NULL,
                           cutoff_prob = 0.5,
                           cutoff_strategy = c("fixed", "prevalence", "youden"),
                           positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute precision.", call. = FALSE)
  }
  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy, positive)

  positives <- truth == positive
  predicted <- scores >= cutoff

  tp <- sum(positives & predicted)
  fp <- sum(!positives & predicted)

  if ((tp + fp) == 0) {
    return(NA_real_)
  }
  tp / (tp + fp)
}

#' Negative Predictive Value (NPV)
#'
#' Computes the proportion of predicted negatives that are true negatives.
#' Important for rule-out diagnostics where a negative prediction should
#' reliably exclude disease.
#'
#' @inheritParams metric_sensitivity
#' @return NPV between 0 and 1, or `NA_real_` if undefined.
#' @export
metric_npv <- function(truth, scores = NULL, selected = NULL,
                     cutoff_prob = 0.5,
                     cutoff_strategy = c("fixed", "prevalence", "youden"),
                     positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute NPV.", call. = FALSE)
  }
  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy, positive)

  negatives <- truth != positive
  predicted_negative <- scores < cutoff

  tn <- sum(negatives & predicted_negative)
  fn <- sum(!negatives & predicted_negative)

  if ((tn + fn) == 0) {
    return(NA_real_)
  }
  tn / (tn + fn)
}

#' F1 Score
#'
#' Computes the harmonic mean of precision and recall (sensitivity).
#' F1 = 2 * (precision * recall) / (precision + recall).
#'
#' @inheritParams metric_sensitivity
#' @return F1 score between 0 and 1, or `NA_real_` if undefined.
#' @export
metric_f1 <- function(truth, scores = NULL, selected = NULL,
                    cutoff_prob = 0.5,
                    cutoff_strategy = c("fixed", "prevalence", "youden"),
                    positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute F1 score.", call. = FALSE)
  }
  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy, positive)

  positives <- truth == positive
  predicted <- scores >= cutoff

  tp <- sum(positives & predicted)
  fp <- sum(!positives & predicted)
  fn <- sum(positives & !predicted)

  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)

  if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    return(NA_real_)
  }
  2 * (precision * recall) / (precision + recall)
}

#' Balanced Accuracy
#'
#' @inheritParams metric_sensitivity
#' @return Mean of sensitivity and specificity.
#' @export
metric_balanced_accuracy <- function(truth, scores = NULL, selected = NULL,
                                   cutoff_prob = 0.5,
                                   cutoff_strategy = c("fixed", "prevalence", "youden"),
                                   positive = "Yes") {
  cutoff_strategy <- match.arg(cutoff_strategy)
  sens <- metric_sensitivity(truth, scores, selected,
                           cutoff_prob = cutoff_prob,
                           cutoff_strategy = cutoff_strategy,
                           positive = positive)
  spec <- metric_specificity(truth, scores, selected,
                           cutoff_prob = cutoff_prob,
                           cutoff_strategy = cutoff_strategy,
                           positive = positive)
  if (is.na(sens) || is.na(spec)) {
    return(NA_real_)
  }
  (sens + spec) / 2
}

#' Compute Classification Cutoff Based on Strategy
#'
#' Internal helper to compute the optimal classification cutoff based on
#' the specified strategy.
#'
#' @param truth Binary response factor.
#' @param scores Numeric prediction scores.
#' @param cutoff_prob Fixed cutoff (used when strategy = "fixed").
#' @param cutoff_strategy Strategy for computing cutoff.
#' @param positive Label for positive class.
#' @return Numeric cutoff value.
#' @keywords internal
.compute_cutoff <- function(truth, scores, cutoff_prob, cutoff_strategy, positive) {
  if (cutoff_strategy == "fixed") {
    return(cutoff_prob)
  }

  if (cutoff_strategy == "prevalence") {
    # Cutoff equals training class prevalence
    return(mean(truth == positive))
  }

  if (cutoff_strategy == "youden") {
    # Optimal Youden's J point from ROC curve
    # Requires pROC package
    if (!requireNamespace("pROC", quietly = TRUE)) {
      stop(
        "The 'pROC' package is required for Youden cutoff strategy. ",
        "Install it via install.packages('pROC') or use cutoff_strategy = 'fixed'.",
        call. = FALSE
      )
    }

    tryCatch({
      roc_obj <- pROC::roc(
        response = truth,
        predictor = scores,
        levels = c(setdiff(levels(truth), positive), positive),
        direction = "<",
        quiet = TRUE
      )

      # Find optimal threshold using Youden's J
      coords <- pROC::coords(roc_obj, "best", ret = "threshold",
                             best.method = "youden")
      if (length(coords) == 0 || is.na(coords[1])) {
        stop(
          "Youden cutoff computation returned no valid threshold. ",
          "This may indicate degenerate ROC data. ",
          "Consider using cutoff_strategy = 'fixed' or 'prevalence' instead.",
          call. = FALSE
        )
      }
      as.numeric(coords[1])
    }, error = function(e) {
      stop(
        "Youden cutoff computation failed: ", conditionMessage(e),
        "\nThis typically indicates insufficient class separation in the data. ",
        "Consider using cutoff_strategy = 'fixed' or 'prevalence' instead.",
        call. = FALSE
      )
    })
  }
}
