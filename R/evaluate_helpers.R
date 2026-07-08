#' Panel Evaluation Helpers
#'
#' Internal helpers used by [evaluate_panel()] to build confusion matrices and
#' compute ROC curves (via the C++ backend). Stored-model prediction lives in
#' [.predict_panel_model()] (`R/model_prediction.R`).
#'
#' @name evaluate_helpers
#' @keywords internal
NULL

#' Compute Confusion Matrix
#'
#' @param truth Binary response factor.
#' @param scores Numeric prediction scores.
#' @param cutoff_prob Classification cutoff.
#' @param positive Label for positive class.
#' @return A 2x2 confusion matrix with cutoff, positive, and negative attributes.
#' @keywords internal
.compute_confusion_matrix <- function(truth, scores, cutoff_prob, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }
  negative <- setdiff(levels_truth, positive)
  if (!length(negative)) {
    stop("Unable to determine negative class for confusion matrix.", call. = FALSE)
  }
  negative <- negative[1]
  predicted_positive <- scores >= cutoff_prob
  predicted <- ifelse(predicted_positive, positive, negative)
  predicted <- factor(predicted, levels = levels_truth)
  table_res <- table(
    truth = truth,
    predicted = predicted
  )
  structure(
    as.matrix(table_res),
    cutoff = cutoff_prob,
    positive = positive,
    negative = negative
  )
}

#' Compute ROC Curve
#'
#' Computes the ROC curve for a set of predictions. Uses C++ for performance.
#'
#' @param truth A binary response vector or factor.
#' @param scores Numeric vector of predicted scores or probabilities.
#' @param positive The label for the positive class.
#' @return A data.frame with columns `threshold`, `tpr`, `fpr`, `sensitivity`,
#'   and `specificity`.
#' @keywords internal
.compute_roc_curve <- function(truth, scores, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }

  # Convert truth to logical vector for C++
  is_positive <- truth == positive

  # Call C++ implementation
  df <- .compute_roc_curve_cpp(as.numeric(scores), is_positive)

  # Sort by fpr, tpr for consistency with pure R version
  df[order(df$fpr, df$tpr), , drop = FALSE]
}

