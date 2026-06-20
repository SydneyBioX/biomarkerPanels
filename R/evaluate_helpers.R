#' Panel Evaluation Helpers
#'
#' Internal helpers used by [evaluate_panel()] to build confusion matrices,
#' compute ROC curves (via the C++ backend), and obtain predictions from a
#' stored glmnet model.
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

#' Predict Using Stored glmnet Model
#'
#' Helper function to make predictions from a cv.glmnet model stored in a
#' BiomarkerPanelResult. Handles cohort dummy variables using the metadata
#' stored during training.
#'
#' @param model A cv.glmnet model object with biomarkerPanels_meta attribute.
#' @param x_selected Matrix of selected features for prediction.
#' @param cohort_vec Factor of cohort membership for samples.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_glmnet_model <- function(model, x_selected, cohort_vec) {
  x_mat <- as.matrix(x_selected)

  # Get metadata stored during training

  meta <- model$biomarkerPanels_meta

  # Add cohort dummies if the model was trained with them.
  # Always zero out — predictions should be cohort-agnostic for

  # transferability. Cohort-aware metrics split by cohort downstream.
  if (!is.null(meta$cohort_info)) {
    n_dummies <- meta$cohort_info$n_dummies
    dummy_cols <- matrix(0, nrow = nrow(x_mat), ncol = n_dummies)
    colnames(dummy_cols) <- paste0(".cohort_", seq_len(n_dummies))
    x_mat <- cbind(x_mat, dummy_cols)
  }

  # Predict using lambda.min
  preds <- stats::predict(model, newx = x_mat, s = "lambda.min", type = "response")[, 1]

  preds
}
