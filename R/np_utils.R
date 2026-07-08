#' Neyman-Pearson Utility Functions
#'
#' Internal utility functions for Neyman-Pearson (NP) penalized regression
#' support via the nproc package.
#'
#' @name np_utils
#' @keywords internal
NULL

#' Check nproc Package Availability
#'
#' Validates that the nproc package is installed and available.
#'
#' @return Invisible TRUE if available.
#' @keywords internal
.check_nproc_available <- function() {
  if (!requireNamespace("nproc", quietly = TRUE)) {
    stop(
      "Package 'nproc' is required for Neyman-Pearson fitting but not installed.\n",
      "Install it via: install.packages('nproc')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Prepare Labels for NP Classification
#'
#' Core label preparation logic that handles the label flipping required
#' for rule-out tests. When `minimize_FPR = FALSE`, labels are flipped so
#' that nproc controls the false negative rate instead of the false positive rate.
#'
#' @param truth Binary response factor with levels c("No", "Yes").
#' @param minimize_FPR Logical; if TRUE (default), controls FPR (standard).
#'   If FALSE, flips labels to control FNR (rule-out test).
#' @return List with:
#'   - `y_nproc`: Integer vector (0/1) for nproc
#'   - `labels_flipped`: Logical indicating if labels were flipped
#'   - `original_positive`: The original positive class label
#' @keywords internal
.prepare_np_labels <- function(truth, minimize_FPR = TRUE) {
  truth <- ensure_binary_response(truth)

  if (minimize_FPR) {
    # Standard encoding: No=0, Yes=1 - nproc controls Type I error (FPR)
    y_nproc <- as.integer(truth) - 1L
    labels_flipped <- FALSE
  } else {
    # Flipped encoding: Yes=0, No=1 - nproc controls FNR (miss rate)
    # This achieves high sensitivity (low FNR) at the expense of specificity
    y_nproc <- 1L - (as.integer(truth) - 1L)
    labels_flipped <- TRUE
  }

  list(
    y_nproc = y_nproc,
    labels_flipped = labels_flipped,
    original_positive = "Yes"
  )
}

#' Invert NP Predictions
#'
#' When labels were flipped for rule-out classification, predictions must be
#' inverted back to the original scale.
#'
#' @param predictions Numeric vector of predicted probabilities from nproc.
#' @param labels_flipped Logical indicating if labels were flipped.
#' @return Numeric vector of predictions on the original scale.
#' @keywords internal
.invert_np_predictions <- function(predictions, labels_flipped) {
  if (labels_flipped) {
    1 - predictions
  } else {
    predictions
  }
}

#' Validate NP Parameters
#'
#' Validates parameters for Neyman-Pearson classification.
#'
#' @param alpha Type I error rate (must be in (0, 1)).
#' @param delta Tolerance parameter (must be in (0, 1)).
#' @param method nproc method name.
#' @return Invisible TRUE if all parameters valid.
#' @keywords internal
.validate_np_params <- function(alpha, delta, method) {
  .validate_probability(alpha, "alpha", bounds = "open")
  .validate_probability(delta, "delta", bounds = "open")

  # Must match nproc::npc's accepted method names exactly (nproc uses
  # "randomforest", not "rf").
  valid_methods <- c("penlog", "logistic", "svm", "lda", "slda",
                     "nb", "nnb", "ada", "tree", "randomforest")
  if (!is.character(method) || length(method) != 1L ||
      !method %in% valid_methods) {
    stop(
      "`method` must be one of: ", paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Predict from NP Classifier
#'
#' Extracts continuous probability scores from an npc model object.
#' Handles label inversion when labels were flipped for rule-out classification.
#'
#' @param model An npc object from nproc::npc().
#' @param x_selected Matrix of selected features for prediction.
#' @return Numeric vector of predicted probabilities on the original scale.
#' @keywords internal
.predict_np_model <- function(model, x_selected) {
  .check_nproc_available()

  x_mat <- as.matrix(x_selected)

  # Extract labels_flipped from model metadata

  labels_flipped <- isTRUE(model$biomarkerPanels_meta$labels_flipped)

  # Get predictions from nproc - use predict method for npc objects
  np_pred <- stats::predict(model, newx = x_mat)

  # Extract continuous scores
  # npc predict returns a list with pred.label and pred.score
  if (is.list(np_pred) && "pred.score" %in% names(np_pred)) {
    scores <- as.numeric(np_pred$pred.score)
  } else if (is.numeric(np_pred)) {
    scores <- np_pred
  } else {
    stop(
      "Unexpected prediction format from npc model. ",
      "Expected list with pred.score or numeric vector.",
      call. = FALSE
    )
  }

  # Invert if labels were flipped
  .invert_np_predictions(scores, labels_flipped)
}
