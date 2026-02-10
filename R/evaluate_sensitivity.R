#' Sensitivity-Based Panel Evaluation
#'
#' Functions for evaluating biomarker panels at target sensitivity operating
#' points, useful for rule-out diagnostic applications.
#'
#' @name evaluate_sensitivity
NULL

#' Find Threshold for Target Sensitivity
#'
#' Computes the classification threshold that achieves a target sensitivity
#' level on the ROC curve. This is essential when using ROC-based optimization
#' objectives like [metric_specificity_at_sensitivity()], where the optimization
#' finds panels that perform well at a specific sensitivity level, but the
#' default evaluation cutoff (0.5) may correspond to a different operating point.
#'
#' @param panel A `BiomarkerPanelResult`.
#' @param x Validation feature matrix, or a list of matrices for multi-cohort
#'   evaluation.
#' @param y Validation binary response factor, or a list of factors when `x` is
#'   a list.
#' @param target_sensitivity The desired sensitivity level (default 0.90).
#'   Must be between 0 and 1.
#' @param positive Label treated as the positive class. Defaults to the value
#'   stored in the panel, or `"Yes"`.
#' @param feature_transform Optional feature transform name. Defaults to the
#'   transform stored in the fitted panel.
#' @param assay For `SummarizedExperiment` inputs, assay name or index.
#' @return A list with:
#'   \describe{
#'     \item{threshold}{The cutoff probability that achieves the target sensitivity}
#'     \item{sensitivity}{The actual sensitivity at this threshold (may differ
#'       slightly from target due to discrete ROC points)}
#'     \item{specificity}{The specificity achieved at this threshold}
#'     \item{target_sensitivity}{The requested target sensitivity}
#'     \item{auc}{The full AUC of the ROC curve}
#'   }
#' @export
#' @seealso [evaluate_panel()] to evaluate with a specific cutoff,
#'   [metric_specificity_at_sensitivity()] for the optimization objective.
#' @examples
#' # After optimizing with specificity_at_sensitivity objective:
#' # result <- optimize_panel(x, y, objectives = build_objectives(
#' #   metrics = c("specificity_at_sensitivity", "num_features")
#' # ))
#' #
#' # Find the threshold for 90% sensitivity:
#' # thresh_info <- find_threshold_for_sensitivity(result, x_val, y_val)
#' #
#' # Evaluate at that threshold:
#' # eval <- evaluate_panel(result, x_val, y_val, cutoff_prob = thresh_info$threshold)
find_threshold_for_sensitivity <- function(panel,
                                            x,
                                            y,
                                            target_sensitivity = 0.90,
                                            positive = NULL,
                                            feature_transform = NULL,
                                            assay = NULL) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required. Install it via install.packages('pROC').",
         call. = FALSE)
  }

  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  if (target_sensitivity < 0 || target_sensitivity > 1) {
    stop("`target_sensitivity` must be between 0 and 1.", call. = FALSE)
  }

  # Use stored positive class from training if not explicitly specified
  if (is.null(positive)) {
    stored_positive <- panel@control$positive_class
    positive <- if (!is.null(stored_positive)) stored_positive else "Yes"
  }

  # Get predictions using evaluate_panel (reuse existing infrastructure)
  eval_result <- evaluate_panel(
    panel = panel,
    x = x,
    y = y,
    cohort = NULL,
    assay = assay,
    feature_transform = feature_transform,
    cutoff_prob = 0.5,  # Doesn't matter for ROC computation
    positive = positive
  )

  scores <- eval_result$scores
  truth <- eval_result$roc$evalm$obs
  truth <- ensure_binary_response(truth, positive = positive)

  # Build ROC curve
  roc_obj <- pROC::roc(
    response = truth,
    predictor = scores,
    levels = c(setdiff(levels(truth), positive), positive),
    direction = "<",
    quiet = TRUE
  )

  # Extract ROC curve data to find actual threshold

  # pROC::coords can return NA for interpolated thresholds, so we find the
  # actual threshold from the ROC curve that achieves at least target sensitivity
  roc_data <- data.frame(
    threshold = roc_obj$thresholds,
    sensitivity = roc_obj$sensitivities,
    specificity = roc_obj$specificities
  )

  # Find all points where sensitivity >= target
  valid_idx <- which(roc_data$sensitivity >= target_sensitivity)

  if (length(valid_idx) == 0L) {
    # If no point achieves target sensitivity, use the maximum sensitivity point
    best_idx <- which.max(roc_data$sensitivity)
    warning(
      sprintf(
        "Target sensitivity %.2f cannot be achieved. Maximum achievable is %.2f.",
        target_sensitivity, roc_data$sensitivity[best_idx]
      ),
      call. = FALSE
    )
  } else {
    # Among points that achieve target sensitivity, pick the one with highest specificity
    # This gives the most conservative threshold that still meets the sensitivity target
    best_idx <- valid_idx[which.max(roc_data$specificity[valid_idx])]
  }

  threshold <- roc_data$threshold[best_idx]
  # Handle edge case where threshold is -Inf (predict all positive)
  if (is.infinite(threshold) && threshold < 0) {
    # Use the minimum actual score as threshold
    threshold <- min(scores) - .Machine$double.eps
  }

  list(
    threshold = threshold,
    sensitivity = roc_data$sensitivity[best_idx],
    specificity = roc_data$specificity[best_idx],
    target_sensitivity = target_sensitivity,
    auc = as.numeric(pROC::auc(roc_obj))
  )
}

#' Evaluate Panel at Target Sensitivity
#'
#' Convenience wrapper that combines [find_threshold_for_sensitivity()] and
#' [evaluate_panel()] to evaluate a panel at the operating point corresponding
#' to a target sensitivity level. This is the recommended way to evaluate
#' panels optimized with [metric_specificity_at_sensitivity()].
#'
#' @inheritParams find_threshold_for_sensitivity
#' @param objectives Optional override for objectives. Defaults to sensitivity,
#'   specificity, and auc.
#' @return Same as [evaluate_panel()], but with `cutoff_prob` set to the
#'   threshold that achieves the target sensitivity. The result also includes
#'   a `target_sensitivity_info` element with details about the threshold
#'   selection.
#' @export
#' @seealso [find_threshold_for_sensitivity()], [evaluate_panel()].
#' @examples
#' # After optimizing with specificity_at_sensitivity:
#' # result <- optimize_panel(x, y, objectives = build_objectives(
#' #   metrics = c("specificity_at_sensitivity", "num_features")
#' # ))
#' #
#' # Evaluate at 90% sensitivity operating point:
#' # eval <- evaluate_panel_at_sensitivity(result, x_val, y_val)
evaluate_panel_at_sensitivity <- function(panel,
                                           x,
                                           y,
                                           target_sensitivity = 0.90,
                                           objectives = define_objectives(
                                             metrics = c("sensitivity", "specificity", "auc")
                                           ),
                                           positive = NULL,
                                           feature_transform = NULL,
                                           assay = NULL) {
  # Find the optimal threshold
  thresh_info <- find_threshold_for_sensitivity(
    panel = panel,
    x = x,
    y = y,
    target_sensitivity = target_sensitivity,
    positive = positive,
    feature_transform = feature_transform,
    assay = assay
  )

  # Evaluate at that threshold
  eval_result <- evaluate_panel(
    panel = panel,
    x = x,
    y = y,
    objectives = objectives,
    cohort = NULL,
    assay = assay,
    feature_transform = feature_transform,
    cutoff_prob = thresh_info$threshold,
    positive = positive
  )

  # Add threshold info to result
  eval_result$target_sensitivity_info <- thresh_info

  eval_result
}
