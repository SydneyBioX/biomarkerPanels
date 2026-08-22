#' Discrimination Metric Functions
#'
#' Threshold-free metric functions that measure classifier discrimination
#' ability without requiring a probability cutoff. Includes AUC, partial AUC,
#' specificity at fixed sensitivity, sensitivity at fixed specificity, and
#' panel size metrics.
#'
#' @name metric_discrimination
NULL

#' Area Under the ROC Curve (AUC)
#'
#' Computes the Area Under the ROC Curve using [pROC::auc()]. AUC represents
#' the probability that a random positive instance receives a higher score
#' than a random negative instance.
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param positive Label treated as the positive ("event") class.
#' @return AUC between 0 and 1, or `NA_real_` if computation fails.
#' @export
metric_auc <- function(truth, scores = NULL, selected = NULL,
                     positive = "Yes") {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required for AUC computation. ",
         "Install it via install.packages('pROC').", call. = FALSE)
  }
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute AUC.", call. = FALSE)
  }

  pos_count <- sum(truth == positive)
  neg_count <- sum(truth != positive)
  if (pos_count == 0L || neg_count == 0L) {
    stop(
      "Cannot compute AUC: ",
      if (pos_count == 0L) "no positive samples" else "no negative samples",
      " in the data. Both classes must be represented.",
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
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) {
    stop(
      "AUC computation failed: ", conditionMessage(e),
      "\nThis typically indicates insufficient data or class imbalance.",
      call. = FALSE
    )
  })
}

#' Partial AUC in a Sensitivity-Focused Region
#'
#' Computes the partial area under the ROC curve restricted to sensitivity
#' values at or above a specified floor. Useful for "rule-out" diagnostic
#' contexts where high sensitivity is mandatory. Uses [pROC::auc()] with
#' `partial.auc.focus = "sens"`.
#'
#' @inheritParams metric_auc
#' @param sens_floor Minimum sensitivity threshold defining the pAUC region
#'   (default 0.90). The partial AUC is computed only where sensitivity >= this
#'   value.
#' @param partial_auc_correct Logical; apply McClish correction to normalize
#'   partial AUC to 0-1 scale (default TRUE).
#' @return Partial AUC value, or `NA_real_` if computation fails.
#' @export
metric_pauc <- function(truth, scores = NULL, selected = NULL,
                      positive = "Yes", sens_floor = 0.90,
                      partial_auc_correct = TRUE) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required for partial AUC. ",
         "Install it via install.packages('pROC').", call. = FALSE)
  }
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute partial AUC.", call. = FALSE)
  }

  if (sens_floor < 0 || sens_floor > 1) {
    stop("`sens_floor` must be between 0 and 1.", call. = FALSE)
  }

  tryCatch({
    # pROC expects: response = truth, predictor = scores
    # partial.auc = c(sens_floor, 1) means sensitivity range [sens_floor, 1]
    # partial.auc.focus = "sens" focuses on sensitivity axis
    roc_obj <- pROC::roc(
      response = truth,
      predictor = scores,
      levels = c(setdiff(levels(truth), positive), positive),
      direction = "<",
      quiet = TRUE
    )

    pauc <- pROC::auc(
      roc_obj,
      partial.auc = c(sens_floor, 1),
      partial.auc.focus = "sens",
      partial.auc.correct = partial_auc_correct
    )

    as.numeric(pauc)
  }, error = function(e) {
    stop(
      "Partial AUC computation failed: ", conditionMessage(e),
      "\nThis typically indicates insufficient data or class imbalance. ",
      "Ensure both positive and negative classes have adequate representation.",
      call. = FALSE
    )
  })
}

#' Specificity at a Fixed Sensitivity Threshold
#'
#' Computes the specificity achievable when the classifier operates at a
#' specified sensitivity level. This is useful for rule-out diagnostic
#' applications where a minimum sensitivity is mandatory, and we want to
#' know the best achievable specificity at that operating point.
#'
#' Uses [pROC::coords()] to find the operating point on the ROC curve
#' corresponding to the target sensitivity, then returns the specificity
#' at that point. If the exact target sensitivity cannot be achieved,
#' pROC interpolates to the nearest achievable point.
#'
#' Degenerate inputs (a single-class subset, or a `coords()` call that yields
#' nothing usable) return `0` — the worst achievable value — rather than
#' erroring or returning `NA_real_`. This keeps the metric safe to use
#' directly as an NSGA objective: an error would abort the whole run, and an
#' `NA` is mapped to the `1e6` infeasibility penalty, which conflates
#' degeneracy with constraint violation, distorts NSGA-III normalization, and
#' still surfaces as `NA` in post-hoc solution re-evaluation. Argument errors
#' (missing `scores`, out-of-range `target_sensitivity`) still signal an
#' error.
#'
#' @inheritParams metric_auc
#' @param target_sensitivity The sensitivity threshold at which to evaluate
#'   specificity (default 0.90). Must be between 0 and 1.
#' @return Specificity at the target sensitivity, between 0 and 1. Returns `0`
#'   for degenerate inputs.
#' @seealso [metric_pauc()] for partial AUC in the high-sensitivity region,
#'   [metric_sensitivity()] and [metric_specificity()] for threshold-based metrics.
#' @export
#' @examples
#' truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
#' scores <- c(runif(50, 0, 0.6), runif(50, 0.4, 1))
#'
#' # Specificity when sensitivity is fixed at 90%
#' metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.90)
#'
#' # Specificity when sensitivity is fixed at 95%
#' metric_specificity_at_sensitivity(truth, scores, target_sensitivity = 0.95)
metric_specificity_at_sensitivity <- function(truth, scores = NULL, selected = NULL,
                                             positive = "Yes",
                                             target_sensitivity = 0.90) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required for this metric. ",
         "Install it via install.packages('pROC').", call. = FALSE)
  }
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute specificity at sensitivity.", call. = FALSE)
  }

  if (target_sensitivity < 0 || target_sensitivity > 1) {
    stop("`target_sensitivity` must be between 0 and 1.", call. = FALSE)
  }

  # Degenerate cases return the worst achievable value rather than erroring
  # (which aborts the run) or NA (which is penalised as 1e6 infeasibility).
  if (sum(truth == positive) == 0L || sum(truth != positive) == 0L) {
    return(0)
  }

  tryCatch({
    roc_obj <- pROC::roc(
      response = truth,
      predictor = scores,
      levels = c(setdiff(levels(truth), positive), positive),
      direction = "<",
      quiet = TRUE
    )

    spec_value <- pROC::coords(
      roc_obj,
      x = target_sensitivity,
      input = "sensitivity",
      ret = "specificity"
    )$specificity

    # coords may return several ties; take the best achievable specificity
    if (length(spec_value) == 0L) {
      return(0)
    }
    spec_value <- suppressWarnings(max(spec_value, na.rm = TRUE))
    if (!is.finite(spec_value)) 0 else as.numeric(spec_value)
  }, error = function(e) 0)
}

#' Sensitivity at a Fixed Specificity Threshold
#'
#' Computes the sensitivity achievable when the classifier operates at a
#' specified specificity level. This is the transpose of
#' [metric_specificity_at_sensitivity()], and is the natural objective for
#' "rule-in" applications where a minimum specificity (i.e. a tolerable false
#' positive rate) is mandatory and we want the best achievable detection rate
#' at that operating point.
#'
#' Uses [pROC::coords()] to find the operating point on the ROC curve
#' corresponding to the target specificity, then returns the sensitivity at
#' that point. If the exact target specificity cannot be achieved, pROC
#' interpolates to the nearest achievable point.
#'
#' Like [metric_specificity_at_sensitivity()], degenerate inputs (a
#' single-class subset, or a `coords()` call that yields nothing usable)
#' return `0` — the worst achievable value — rather than erroring or returning
#' `NA_real_`. This keeps the metric safe to use directly as an NSGA
#' objective: an error would abort the whole run, and an `NA` is mapped to the
#' `1e6` infeasibility penalty, which conflates degeneracy with constraint
#' violation, distorts NSGA-III normalization, and still surfaces as `NA` in
#' post-hoc solution re-evaluation. Argument errors (missing `scores`,
#' out-of-range `target_specificity`) still signal an error.
#'
#' @inheritParams metric_auc
#' @param target_specificity The specificity threshold at which to evaluate
#'   sensitivity (default 0.90). Must be between 0 and 1.
#' @return Sensitivity at the target specificity, between 0 and 1. Returns `0`
#'   for degenerate inputs.
#' @seealso [metric_specificity_at_sensitivity()] for the transpose,
#'   [metric_pauc()] for partial AUC in the high-sensitivity region,
#'   [metric_sensitivity()] and [metric_specificity()] for threshold-based metrics.
#' @export
#' @examples
#' truth <- factor(c(rep("No", 50), rep("Yes", 50)), levels = c("No", "Yes"))
#' scores <- c(runif(50, 0, 0.6), runif(50, 0.4, 1))
#'
#' # Sensitivity when specificity is fixed at 90%
#' metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.90)
#'
#' # Sensitivity when specificity is fixed at 95%
#' metric_sensitivity_at_specificity(truth, scores, target_specificity = 0.95)
metric_sensitivity_at_specificity <- function(truth, scores = NULL, selected = NULL,
                                              positive = "Yes",
                                              target_specificity = 0.90) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required for this metric. ",
         "Install it via install.packages('pROC').", call. = FALSE)
  }
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute sensitivity at specificity.", call. = FALSE)
  }

  if (target_specificity < 0 || target_specificity > 1) {
    stop("`target_specificity` must be between 0 and 1.", call. = FALSE)
  }

  # Degenerate cases return the worst achievable value rather than erroring
  # (which aborts the run) or NA (which is penalised as 1e6 infeasibility).
  if (sum(truth == positive) == 0L || sum(truth != positive) == 0L) {
    return(0)
  }

  tryCatch({
    roc_obj <- pROC::roc(
      response = truth,
      predictor = scores,
      levels = c(setdiff(levels(truth), positive), positive),
      direction = "<",
      quiet = TRUE
    )

    sens_value <- pROC::coords(
      roc_obj,
      x = target_specificity,
      input = "specificity",
      ret = "sensitivity"
    )$sensitivity

    # coords may return several ties; take the best achievable sensitivity
    if (length(sens_value) == 0L) {
      return(0)
    }
    sens_value <- suppressWarnings(max(sens_value, na.rm = TRUE))
    if (!is.finite(sens_value)) 0 else as.numeric(sens_value)
  }, error = function(e) 0)
}

#' Panel Size Penalty
#'
#' Counts the number of base features (individual genes) in a panel, not
#' transformed features such as pairwise ratios.
#'
#' @param truth Ignored; kept for signature compatibility.
#' @param scores Ignored; kept for signature compatibility.
#' @param selected Vector of selected base features (character, numeric, or logical).
#' @param ... Additional arguments (ignored).
#' @return Count of selected base biomarkers (genes).
#' @export
metric_num_features <- function(truth = NULL, scores = NULL, selected = NULL, ...) {
  if (is.null(selected)) {
    return(0)
  }
  if (is.logical(selected)) {
    return(sum(selected))
  }
  length(selected)
}
