#' Default loss functions for biomarker panel optimization.
#'
#' These functions compute scalar metrics such as sensitivity, specificity,
#' AUC, and panel size. They share a common signature so that the package can
#' compose them into multi-objective optimization routines.
#'
#' Each loss accepts a binary response (`truth`), numeric prediction scores
#' (`scores`), and optionally the currently selected feature set (`selected`).
#' Additional parameters (for example, classification cutoffs) can be supplied
#' via [`build_objectives()`].
#'
#' @name loss_functions
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
loss_sensitivity <- function(truth, scores = NULL, selected = NULL,
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
#' @inheritParams loss_sensitivity
#' @return Specificity between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_specificity <- function(truth, scores = NULL, selected = NULL,
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

#' Area Under the ROC Curve (AUC)
#'
#' Computes the Area Under the ROC Curve using [pROC::auc()]. AUC represents
#' the probability that a random positive instance receives a higher score
#' than a random negative instance.
#'
#' @inheritParams loss_sensitivity
#' @return AUC between 0 and 1, or `NA_real_` if computation fails.
#' @export
loss_auc <- function(truth, scores = NULL, selected = NULL,
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

#' F1 Score
#'
#' Computes the harmonic mean of precision and recall (sensitivity).
#' F1 = 2 * (precision * recall) / (precision + recall).
#'
#' @inheritParams loss_sensitivity
#' @return F1 score between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_f1 <- function(truth, scores = NULL, selected = NULL,
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

#' Precision (Positive Predictive Value)
#'
#' Computes the proportion of predicted positives that are true positives.
#' Also known as Positive Predictive Value (PPV).
#'
#' @inheritParams loss_sensitivity
#' @return Precision between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_precision <- function(truth, scores = NULL, selected = NULL,
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
#' @inheritParams loss_sensitivity
#' @return NPV between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_npv <- function(truth, scores = NULL, selected = NULL,
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

#' Panel Size Penalty
#'
#' @inheritParams loss_sensitivity
#' @param selected Vector of selected features (character, numeric, or logical).
#' @return Count of selected biomarkers.
#' @export
loss_num_features <- function(truth = NULL, scores = NULL, selected = NULL, ...) {
  if (is.null(selected)) {
    return(0)
  }
  if (is.logical(selected)) {
    return(sum(selected))
  }
  length(selected)
}

#' Balanced Accuracy
#'
#' @inheritParams loss_sensitivity
#' @return Mean of sensitivity and specificity.
#' @export
loss_balanced_accuracy <- function(truth, scores = NULL, selected = NULL,
                                   cutoff_prob = 0.5,
                                   cutoff_strategy = c("fixed", "prevalence", "youden"),
                                   positive = "Yes") {
  cutoff_strategy <- match.arg(cutoff_strategy)
  sens <- loss_sensitivity(truth, scores, selected,
                           cutoff_prob = cutoff_prob,
                           cutoff_strategy = cutoff_strategy,
                           positive = positive)
  spec <- loss_specificity(truth, scores, selected,
                           cutoff_prob = cutoff_prob,
                           cutoff_strategy = cutoff_strategy,
                           positive = positive)
  if (is.na(sens) || is.na(spec)) {
    return(NA_real_)
  }
  (sens + spec) / 2
}

# -----------------------------------------------------------------------------
# Factory for Cohort-Aware Loss Functions
# -----------------------------------------------------------------------------
# This factory pattern eliminates code duplication across cohort-aware losses.
# Each cohort-aware loss applies a base metric to each cohort and aggregates.

#' Create a Cohort-Aware Loss Function
#'
#' Factory function that wraps a base loss function to compute per-cohort
#' values and aggregate them. Reduces code duplication for cohort-aware metrics.
#'
#' @param base_loss The base loss function (e.g., loss_sensitivity).
#' @param aggregator Function to aggregate cohort values (e.g., min, max).
#' @param metric_name Name of the metric for error messages.
#' @param single_cohort_fallback Value to return for single cohort with gap
#'   aggregation (default NULL uses base_loss result).
#' @return A cohort-aware loss function.
#' @keywords internal
.make_cohort_aware_loss <- function(base_loss, aggregator, metric_name,
                                     single_cohort_fallback = NULL) {
  function(truth, scores = NULL, selected = NULL,
           cutoff_prob = 0.5,
           cutoff_strategy = c("fixed", "prevalence", "youden"),
           positive = "Yes",
           cohort = NULL) {
    # Standardize truth once at the start to get consistent factor levels
    truth <- ensure_binary_response(truth)
    truth_levels <- levels(truth)

    if (is.null(scores)) {
      stop("`scores` must be supplied to compute ", metric_name, ".", call. = FALSE)
    }
    cutoff_strategy <- match.arg(cutoff_strategy)

    # Fall back to base loss when no cohort provided
    if (is.null(cohort)) {
      base_result <- base_loss(truth, scores, selected,
                               cutoff_prob = cutoff_prob,
                               cutoff_strategy = cutoff_strategy,
                               positive = positive)
      if (!is.null(single_cohort_fallback)) {
        return(single_cohort_fallback * base_result)
      }
      return(base_result)
    }

    cohort <- factor(cohort)
    if (length(cohort) != length(truth)) {
      stop("Length of `cohort` must match `truth`.", call. = FALSE)
    }

    # Compute per-cohort values
    # Preserve factor levels when subsetting to avoid re-inferring classes
    values <- vapply(levels(cohort), function(level) {
      idx <- !is.na(cohort) & cohort == level
      if (!any(idx)) return(NA_real_)
      # Subset and preserve factor levels
      truth_subset <- factor(truth[idx], levels = truth_levels)
      base_loss(truth_subset, scores[idx], selected = selected,
                cutoff_prob = cutoff_prob,
                cutoff_strategy = cutoff_strategy,
                positive = positive)
    }, numeric(1))

    if (all(is.na(values))) return(NA_real_)
    aggregator(values, na.rm = TRUE)
  }
}

# Gap aggregator: max - min
.gap_aggregator <- function(values, na.rm = TRUE) {
  max(values, na.rm = na.rm) - min(values, na.rm = na.rm)
}

#' Minimum Cohort Sensitivity
#'
#' Computes sensitivity within each cohort and returns the minimum value to
#' capture worst-case performance.
#'
#' @inheritParams loss_sensitivity
#' @param cohort Factor indicating cohort membership.
#' @return Sensitivity of the weakest cohort.
#' @export
loss_min_cohort_sensitivity <- .make_cohort_aware_loss(
  base_loss = loss_sensitivity,
  aggregator = min,
  metric_name = "sensitivity"
)

#' Minimum Cohort Specificity
#'
#' Computes specificity within each cohort and returns the minimum value.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Specificity of the weakest cohort.
#' @export
loss_min_cohort_specificity <- .make_cohort_aware_loss(
  base_loss = loss_specificity,
  aggregator = min,
  metric_name = "specificity"
)

#' Cohort Sensitivity Range
#'
#' Difference between maximum and minimum cohort sensitivities. Smaller values
#' indicate more uniform transfer across cohorts.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Sensitivity range across cohorts.
#' @export
loss_cohort_sensitivity_gap <- .make_cohort_aware_loss(
  base_loss = loss_sensitivity,
  aggregator = .gap_aggregator,
  metric_name = "sensitivity",
  single_cohort_fallback = 0
)

#' Maximum Cohort Brier Score
#'
#' Computes the Brier score (mean squared error on probabilities) within each
#' cohort and returns the maximum, highlighting the worst calibrated cohort.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Maximum Brier score across cohorts.
#' @export
loss_max_cohort_brier <- function(truth, scores = NULL, selected = NULL,
                                  positive = "Yes", cohort = NULL) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute the Brier score.", call. = FALSE)
  }
  target <- as.numeric(truth == positive)
  cohort_vec <- factor(if (is.null(cohort)) rep("cohort_01", length(truth)) else cohort)
  if (length(cohort_vec) != length(truth)) {
    stop("Length of `cohort` must match `truth`.", call. = FALSE)
  }
  values <- vapply(levels(cohort_vec), function(level) {
    idx <- !is.na(cohort_vec) & cohort_vec == level
    if (!any(idx)) {
      return(NA_real_)
    }
    mean((scores[idx] - target[idx])^2)
  }, numeric(1))
  max(values, na.rm = TRUE)
}

#' Partial AUC in a Sensitivity-Focused Region
#'
#' Computes the partial area under the ROC curve restricted to sensitivity
#' values at or above a specified floor. Useful for "rule-out" diagnostic
#' contexts where high sensitivity is mandatory. Uses [pROC::auc()] with
#' `partial.auc.focus = "sens"`.
#'
#' @inheritParams loss_sensitivity
#' @param sens_floor Minimum sensitivity threshold defining the pAUC region
#'   (default 0.90). The partial AUC is computed only where sensitivity >= this
#'   value.
#' @param partial_auc_correct Logical; apply McClish correction to normalize
#'   partial AUC to 0-1 scale (default TRUE).
#' @return Partial AUC value, or `NA_real_` if computation fails.
#' @export
loss_pauc <- function(truth, scores = NULL, selected = NULL,
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

#' Maximum Mean Shift Across Cohorts
#'
#' Computes pairwise distances between cohort-specific mean expression vectors
#' for the selected features and returns the maximum distance. Encourages panels
#' whose selected biomarkers exhibit similar distributions across cohorts.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @param x Matrix of selected feature values (samples x features).
#' @return Maximum pairwise distance between cohort means.
#' @export
loss_max_cohort_mean_shift <- function(truth, scores = NULL, selected = NULL,
                                       cohort = NULL, x = NULL) {
  # Validate inputs - x must be provided with features
  if (is.null(x)) {
    stop("`x` (feature matrix) must be provided to compute cohort mean shift.", call. = FALSE)
  }
  if (ncol(x) == 0L) {
    stop("Feature matrix `x` has no columns. Cannot compute cohort mean shift.", call. = FALSE)
  }
  # Single cohort is valid - return 0 as there's no shift to measure
  if (is.null(cohort)) {
    return(0)
  }
  cohort <- factor(cohort)
  if (length(cohort) != nrow(x)) {
    stop("Length of `cohort` must match the number of rows in `x`.", call. = FALSE)
  }
  if (nlevels(cohort) <= 1L) {
    return(0)
  }
  means <- lapply(levels(cohort), function(level) {
    idx <- !is.na(cohort) & cohort == level
    if (!any(idx)) {
      return(NULL)
    }
    colMeans(x[idx, , drop = FALSE])
  })
  means <- Filter(Negate(is.null), means)
  if (length(means) <= 1L) {
    return(0)
  }
  mean_matrix <- do.call(rbind, means)
  if (nrow(mean_matrix) <= 1L) {
    return(0)
  }
  max(stats::dist(mean_matrix))
}

.loss_registry <- new.env(parent = emptyenv())

.register_default_loss <- function(name, fun, direction, label) {
  assign(name, list(fun = fun, direction = direction, label = label),
         envir = .loss_registry)
}

.register_default_loss(
  "sensitivity", loss_sensitivity, "maximize", "Sensitivity"
)
.register_default_loss(
  "specificity", loss_specificity, "maximize", "Specificity"
)
.register_default_loss(
  "auc", loss_auc, "maximize", "Area Under ROC Curve"
)
.register_default_loss(
  "f1", loss_f1, "maximize", "F1 Score"
)
.register_default_loss(
  "precision", loss_precision, "maximize", "Precision"
)
.register_default_loss(
  "npv", loss_npv, "maximize", "Negative Predictive Value"
)
.register_default_loss(
  "num_features", loss_num_features, "minimize", "Number of Features"
)
.register_default_loss(
  "balanced_accuracy", loss_balanced_accuracy, "maximize", "Balanced Accuracy"
)
.register_default_loss(
  "min_cohort_sensitivity", loss_min_cohort_sensitivity, "maximize",
  "Minimum Cohort Sensitivity"
)
.register_default_loss(
  "min_cohort_specificity", loss_min_cohort_specificity, "maximize",
  "Minimum Cohort Specificity"
)
.register_default_loss(
  "cohort_sensitivity_gap", loss_cohort_sensitivity_gap, "minimize",
  "Cohort Sensitivity Gap"
)
.register_default_loss(
  "max_cohort_brier", loss_max_cohort_brier, "minimize",
  "Maximum Cohort Brier Score"
)
.register_default_loss(
  "max_cohort_mean_shift", loss_max_cohort_mean_shift, "minimize",
  "Maximum Cohort Mean Shift"
)
.register_default_loss(
  "pauc", loss_pauc, "maximize", "Partial AUC (High Sensitivity)"
)

#' Register a loss function.
#'
#' @param name Unique identifier for the loss.
#' @param fun Function implementing the loss. Must accept at least the arguments
#'   `truth`, `scores`, and `selected`.
#' @param direction Either `"maximize"` or `"minimize"`.
#' @param label Human-readable name.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#' @export
register_loss_function <- function(name, fun,
                                   direction = c("maximize", "minimize"),
                                   label = name, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  direction <- match.arg(direction)
  if (!overwrite && exists(name, envir = .loss_registry, inherits = FALSE)) {
    stop(sprintf("Loss function '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(fun = fun, direction = direction, label = label),
         envir = .loss_registry)
  invisible(name)
}

#' List registered loss functions.
#'
#' @return Named list of loss registrations (`fun`, `direction`, `label`).
#' @export
loss_registry <- function() {
  if (length(ls(.loss_registry)) == 0L) {
    return(list())
  }
  mget(ls(.loss_registry), envir = .loss_registry)
}

#' Build objective descriptors from registered losses.
#'
#' Creates a list compatible with [define_objectives()] where each entry contains
#' the loss label, optimization direction, and an evaluation function with the
#' signature `(truth, estimate, selected = NULL)`. Additional parameters may be
#' supplied per loss via the `params` list (e.g., custom cutoffs).
#'
#' @param losses Character vector of registered loss names.
#' @param params Named list mapping loss names to argument lists.
#' @param directions Optional named character vector overriding the optimization
#'   direction for subsets of `losses`.
#' @return Named list of objective descriptors.
#' @export
build_objectives <- function(losses,
                             params = list(),
                             directions = NULL) {
  stopifnot(is.character(losses), length(losses) >= 1L)
  registry <- loss_registry()
  missing <- setdiff(losses, names(registry))
  if (length(missing)) {
    stop(sprintf("Loss function(s) not registered: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  objs <- lapply(losses, function(name) {
    entry <- registry[[name]]
    extras <- params[[name]]
    if (is.null(extras)) {
      extras <- list()
    } else if (!is.list(extras)) {
      extras <- as.list(extras)
    }
    direction <- entry$direction
    if (!is.null(directions) && !is.null(directions[[name]])) {
      direction <- match.arg(directions[[name]], c("maximize", "minimize"))
    }
    fun <- entry$fun
    fun_formals <- names(formals(fun))
    has_dots <- any(fun_formals == "...")
    wrapper <- function(truth, estimate, selected = NULL, ...) {
      dots <- list(...)
      extras_filtered <- extras
      if (!has_dots) {
        if (length(extras_filtered)) {
          extras_filtered <- extras_filtered[names(extras_filtered) %in% fun_formals]
        }
        if (length(dots)) {
          dots <- dots[names(dots) %in% fun_formals]
        }
      }
      args <- c(
        list(truth = truth, scores = estimate, selected = selected),
        extras_filtered,
        dots
      )
      do.call(fun, args)
    }
    list(
      label = entry$label,
      direction = direction,
      fun = wrapper
    )
  })
  names(objs) <- losses
  objs
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
