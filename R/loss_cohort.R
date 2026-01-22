#' Cohort-Aware Loss Functions
#'
#' Loss functions that evaluate performance across multiple cohorts, capturing
#' worst-case performance, cohort gaps, and transferability metrics.
#'
#' @name loss_cohort
NULL

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

#' Gap Aggregator for Cohort Metrics
#'
#' Computes the range (max - min) of values across cohorts. Used internally
#' by cohort gap loss functions.
#'
#' @param values Numeric vector of per-cohort metric values.
#' @param na.rm Logical; whether to remove NA values before computation.
#' @return The difference between max and min values.
#' @keywords internal
.gap_aggregator <- function(values, na.rm = TRUE) {
  max(values, na.rm = na.rm) - min(values, na.rm = na.rm)
}

#' Minimum Cohort Sensitivity
#'
#' Computes sensitivity within each cohort and returns the minimum value to
#' capture worst-case performance.
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param cutoff_prob Classification probability cutoff applied to `scores`.
#' @param cutoff_strategy Strategy for computing cutoff.
#' @param positive Label treated as the positive ("event") class.
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
