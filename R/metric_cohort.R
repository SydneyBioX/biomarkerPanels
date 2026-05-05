#' Cohort-Aware Metric Functions
#'
#' Metric functions that evaluate performance across multiple cohorts, capturing
#' worst-case performance, cohort gaps, and transferability metrics.
#'
#' @name metric_cohort
NULL

# -----------------------------------------------------------------------------
# Factory for Cohort-Aware Metric Functions
# -----------------------------------------------------------------------------
# This factory pattern eliminates code duplication across cohort-aware metrics.
# Each cohort-aware metric applies a base metric to each cohort and aggregates.

#' Create a Cohort-Aware Metric Function
#'
#' Factory function that wraps a base metric function to compute per-cohort
#' values and aggregate them. Reduces code duplication for cohort-aware metrics.
#'
#' @param base_metric The base metric function (e.g., metric_sensitivity).
#' @param aggregator Function to aggregate cohort values (e.g., min, max).
#' @param metric_name Name of the metric for error messages.
#' @param single_cohort_fallback Value to return for single cohort with gap
#'   aggregation (default NULL uses base metric result).
#' @return A cohort-aware metric function.
#' @keywords internal
.make_cohort_aware_metric <- function(base_metric, aggregator, metric_name,
                                     single_cohort_fallback = NULL) {
  function(truth, scores = NULL, selected = NULL,
           positive = "Yes",
           cohort = NULL, ...) {
    # Standardize truth once at the start to get consistent factor levels
    truth <- ensure_binary_response(truth)
    truth_levels <- levels(truth)

    if (is.null(scores)) {
      stop("`scores` must be supplied to compute ", metric_name, ".", call. = FALSE)
    }

    # Fall back to base metric when no cohort provided
    if (is.null(cohort)) {
      base_result <- base_metric(truth, scores, selected,
                               positive = positive, ...)
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
      # Return NA when a cohort lacks both classes (e.g., all positive),
      # but let genuine errors propagate
      pos_n <- sum(truth_subset == positive)
      neg_n <- sum(truth_subset != positive)
      if (pos_n == 0L || neg_n == 0L) return(NA_real_)
      base_metric(truth_subset, scores[idx], selected = selected,
                  positive = positive)
    }, numeric(1))

    if (all(is.na(values))) return(NA_real_)
    aggregator(values, na.rm = TRUE)
  }
}

#' Gap Aggregator for Cohort Metrics
#'
#' Computes the range (max - min) of values across cohorts. Used internally
#' by cohort gap metric functions.
#'
#' @param values Numeric vector of per-cohort metric values.
#' @param na.rm Logical; whether to remove NA values before computation.
#' @return The difference between max and min values.
#' @keywords internal
.gap_aggregator <- function(values, na.rm = TRUE) {
  max(values, na.rm = na.rm) - min(values, na.rm = na.rm)
}

#' Variance Aggregator for Cohort Metrics
#'
#' Computes the sample variance of values across cohorts. Returns 0 when
#' fewer than 2 non-NA values are available.
#'
#' @param values Numeric vector of per-cohort metric values.
#' @param na.rm Logical; whether to remove NA values before computation.
#' @return Variance of the values, or 0 for single-cohort input.
#' @keywords internal
.variance_aggregator <- function(values, na.rm = TRUE) {
  if (na.rm) values <- values[!is.na(values)]
  if (length(values) <= 1L) return(0)
  stats::var(values)
}

#' Minimum Cohort AUC
#'
#' Computes AUC within each cohort and returns the minimum value to capture
#' worst-case discrimination. Unlike cutoff-dependent cohort metrics, this
#' metric is threshold-free.
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param positive Label treated as the positive ("event") class.
#' @param cohort Factor indicating cohort membership.
#' @param ... Additional arguments forwarded to [metric_auc()].
#' @return AUC of the weakest cohort.
#' @note Per-cohort AUC is noisy when cohorts contain fewer than ~20 samples
#'   of each class. Consider using `metric_auc` as the primary objective and
#'   this metric for monitoring.
#' @export
metric_min_cohort_auc <- .make_cohort_aware_metric(
  base_metric = metric_auc,
  aggregator = min,
  metric_name = "AUC"
)

#' Cohort AUC Range
#'
#' Difference between maximum and minimum per-cohort AUC values. Smaller
#' values indicate more uniform discrimination across cohorts.
#'
#' @inheritParams metric_min_cohort_auc
#' @return AUC range across cohorts.
#' @note Per-cohort AUC is noisy when cohorts contain fewer than ~20 samples
#'   of each class.
#' @export
metric_cohort_auc_gap <- .make_cohort_aware_metric(
  base_metric = metric_auc,
  aggregator = .gap_aggregator,
  metric_name = "AUC",
  single_cohort_fallback = 0
)

#' Cohort AUC Variance
#'
#' Sample variance of per-cohort AUC values. Penalises uneven discrimination
#' across cohorts without being dominated by extreme cohorts the way the gap
#' metric can be.
#'
#' @inheritParams metric_min_cohort_auc
#' @return Variance of per-cohort AUC values.
#' @note Per-cohort AUC is noisy when cohorts contain fewer than ~20 samples
#'   of each class.
#' @export
metric_cohort_auc_var <- .make_cohort_aware_metric(
  base_metric = metric_auc,
  aggregator = .variance_aggregator,
  metric_name = "AUC",
  single_cohort_fallback = 0
)

#' Maximum Cohort Brier Score
#'
#' Computes the Brier score (mean squared error on probabilities) within each
#' cohort and returns the maximum, highlighting the worst calibrated cohort.
#'
#' @inheritParams metric_min_cohort_auc
#' @return Maximum Brier score across cohorts.
#' @export
metric_max_cohort_brier <- function(truth, scores = NULL, selected = NULL,
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
