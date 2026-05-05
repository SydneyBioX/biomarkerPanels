#' Cohort Transferability Metric Functions
#'
#' Metric functions that detect cohort-specific contamination in model scores.
#' Unlike the performance-based cohort metrics in \code{\link{metric_cohort}},
#' these metrics measure whether score distributions are systematically shifted
#' by cohort membership within the same class — a direct proxy for batch-effect
#' leakage.
#'
#' @name metric_cohort_transferability
NULL

# -----------------------------------------------------------------------------
# Internal Helpers
# -----------------------------------------------------------------------------

#' 1-D Wasserstein Distance via Quantile Approximation
#'
#' Computes the W1 (earth mover's) distance between two univariate samples
#' using a quantile-based approximation.
#'
#' @param a,b Numeric vectors of observations.
#' @param n Number of quantile evaluation points (default 200).
#' @return W1 distance, or \code{NA_real_} if either sample has fewer than 10
#'   finite values.
#' @keywords internal
.w1 <- function(a, b, n = 200L) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 10L || length(b) < 10L) return(NA_real_)
  probs <- seq(0, 1, length.out = n)
  qa <- as.numeric(quantile(a, probs, type = 8, na.rm = TRUE))
  qb <- as.numeric(quantile(b, probs, type = 8, na.rm = TRUE))
  mean(abs(qa - qb), na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Cohort Leakage (R-squared)
# -----------------------------------------------------------------------------

#' Cohort Leakage
#'
#' Measures how much of the within-class score variance is explained by cohort
#' membership, using adjusted R-squared from a one-way ANOVA
#' (\code{lm(scores ~ cohort)}). For each class (positives and negatives), the
#' adjusted R-squared is computed separately; the metric returns the maximum
#' across both classes.
#'
#' A value near 0 indicates scores are independent of cohort within each class
#' (good transferability). A value near 1 indicates scores are almost entirely
#' determined by cohort (batch-effect leakage).
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param positive Label treated as the positive ("event") class.
#' @param cohort Factor indicating cohort membership.
#' @param ... Additional arguments (ignored).
#' @return Maximum adjusted R-squared across classes, clamped to \code{[0, 1]}.
#'   Returns 0 when \code{cohort} is \code{NULL} or contains a single level.
#' @export
metric_cohort_leakage <- function(truth, scores = NULL, selected = NULL,
                                  positive = "Yes", cohort = NULL, ...) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute cohort leakage.", call. = FALSE)
  }
  if (is.null(cohort)) return(0)

  cohort <- droplevels(factor(cohort))
  y <- truth == positive

  leakage_one_class <- function(idx) {
    c_sub <- droplevels(cohort[idx])
    s_sub <- scores[idx]
    if (nlevels(c_sub) < 2L || length(s_sub) < 10L) return(0)
    adj_r2 <- summary(lm(s_sub ~ c_sub))$adj.r.squared
    max(0, adj_r2)
  }

  leak_case    <- leakage_one_class(which(y))
  leak_control <- leakage_one_class(which(!y))
  max(leak_case, leak_control, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Class-Conditional Score Shift (Wasserstein)
# -----------------------------------------------------------------------------

#' Class-Conditional Score Shift
#'
#' Measures the maximum pairwise Wasserstein (W1) distance between cohort-level
#' score distributions within each class. For each class, all cohort pairs are
#' compared; the metric returns the worst-case shift across both classes.
#'
#' Unlike \code{\link{metric_cohort_leakage}}, this metric is fully
#' non-parametric and captures any distributional shift (location, scale, or
#' shape), not only mean differences.
#'
#' @inheritParams metric_cohort_leakage
#' @return Maximum pairwise W1 distance across classes, or 0 when \code{cohort}
#'   is \code{NULL} or contains a single level. Returns \code{NA_real_} when all
#'   cohort subgroups are too small (< 10 samples) for reliable comparison.
#' @export
metric_conditional_score_shift <- function(truth, scores = NULL, selected = NULL,
                                           positive = "Yes", cohort = NULL, ...) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute conditional score shift.",
         call. = FALSE)
  }
  if (is.null(cohort)) return(0)

  cohort <- droplevels(factor(cohort))
  y <- truth == positive

  shift_one_class <- function(idx) {
    c_sub <- droplevels(cohort[idx])
    s_sub <- scores[idx]
    lev <- levels(c_sub)
    if (length(lev) < 2L) return(0)

    d <- combn(lev, 2, function(p) {
      .w1(s_sub[c_sub == p[1]], s_sub[c_sub == p[2]])
    })
    if (all(is.na(d))) return(NA_real_)
    max(d, na.rm = TRUE)
  }

  shift_case    <- shift_one_class(which(y))
  shift_control <- shift_one_class(which(!y))

  vals <- c(shift_case, shift_control)
  if (all(is.na(vals))) return(NA_real_)
  max(vals, na.rm = TRUE)
}
