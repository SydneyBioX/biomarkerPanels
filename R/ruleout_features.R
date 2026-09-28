#' Rule-Out Feature Selection
#'
#' Functions for ranking features by the univariate specificity they achieve at
#' a fixed high-sensitivity operating point, the quantity that matters for a
#' rule-out diagnostic.
#'
#' @name ruleout_features
NULL

#' Select Rule-Out Features by Specificity at a Sensitivity Target
#'
#' Ranks each feature by the specificity it achieves as a single-marker
#' classifier when its threshold is set to capture at least
#' `target_sensitivity` of cases. Both directions are tried (high values call
#' positive, or low values call positive) and the better one is kept, so the
#' score is direction-agnostic. With multiple cohorts the per-cohort
#' specificities are combined by `cohort_aggregation` (the default `"min"`
#' rewards features that hold up in every cohort) *within* each direction, and
#' only then is the better direction kept. The direction is therefore shared
#' by every cohort: a feature that is high in cases in one cohort and low in
#' another cannot score well in both, which matters for ratio transforms, where
#' a direction flip cancels rather than adds signal.
#'
#' This differs from AUC-based ranking (e.g. [select_discriminative_features()])
#' in that only the high-sensitivity tail of the ROC curve counts: a feature
#' whose cases are mostly well separated but with a minority buried inside the
#' control distribution has a good AUC yet a poor rule-out specificity. It is
#' the univariate, per-feature analogue of the
#' `"specificity_at_sensitivity"` objective used to score fitted panels during
#' optimisation.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Provide a list when
#'   `x_list` is a list.
#' @param n_features Maximum number of features to return (default `50`).
#' @param target_sensitivity Minimum fraction of cases a feature's threshold
#'   must capture (default `0.90`). Must lie in `(0, 1]`.
#' @param cohort_aggregation How per-cohort specificities are combined into a
#'   single ranking score: `"min"` (default, worst cohort) or `"mean"`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to analysis.
#' @return Character vector of feature identifiers ordered by decreasing
#'   specificity at the sensitivity target. A feature whose score is undefined
#'   in any cohort (no non-missing case values there) is dropped.
#' @details Thresholds are observed case values (order statistics): with `n`
#'   non-missing cases the smallest `k = ceiling(n * target_sensitivity)` cases
#'   that still meet the target are called positive, so the positive side
#'   always contains at least `target_sensitivity` of the cases and never more
#'   than needed. Missing values are ignored when locating thresholds and
#'   computing specificities. Every cohort must contain at least one case and
#'   one control.
#' @seealso [select_de_features()], [select_discriminative_features()],
#'   [metric_specificity_at_sensitivity()]
#' @export
select_ruleout_features <- function(x_list,
                                    y_list,
                                    n_features = 50L,
                                    target_sensitivity = 0.90,
                                    cohort_aggregation = c("min", "mean"),
                                    assay = NULL) {
  n_features <- .validate_positive_integer(n_features, "n_features")
  target_sensitivity <- .validate_probability(
    target_sensitivity, "target_sensitivity", bounds = "closed"
  )
  if (target_sensitivity == 0) {
    stop("`target_sensitivity` must be greater than 0.", call. = FALSE)
  }
  cohort_aggregation <- match.arg(cohort_aggregation)

  prepared <- .prepare_selection_inputs(x_list, y_list, assay = assay)
  matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  feature_names <- prepared$feature_names

  spec_by_cohort <- lapply(seq_along(matrices), function(k) {
    y <- responses[[k]]
    if (!any(y == 1L) || !any(y == 0L)) {
      stop(
        "Cohort `", cohort_names[[k]], "` must contain at least one case and ",
        "one control to compute specificity at sensitivity.",
        call. = FALSE
      )
    }
    .spec_at_sensitivity_by_feature(matrices[[k]], y, target_sensitivity)
  })

  # Aggregate across cohorts within each direction, then keep the better
  # direction, so one direction has to hold in every cohort.
  aggregate_cohorts <- function(direction) {
    spec_matrix <- vapply(spec_by_cohort, function(spec) spec[, direction],
                          numeric(length(feature_names)))
    if (is.null(dim(spec_matrix))) {
      spec_matrix <- matrix(spec_matrix, ncol = length(matrices))
    }
    switch(
      cohort_aggregation,
      min = apply(spec_matrix, 1, min),
      mean = rowMeans(spec_matrix)
    )
  }
  scores <- pmax(aggregate_cohorts("up"), aggregate_cohorts("down"))
  names(scores) <- feature_names

  scores <- scores[!is.na(scores)]
  if (!length(scores)) {
    return(character())
  }

  ordered <- names(scores)[order(scores, decreasing = TRUE)]
  head(ordered, n = min(n_features, length(ordered)))
}

#' Per-Feature Specificity at a Sensitivity Target
#'
#' For each column, sets the threshold at the observed case value that captures
#' the smallest number of cases still meeting `target_sensitivity`, and reports
#' the fraction of controls falling on the negative side. Both call directions
#' are evaluated and returned separately, so callers can hold the direction
#' fixed across cohorts.
#'
#' With `n` non-missing case values, `k = ceiling(n * target_sensitivity)`
#' cases must be called positive. The "up" threshold is the `(n - k + 1)`-th
#' smallest case value (positive if `x >= t`), the "down" threshold is the
#' `k`-th smallest (positive if `x <= t`). Using order statistics rather than an
#' interpolated quantile guarantees the positive side holds exactly `k` cases
#' (at least `target_sensitivity`), and treats both directions symmetrically.
#'
#' @param x Numeric matrix (samples x features).
#' @param y Integer 0/1 vector of labels (1 = case).
#' @param target_sensitivity Sensitivity floor in `(0, 1]`.
#' @return Numeric matrix of specificities with one row per column of `x` and
#'   columns `up` (positive if `x >= t`) and `down` (positive if `x <= t`).
#'   `NA` where a column has no non-missing case or control values.
#' @noRd
.spec_at_sensitivity_by_feature <- function(x, y, target_sensitivity) {
  x_pos <- x[y == 1L, , drop = FALSE]
  x_neg <- x[y == 0L, , drop = FALSE]

  thresholds <- apply(x_pos, 2, .ruleout_thresholds, target = target_sensitivity)
  t_up <- thresholds[1L, ]
  t_down <- thresholds[2L, ]

  # "Up": positive if x >= t_up, so a control is correctly negative when it
  # falls strictly below the threshold. "Down": positive if x <= t_down.
  spec_up <- colMeans(sweep(x_neg, 2, t_up, "<"), na.rm = TRUE)
  spec_down <- colMeans(sweep(x_neg, 2, t_down, ">"), na.rm = TRUE)

  spec <- cbind(up = spec_up, down = spec_down)
  spec[is.na(t_up) | is.nan(spec)] <- NA_real_
  spec
}

#' Rule-Out Thresholds for One Feature
#'
#' @param v Numeric vector of case values (may contain `NA`).
#' @param target Sensitivity floor in `(0, 1]`.
#' @return Length-2 numeric vector `c(up, down)` of observed case values, or
#'   `NA`s when no non-missing case values exist.
#' @noRd
.ruleout_thresholds <- function(v, target) {
  v <- sort(v[!is.na(v)])
  n <- length(v)
  if (!n) {
    return(c(NA_real_, NA_real_))
  }
  # Small tolerance guards against 40 * 0.9 evaluating to 36.000000000000007.
  k <- max(1L, min(n, as.integer(ceiling(n * target - 1e-8))))
  c(v[n - k + 1L], v[k])
}
