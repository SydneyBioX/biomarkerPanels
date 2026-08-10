#' Deprecated Functions
#'
#' Functions retained for backward compatibility. Each forwards to its
#' replacement and signals a deprecation warning. They will be removed in a
#' future release.
#'
#' @name biomarkerPanels-deprecated
NULL

#' @describeIn biomarkerPanels-deprecated Renamed to [select_de_features()] so
#'   that every feature-selection entry point shares the `select_*_features()`
#'   convention.
#' @inheritParams select_de_features
#' @return Character vector of feature identifiers ordered by significance.
#' @export
get_top_de_features <- function(x_list,
                                y_list,
                                contrast = NULL,
                                n_features = 50L,
                                combination_method = c("OSP", "Stouffer", "Fisher", "maxP"),
                                assay = NULL) {
  .Deprecated("select_de_features")
  select_de_features(
    x_list = x_list,
    y_list = y_list,
    contrast = contrast,
    n_features = n_features,
    combination_method = match.arg(combination_method),
    assay = assay
  )
}

#' @describeIn biomarkerPanels-deprecated Renamed to
#'   [select_batch_associated_features()]. The old name described the genes'
#'   downstream role as log-ratio denominators rather than the selection
#'   criterion (association with batch, not biology). The argument order is
#'   unchanged, so positional calls still work; `n_denominators` is now
#'   `n_features` and the returned `denominators` element is now `features`.
#' @inheritParams select_batch_associated_features
#' @param n_denominators Number of genes to select. Renamed to `n_features` in
#'   [select_batch_associated_features()].
#' @return A list as returned by [select_batch_associated_features()], with the
#'   `features` element additionally exposed under its former name
#'   `denominators`.
#' @export
select_denominator_features <- function(
    x,
    cohort,
    y = NULL,
    n_denominators = 20L,
    n_pcs = 50L,
    batch_pvalue_threshold = 0.10,
    biology_pvalue_threshold = 0.20,
    biology_penalty_quantile = 0.5) {
  .Deprecated("select_batch_associated_features")
  result <- select_batch_associated_features(
    x = x,
    cohort = cohort,
    y = y,
    n_features = n_denominators,
    n_pcs = n_pcs,
    batch_pvalue_threshold = batch_pvalue_threshold,
    biology_pvalue_threshold = biology_pvalue_threshold,
    biology_penalty_quantile = biology_penalty_quantile
  )
  c(list(denominators = result$features), result)
}
