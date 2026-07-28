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
