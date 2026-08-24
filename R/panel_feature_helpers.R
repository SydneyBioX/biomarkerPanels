#' Panel Feature Preparation Helpers
#'
#' Internal helpers that centralize how a fitted panel's selected features are
#' prepared from a raw feature matrix before scoring. Keeping this in one place
#' ensures [evaluate_panel()], [evaluate_panel_by_cohort()], and
#' [calibrate_panel()] apply identical base-feature validation and feature
#' transforms, which is statistically load-bearing (the model was fit on these
#' exact transformed columns).
#'
#' @name panel_feature_helpers
#' @noRd
NULL

#' Validate and Transform Selected Base Features
#'
#' Checks that a panel's base features are present in `x_raw`, subsets to them,
#' and applies the panel's feature transform on the fly. Callers are responsible
#' for resolving `base_features` / `feature_transform` from the panel (some
#' allow an explicit transform override) and for any downstream subsetting to a
#' model's expected column subset.
#'
#' @param x_raw Raw feature matrix (samples x base features) with column names.
#' @param base_features Character vector of base feature names to select.
#' @param feature_transform Registered transform name, or `"none"`.
#' @param context Human-readable location for the not-found error message
#'   (e.g. `"validation data"`, `"held-out data"`).
#' @return List with `x_selected` (transformed matrix) and `selected`
#'   (transformed feature names).
#' @noRd
.prepare_scoring_matrix <- function(x_raw, base_features, feature_transform,
                                    context = "data") {
  if (!all(base_features %in% colnames(x_raw))) {
    missing <- setdiff(base_features, colnames(x_raw))
    stop("Base feature(s) not found in ", context, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  x_base <- x_raw[, base_features, drop = FALSE]

  if (feature_transform != "none" && length(base_features) >= 2L) {
    x_selected <- .apply_feature_transform_single(x_base, feature_transform)
    selected <- colnames(x_selected)
  } else {
    x_selected <- x_base
    selected <- base_features
  }

  list(x_selected = x_selected, selected = selected)
}
