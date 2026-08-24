#' Panel Selection from NSGA Decision Vectors
#'
#' Maps continuous NSGA decision vectors onto discrete base-feature panels.
#' The adaptive strategy finds a natural breakpoint (largest gap) in the
#' sorted weights so panel sizes can vary across the Pareto front; the fixed
#' strategy thresholds weights directly.
#'
#' @name panel_selection
#' @keywords internal
NULL

.select_panel_indices <- function(decision_vec, feature_pool, max_features,
                                  min_features_required,
                                  selection_threshold = "adaptive") {
  n_pool <- length(feature_pool)
  if (identical(selection_threshold, "adaptive")) {
    ord <- order(decision_vec, feature_pool,
                 decreasing = c(TRUE, FALSE), method = "radix")
    sorted_weights <- decision_vec[ord]

    if (n_pool > max_features) {
      top_n <- min(max_features + 5L, n_pool)
      weight_diffs <- -diff(sorted_weights[seq_len(top_n)])
      search_range <- seq(min_features_required, top_n - 1L)
      if (length(search_range) > 0L) {
        gap_idx <- which.max(weight_diffs[search_range])
        natural_cutoff <- min_features_required + gap_idx
      } else {
        natural_cutoff <- min_features_required
      }
      n_selected <- min(natural_cutoff, max_features)
      n_selected <- max(n_selected, min_features_required)
    } else {
      above_median <- sum(decision_vec > 0.5)
      n_selected <- max(min_features_required, min(above_median, max_features))
    }
    selected_idx <- ord[seq_len(n_selected)]
  } else {
    threshold <- as.numeric(selection_threshold)
    above_threshold <- which(decision_vec > threshold)

    if (length(above_threshold) < min_features_required) {
      ord <- order(decision_vec, feature_pool,
                   decreasing = c(TRUE, FALSE), method = "radix")
      selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
    } else if (length(above_threshold) > max_features) {
      weights_above <- decision_vec[above_threshold]
      names_above <- feature_pool[above_threshold]
      ord <- order(weights_above, names_above,
                   decreasing = c(TRUE, FALSE), method = "radix")
      selected_idx <- above_threshold[ord[seq_len(max_features)]]
    } else {
      selected_idx <- above_threshold
    }
  }

  selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
}

.make_panel_selector <- function(feature_pool, max_features,
                                 min_features_required,
                                 selection_threshold = "adaptive") {
  force(feature_pool)
  force(max_features)
  force(min_features_required)
  force(selection_threshold)
  function(decision_vec) {
    selected_idx <- .select_panel_indices(
      decision_vec = decision_vec,
      feature_pool = feature_pool,
      max_features = max_features,
      min_features_required = min_features_required,
      selection_threshold = selection_threshold
    )
    base_features <- feature_pool[selected_idx]
    list(
      indices = selected_idx,
      base_features = base_features,
      key = .panel_key(base_features)
    )
  }
}
