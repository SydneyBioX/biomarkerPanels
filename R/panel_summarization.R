#' Pareto solution utilities and selection strategies
#'
#' These helpers support downstream triage of Pareto-optimal solutions produced
#' by [`optimize_panel()`]: inspecting solutions, computing biomarker inclusion
#' frequencies across solutions, and assembling a panel from the most
#' frequently selected biomarkers.

#' Summarize Pareto-Optimal Solutions
#'
#' Returns a wide-format data frame with one row per Pareto solution,
#' showing all objective values and feature count for easy inspection.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @return A data frame with columns:
#'   \describe{
#'     \item{solution_id}{Integer ID of the solution.}
#'     \item{n_features}{Number of base features (individual genes) in the solution, not ratio pairs.}
#'     \item{...}{One column per objective with numeric values.}
#'   }
#' @export
#' @seealso [optimize_panel()], [fit_panel()], [get_solution_features()]
#' @examples
#' \dontrun{
#' opt <- optimize_panel(x, y, objectives = define_objectives())
#' summarize_solutions(opt)
#' }
summarize_solutions <- function(optimization_result) {
  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult from optimize_panel().",
         call. = FALSE)
  }

  solutions_df <- optimization_result@solutions
  if (nrow(solutions_df) == 0L) {
    return(data.frame(
      solution_id = integer(),
      n_features = integer(),
      stringsAsFactors = FALSE
    ))
  }

  # Count base features (genes to measure), not transformed features (ratios)
  feature_col <- if ("base_features" %in% names(solutions_df)) "base_features" else "features"
  n_features <- vapply(solutions_df[[feature_col]], length, integer(1))

  # Get objective columns (everything except solution_id and features)
  objective_cols <- setdiff(names(solutions_df), c("solution_id", "features", "base_features"))

  # Build result data frame
  result <- data.frame(
    solution_id = solutions_df$solution_id,
    n_features = n_features,
    stringsAsFactors = FALSE
  )

  # Add objective columns (num_features from the metric also counts base features now)
  for (col in objective_cols) {
    result[[col]] <- solutions_df[[col]]
  }

  result
}

#' Extract solution-level feature sets from a `BiomarkerPanelResult`.
#'
#' @param panel_result A [`BiomarkerPanelResult`].
#' @return A list where each element corresponds to a Pareto solution and
#'   contains the unique biomarker identifiers selected in that solution. The
#'   list is named by `solution_id`.
#' @noRd
.extract_solution_features <- function(panel_result) {
  stopifnot(inherits(panel_result, "BiomarkerPanelResult"))
  objectives <- panel_result@objectives
  if (!nrow(objectives)) {
    return(list())
  }

  split_idx <- split(seq_len(nrow(objectives)), objectives$solution_id)
  lapply(split_idx, function(idx) {
    feats <- objectives$features[idx]
    unique(unlist(feats, use.names = FALSE))
  })
}

#' Compute biomarker inclusion frequencies across optimisation runs.
#'
#' @param panels List of `OptimizationResult`, [`BiomarkerPanelResult`]
#'   objects, character vectors of biomarkers, or a mixture of these. A single
#'   (unwrapped) `OptimizationResult` is also accepted.
#' @param feature_type One of `"features"` (default) or `"base_features"`.
#'   Controls whether transformed feature names (e.g. `"A--B"` for pairwise
#'   ratios) or base feature names (genes) are counted. Only relevant for
#'   `OptimizationResult` entries; `BiomarkerPanelResult` and character vector
#'   entries always use their stored feature names.
#' @return Data frame with columns `feature`, `count`, and `proportion`.
#' @export
compute_inclusion_frequencies <- function(panels,
                                          feature_type = c("features",
                                                           "base_features")) {
  feature_type <- match.arg(feature_type)

  # Accept a single unwrapped OptimizationResult
  if (inherits(panels, "OptimizationResult")) {
    panels <- list(panels)
  }

  if (is.null(panels) || !length(panels)) {
    return(data.frame(
      feature = character(),
      count = integer(),
      proportion = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  collected <- lapply(panels, function(p) {
    if (inherits(p, "OptimizationResult")) {
      sol_df <- p@solutions
      if (!nrow(sol_df)) {
        return(list(features = character(), solutions = 0L))
      }
      feat_list <- sol_df[[feature_type]]
      feat_list <- lapply(feat_list, unique)
      list(
        features = unlist(feat_list, use.names = FALSE),
        solutions = nrow(sol_df)
      )
    } else if (inherits(p, "BiomarkerPanelResult")) {
      sol <- .extract_solution_features(p)
      sol <- lapply(sol, unique)
      list(
        features = unlist(sol, use.names = FALSE),
        solutions = length(sol)
      )
    } else if (is.character(p)) {
      list(
        features = unique(p),
        solutions = 1L
      )
    } else {
      stop("Unsupported panel entry: must be OptimizationResult, ",
           "BiomarkerPanelResult, or character vector.", call. = FALSE)
    }
  })

  flattened <- unlist(lapply(collected, `[[`, "features"), use.names = FALSE)

  if (!length(flattened)) {
    return(data.frame(
      feature = character(),
      count = integer(),
      proportion = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  freq_table <- sort(table(flattened), decreasing = TRUE)
  total_solutions <- sum(vapply(collected, `[[`, integer(1), "solutions"))
  df <- data.frame(
    feature = names(freq_table),
    count = as.integer(freq_table),
    stringsAsFactors = FALSE
  )
  df$proportion <- df$count / total_solutions
  df
}

#' Assemble a panel by inclusion frequency.
#'
#' @param frequency_df Data frame returned by
#'   [compute_inclusion_frequencies()].
#' @param panel_size Target number of biomarkers.
#' @return Character vector of selected biomarker identifiers (length up to
#'   `panel_size`).
#' @export
select_panel_inclusion_frequency <- function(frequency_df, panel_size = 4L) {
  if (!nrow(frequency_df)) {
    return(character())
  }
  panel_size <- .validate_positive_integer(panel_size, "panel_size")
  head(frequency_df$feature, n = min(panel_size, nrow(frequency_df)))
}

