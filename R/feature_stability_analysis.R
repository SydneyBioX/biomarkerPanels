#' Feature Selection Stability Analysis
#'
#' Quantifies how consistently features are chosen across the Pareto-optimal
#' solutions returned by [optimize_panel()], via inclusion frequencies and
#' pairwise Jaccard similarity. Complements the selection helpers in
#' `panel_summarization.R`.
#'
#' @name feature_stability_analysis
NULL

#' Analyze Feature Selection Stability Across Pareto Solutions
#'
#' Computes feature inclusion frequencies and pairwise Jaccard similarity
#' between all Pareto-optimal solutions. Accepts either an
#' [`OptimizationResult`] (preferred — contains all Pareto solutions) or a
#' [`BiomarkerPanelResult`] (legacy path).
#'
#' @param result An [`OptimizationResult`] or [`BiomarkerPanelResult`].
#' @param feature_type One of `"features"` (default) or `"base_features"`.
#'   When `"base_features"`, counts original gene names rather than transformed
#'   feature names (e.g. pairwise ratios). Only supported for
#'   `OptimizationResult` input.
#' @return A list with class `"FeatureStabilityResult"` containing:
#'   \describe{
#'     \item{frequencies}{Data frame with columns `feature`, `count`, `proportion`
#'       ordered by descending count.}
#'     \item{jaccard_matrix}{Symmetric matrix of pairwise Jaccard similarity
#'       coefficients between solutions. Row/column names are `solution_id` values.}
#'     \item{n_solutions}{Integer count of Pareto-optimal solutions analyzed.}
#'     \item{mean_jaccard}{Mean pairwise Jaccard similarity (excluding diagonal).}
#'   }
#' @export
#' @seealso [compute_inclusion_frequencies()], [plot_feature_stability()]
#' @examples
#' \dontrun{
#' # From an OptimizationResult (preferred):
#' stability <- analyze_feature_stability(opt_result)
#' stability$frequencies
#' plot_feature_stability(stability)
#'
#' # Gene-level frequencies with pairwise transforms:
#' stability_genes <- analyze_feature_stability(opt_result,
#'                                              feature_type = "base_features")
#' }
analyze_feature_stability <- function(result,
                                      feature_type = c("features",
                                                       "base_features")) {
  feature_type <- match.arg(feature_type)

  if (inherits(result, "OptimizationResult")) {
    solution_features <- .extract_optimization_features(result, feature_type)
  } else if (inherits(result, "BiomarkerPanelResult")) {
    if (feature_type == "base_features") {
      stop("feature_type = \"base_features\" is only supported for ",
           "OptimizationResult input.", call. = FALSE)
    }
    solution_features <- .extract_solution_features(result)
  } else {
    stop("`result` must be an OptimizationResult or BiomarkerPanelResult.",
         call. = FALSE)
  }

  n_solutions <- length(solution_features)

  # Handle empty result
  if (n_solutions == 0L) {
    return(structure(
      list(
        frequencies = data.frame(
          feature = character(),
          count = integer(),
          proportion = numeric(),
          stringsAsFactors = FALSE
        ),
        jaccard_matrix = matrix(nrow = 0L, ncol = 0L),
        n_solutions = 0L,
        mean_jaccard = NA_real_
      ),
      class = "FeatureStabilityResult"
    ))
  }

  # Compute frequencies
  frequencies <- compute_inclusion_frequencies(list(result),
                                               feature_type = feature_type)

  # Compute Jaccard similarity matrix
  solution_ids <- names(solution_features)
  jaccard_matrix <- matrix(
    NA_real_,
    nrow = n_solutions,
    ncol = n_solutions,
    dimnames = list(solution_ids, solution_ids)
  )

  for (i in seq_len(n_solutions)) {
    for (j in seq_len(i)) {
      set_i <- solution_features[[i]]
      set_j <- solution_features[[j]]

      if (length(set_i) == 0L && length(set_j) == 0L) {
        jaccard <- 1.0
      } else if (length(set_i) == 0L || length(set_j) == 0L) {
        jaccard <- 0.0
      } else {
        intersection_size <- length(intersect(set_i, set_j))
        union_size <- length(union(set_i, set_j))
        jaccard <- intersection_size / union_size
      }

      jaccard_matrix[i, j] <- jaccard
      jaccard_matrix[j, i] <- jaccard
    }
  }

  # Calculate mean Jaccard (off-diagonal)
  if (n_solutions > 1L) {
    off_diag <- jaccard_matrix[lower.tri(jaccard_matrix)]
    mean_jaccard <- mean(off_diag)
  } else {
    mean_jaccard <- NA_real_
  }

  structure(
    list(
      frequencies = frequencies,
      jaccard_matrix = jaccard_matrix,
      n_solutions = n_solutions,
      mean_jaccard = mean_jaccard
    ),
    class = "FeatureStabilityResult"
  )
}

#' Extract per-solution feature lists from an OptimizationResult
#'
#' @param opt An `OptimizationResult`.
#' @param feature_type One of `"features"` or `"base_features"`.
#' @return Named list of character vectors (one per solution).
#' @keywords internal
.extract_optimization_features <- function(opt, feature_type = "features") {
  sol_df <- opt@solutions
  if (!nrow(sol_df)) return(list())

  feat_list <- sol_df[[feature_type]]
  names(feat_list) <- as.character(sol_df$solution_id)
  lapply(feat_list, unique)
}
