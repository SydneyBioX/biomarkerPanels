#' @useDynLib biomarkerPanels, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Compute pairwise difference between matrix columns
#' @param x A data matrix of size n times p. Where rows are observations and
#' columns are features.
#' @return A matrix of size n times (p choose 2), where each column is the
#' difference between two of the original columns.
#' @export
#' @examples
#' n = 1
#' p = 4
#' x = matrix(rep(1:p, n), nrow = n, ncol = p, byrow = TRUE)
#' colnames(x) = paste0("X", 1:p)
#' pairwise_col_diff(x)
pairwise_col_diff = function(x){
  if (is.null(colnames(x))) {
    stop("x must have column names")
  }

  p <- ncol(x)

  # Check for integer overflow in n_pairs calculation
  # p * (p-1) / 2 overflows 32-bit int when p > ~65,536
  if (p > 46340) {  # sqrt(.Machine$integer.max) ≈ 46340
    stop(
      "Too many features (", p, ") for pairwise aggregation. ",
      "Maximum supported is 46,340 features to avoid integer overflow. ",
      "Consider pre-filtering features before aggregation.",
      call. = FALSE
    )
  }

  # Warn about memory usage for large feature sets
  n_pairs <- (as.numeric(p) * (p - 1)) / 2
  n_samples <- nrow(x)
  estimated_mb <- (n_samples * n_pairs * 8) / (1024^2)
  if (estimated_mb > 1000) {
    warning(
      sprintf(
        "Pairwise aggregation will create %.0f pairs (%.1f GB). ",
        n_pairs, estimated_mb / 1024
      ),
      "Consider reducing the feature pool.",
      call. = FALSE
    )
  }

  # Sort columns
  col_order = order(colnames(x))
  x = x[, col_order, drop = FALSE]

  # Call Rcpp function
  .pairwise_col_diff_cpp(x, colnames(x))
}

#' Remove dominated solutions from a metric matrix
#'
#' Re-evaluation of Pareto-front candidates can shift metric values (e.g. due
#' to non-deterministic inner CV in glmnet), introducing dominated solutions.
#' This function returns the indices of non-dominated rows.
#'
#' @param metric_matrix Numeric matrix (rows = solutions, cols = objectives).
#' @param directions Named character vector with "maximize" or "minimize" per objective.
#' @return Integer vector of row indices that are non-dominated.
#' @keywords internal
.filter_dominated <- function(metric_matrix, directions) {
  n <- nrow(metric_matrix)
  if (n <= 1L) return(seq_len(n))

  dominated <- logical(n)
  for (i in seq_len(n)) {
    if (dominated[i]) next
    for (j in seq_len(n)) {
      if (i == j || dominated[j]) next
      # Check if j dominates i: j is >= on all objectives and > on at least one
      j_at_least <- TRUE
      j_strictly <- FALSE
      for (k in seq_len(ncol(metric_matrix))) {
        vi <- metric_matrix[i, k]
        vj <- metric_matrix[j, k]
        if (directions[k] == "maximize") {
          if (vj < vi) { j_at_least <- FALSE; break }
          if (vj > vi) j_strictly <- TRUE
        } else {
          if (vj > vi) { j_at_least <- FALSE; break }
          if (vj < vi) j_strictly <- TRUE
        }
      }
      if (j_at_least && j_strictly) {
        dominated[i] <- TRUE
        break
      }
    }
  }
  which(!dominated)
}
