#' @useDynLib biomarkerPanels, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Rcpp implementations (compiled on package load)
# These are internal functions used by the exported R wrappers

#' @title Compute pairwise difference between matrix columns (C++ implementation)
#' @param x A numeric matrix
#' @param col_names Character vector of column names
#' @return A matrix of pairwise differences
#' @keywords internal
# This will be compiled from src/pairwise.cpp in your package

#' @title Compute pairwise ratios to a specific feature (C++ implementation)
#' @param x A numeric matrix
#' @param feature_col The feature column to divide by
#' @param feature_name Name of the feature
#' @param other_names Names of other columns
#' @return A matrix of pairwise ratios
#' @keywords internal
# This will be compiled from src/pairwise.cpp in your package

#' @title Compute pairwise difference between matrix columns
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
  
  # Sort columns
  col_order = order(colnames(x))
  x = x[, col_order, drop = FALSE]
  
  # Call Rcpp function
  .pairwise_col_diff_cpp(x, colnames(x))
}
