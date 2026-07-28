#' Dual-Stream Ratio Feature Construction
#'
#' Functions for constructing ratio features from separate numerator (signal)
#' and denominator (batch-proxy) gene sets.
#'
#' @name dual_stream_ratios
NULL

#' Construct Dual-Stream Ratio Features
#'
#' Creates pairwise ratio features by combining numerator genes (signal genes,
#' typically differentially expressed) with denominator genes (batch-proxy genes
#' from [select_denominator_features()]).
#'
#' For log-transformed expression data, log-ratios are computed as:
#' `log(Numerator) - log(Denominator) = log(Numerator/Denominator)`
#'
#' This normalization strategy exploits the fact that batch effects often
#' manifest as multiplicative scaling factors. By constructing ratios with
#' batch-associated denominators, the batch-specific scaling cancels out,
#' improving cross-cohort transferability.
#'
#' @param x Expression matrix (samples x genes). Should contain log-transformed
#'   expression values when `log_transform = TRUE`.
#' @param numerators Character vector of numerator gene names (signal genes).
#' @param denominators Character vector of denominator gene names (batch-proxy
#'   genes, typically from [select_denominator_features()]).
#' @param log_transform Logical. If `TRUE` (default), compute log-ratios by

#'   subtraction (appropriate for log-transformed expression). If `FALSE`,
#'   compute raw ratios by division.
#'
#' @return A matrix of ratio features with dimensions (n_samples x
#'   (n_numerators * n_denominators)). Column names are formatted as
#'   `"Numerator/Denominator"`.
#'
#' @details
#' The function creates all pairwise combinations of numerator and denominator
#' genes. For n numerators and d denominators, this produces n*d ratio features.
#'
#' Missing genes (not found in the input matrix) are silently filtered out.
#' The function will error if no valid numerators or denominators remain after
#' filtering.
#'
#' @examples
#' \dontrun{
#' # First select denominators
#' denom_result <- select_denominator_features(x, cohort, y)
#'
#' # Get signal genes (e.g., from DE analysis)
#' numerators <- select_de_features(
#'   x_list = list(x1, x2),
#'   y_list = list(y1, y2),
#'   n_features = 30
#' )
#'
#' # Construct ratio features
#' x_ratios <- construct_dual_stream_ratios(
#'   x = x,
#'   numerators = numerators,
#'   denominators = denom_result$denominators
#' )
#'
#' # Use ratios for optimization
#' result <- optimize_panel(x_ratios, y, ...)
#' }
#'
#' @seealso [select_denominator_features()] for selecting batch-proxy genes,
#'   [select_de_features()] for selecting signal genes,
#'   [pairwise_col_diff()] for an alternative all-pairs ratio construction.
#'
#' @export
construct_dual_stream_ratios <- function(x,
                                         numerators,
                                         denominators,
                                         log_transform = TRUE) {
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be a matrix or data.frame.", call. = FALSE)
  }
  x <- as.matrix(x)

  if (is.null(colnames(x))) {
    stop("`x` must have column names (gene identifiers).", call. = FALSE)
  }

  if (!is.character(numerators) || length(numerators) == 0L) {
    stop("`numerators` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(denominators) || length(denominators) == 0L) {
    stop("`denominators` must be a non-empty character vector.", call. = FALSE)
  }

  # Filter to genes present in data
  numerators_present <- intersect(numerators, colnames(x))
  denominators_present <- intersect(denominators, colnames(x))

  if (length(numerators_present) == 0L) {
    stop("No numerator genes found in expression matrix.", call. = FALSE)
  }
  if (length(denominators_present) == 0L) {
    stop("No denominator genes found in expression matrix.", call. = FALSE)
  }

  n_missing_num <- length(numerators) - length(numerators_present)
  n_missing_denom <- length(denominators) - length(denominators_present)
  if (n_missing_num > 0L || n_missing_denom > 0L) {
    warning(
      sprintf(
        paste0("Filtered genes not in matrix: %d numerators, %d denominators. ",
               "Ratios are built from the remaining genes only."),
        n_missing_num, n_missing_denom
      ),
      call. = FALSE
    )
  }

  # Create all pairwise ratio features
  n_ratios <- length(numerators_present) * length(denominators_present)
  ratio_mat <- matrix(NA_real_, nrow = nrow(x), ncol = n_ratios)
  ratio_names <- character(n_ratios)

  idx <- 1L
  for (num in numerators_present) {
    for (denom in denominators_present) {
      if (log_transform) {
        # Assuming data is already log-transformed, subtract directly
        ratio_mat[, idx] <- x[, num] - x[, denom]
      } else {
        # Raw ratio with small constant to avoid division by zero
        ratio_mat[, idx] <- x[, num] / (x[, denom] + 1e-10)
      }
      ratio_names[idx] <- paste0(num, "/", denom)
      idx <- idx + 1L
    }
  }

  colnames(ratio_mat) <- ratio_names
  rownames(ratio_mat) <- rownames(x)

  ratio_mat
}
