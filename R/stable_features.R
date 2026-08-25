#' Feature Stability Selection
#'
#' Functions for identifying stable features suitable for ratio construction,
#' using metrics based on t-statistic variability across cohorts.
#'
#' @name feature_stability
NULL

#' Select Features for Ratio Construction
#'
#' Identifies both highly informative (differentially expressed) features and
#' features with stable expression across cohorts, supporting ratio-based panel
#' design. Stability can be measured via precision-weighted inverse t-statistics,
#' coefficient-of-variation penalties, or inverse t-statistic standard errors.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Must be a factor with
#'   levels `c("No", "Yes")`. Provide a list when `x_list` is a list.
#' @param contrast Optional contrast string (or list of strings) passed to
#'   [limma::makeContrasts()]. Defaults to `"Yes-No"`.
#' @param n_stable Number of approximately invariant features to return
#'   (default `50`).
#' @param n_informative Number of informative features to return (default `50`).
#' @param stability_method Method used to score stability. One of
#'   `"precision_weighted"`, `"cv_t_stats"`, or `"inverse_t_se"`.
#' @param combination_method Method for combining cohort-specific z-scores.
#'   One of `"OSP"`, `"Stouffer"`, `"Fisher"`, or `"maxP"`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to modelling.
#' @return A list with elements `stable` and `informative`, each containing a
#'   character vector of feature identifiers.
#' @export
select_features_for_ratios <- function(x_list,
                                       y_list,
                                       contrast = NULL,
                                       n_stable = 50L,
                                       n_informative = 50L,
                                       stability_method = c("precision_weighted", "cv_t_stats", "inverse_t_se"),
                                       combination_method = c("OSP", "Stouffer", "Fisher", "maxP"),
                                       assay = NULL) {
  stability_method <- match.arg(stability_method)
  combination_method <- match.arg(combination_method)
  n_stable <- .validate_positive_integer(n_stable, "n_stable")
  n_informative <- .validate_positive_integer(n_informative, "n_informative")

  limma_stats <- .compute_limma_statistics(
    x_list = x_list,
    y_list = y_list,
    contrast = contrast,
    assay = assay
  )

  stable_genes <- .select_stable_genes(
    t_matrix = limma_stats$t_matrix,
    se_matrix = limma_stats$se_matrix,
    method = stability_method,
    top_n = n_stable
  )

  informative_order <- .aggregate_de_pvalues(limma_stats$t_matrix, combination_method)
  informative_genes <- head(
    names(informative_order),
    n = min(n_informative, length(informative_order))
  )

  list(
    stable = stable_genes,
    informative = informative_genes
  )
}

#' Select Stable Genes Based on T-statistic Variability
#'
#' Scores genes by stability across cohorts using the specified method and
#' returns the top-scoring features. Uses C++ for performance.
#'
#' @param t_matrix A matrix of t-statistics (features x cohorts).
#' @param se_matrix A matrix of standard errors (features x cohorts).
#' @param method One of `"precision_weighted"`, `"cv_t_stats"`, or
#'   `"inverse_t_se"`.
#' @param top_n Number of top stable features to return.
#' @return A character vector of feature names.
#' @noRd
.select_stable_genes <- function(t_matrix, se_matrix, method, top_n) {
  if (length(t_matrix) == 0L || nrow(t_matrix) == 0L) {
    return(character())
  }

  abs_t <- abs(t_matrix)

  # Uses C++ for speeeeedddd
  scores <- .select_stable_genes_cpp(abs_t, se_matrix, method)

  # Post-processing same as R version
  valid_idx <- which(!is.na(scores))
  if (!length(valid_idx)) {
    return(character())
  }

  ordered <- valid_idx[order(scores[valid_idx], decreasing = TRUE)]
  top_k <- head(ordered, n = min(top_n, length(ordered)))
  rownames(abs_t)[top_k]
}
