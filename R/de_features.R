#' Differential Expression Feature Selection
#'
#' Differential-expression-based selection of candidate features to seed
#' downstream optimization and ratio-construction workflows.
#'
#' @name de_features
NULL

#' Select Top Differentially Expressed Features
#'
#' Applies moderated t-statistics via `limma` across one or more cohorts, then
#' aggregates the resulting evidence into a ranked list of candidate features.
#' The top entries can be used to seed downstream optimization or ratio
#' construction workflows.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features. When a list
#'   is supplied each element is treated as an independent cohort.
#' @param y_list A binary response aligned with `x_list`. Must be a factor with
#'   levels `c("No", "Yes")`. Provide a list when `x_list` is a list.
#' @param contrast Optional contrast string (or list of strings) passed to
#'   [limma::makeContrasts()]. Defaults to `"Yes-No"` to align with
#'   [ensure_binary_response()].
#' @param n_features Number of features to return after ranking (default `50`).
#' @param combination_method Method for combining cohort-specific z-scores.
#'   One of `"OSP"`, `"Stouffer"`, `"Fisher"`, or `"maxP"`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to modelling.
#' @return Character vector of feature identifiers ordered by significance.
#' @seealso [select_transferable_features()], [select_discriminative_features()],
#'   [select_ruv_features()], [select_cpop_features()]
#' @export
#'
#' @importFrom stats model.matrix pnorm qnorm weighted.mean sd
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
select_de_features <- function(x_list,
                               y_list,
                               contrast = NULL,
                               n_features = 50L,
                               combination_method = c("OSP", "Stouffer", "Fisher", "maxP"),
                               assay = NULL) {
  combination_method <- match.arg(combination_method)
  n_features <- .validate_positive_integer(n_features, "n_features")

  limma_stats <- .compute_limma_statistics(
    x_list = x_list,
    y_list = y_list,
    contrast = contrast,
    assay = assay
  )

  ordered_p <- .aggregate_de_pvalues(limma_stats$t_matrix, combination_method)
  if (!length(ordered_p)) {
    stop(
      "Differential expression analysis produced no ranked features. ",
      "This may indicate that all features have identical expression across groups.",
      call. = FALSE
    )
  }

  head(names(ordered_p), n = min(n_features, length(ordered_p)))
}
