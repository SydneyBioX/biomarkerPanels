#' RUV-4 Dual-Objective Feature Selection
#'
#' Functions for selecting features based on RUV-4 (Remove Unwanted Variation)
#' dual-objective scoring: rewarding discrimination strength while penalising
#' sensitivity to unwanted variation factors.
#'
#' @name ruv_features
NULL

#' Select Features via RUV-4 Dual-Objective Scoring
#'
#' Fits RUV-4 to jointly estimate outcome-associated signal and unwanted
#' variation per feature, then ranks features by the ratio (or composite
#' z-score) of discrimination strength to unwanted-variation sensitivity.
#' Unlike cohort-label-based methods, RUV-4 estimates latent batch factors
#' from negative control genes, making it suitable when cohort structure is
#' ambiguous or confounders are continuous.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Must be a factor with
#'   levels `c("No", "Yes")`. Provide a list when `x_list` is a list.
#' @param n_features Maximum number of features to return (default `50`).
#' @param k Number of unwanted variation factors to estimate. When `NULL`
#'   (default), auto-selected via PCA on the negative control submatrix.
#' @param neg_control_method Method for identifying negative control genes.
#'   `"empirical"` (default) selects non-significant genes from a pooled limma
#'   fit; `"list"` uses genes supplied via `neg_control_genes`.
#' @param neg_control_genes Character vector of negative control gene names.
#'   Required when `neg_control_method = "list"`.
#' @param neg_control_quantile Quantile threshold for empirical control
#'   selection (default `0.75`). Features with p-value above this quantile are
#'   selected as controls.
#' @param iterate_controls If `TRUE`, iteratively refine the negative control
#'   set by re-running limma on RUV-corrected expression. Default `FALSE`.
#' @param max_iterations Maximum number of control refinement iterations
#'   (default `3`).
#' @param scoring Scoring method. `"ratio"` (default) uses
#'   `discrimination / (alpha_norm + epsilon)`; `"composite_z"` uses
#'   `z(discrimination) + z(-alpha_norm)`.
#' @param epsilon Small constant added to alpha_norm denominator in ratio
#'   scoring (default `0.01`).
#' @param ebayes If `TRUE` (default), apply empirical Bayes variance shrinkage
#'   via [ruv::variance_adjust()] for better-calibrated t-statistics.
#' @param return_corrected If `TRUE`, include the RUV-4-corrected expression
#'   matrix in the output. Default `FALSE`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to analysis.
#'
#' @return A list with components:
#' \describe{
#'   \item{features}{Character vector of selected feature names.}
#'   \item{scores}{Data frame with per-feature scoring metadata: `feature`,
#'     `discrimination`, `alpha_norm`, `score`, `rank`.}
#'   \item{ruv_summary}{List of RUV-4 diagnostics: `k`, `n_controls`,
#'     `control_genes`, `W_outcome_cor`, `W_cohort_cor`.}
#'   \item{corrected}{RUV-4-corrected expression matrix (only when
#'     `return_corrected = TRUE`).}
#'   \item{settings}{List of parameter values used.}
#' }
#'
#' @details
#' This method complements the existing cohort-label-based selectors
#' ([select_transferable_features()], [select_discriminative_features()]) by
#' estimating latent unwanted factors from negative control genes rather than
#' requiring explicit cohort labels. It is particularly useful when batch
#' structure is ambiguous, confounders are continuous, or multiple sources of
#' unwanted variation overlap.
#'
#' The discrimination axis measures how strongly each feature associates with
#' the outcome after accounting for unwanted variation (via the RUV-4
#' t-statistic). The alpha_norm axis measures each feature's total loading on
#' the estimated unwanted factors (L2 norm across k factors). Features with
#' high discrimination and low alpha_norm are prioritised.
#'
#' @seealso [select_transferable_features()], [select_discriminative_features()],
#'   [select_de_features()], [select_batch_associated_features()] for the
#'   label-based counterpart that targets batch-associated genes directly.
#'
#' @export
select_ruv_features <- function(x_list,
                                y_list,
                                n_features = 50L,
                                k = NULL,
                                neg_control_method = c("empirical", "list"),
                                neg_control_genes = NULL,
                                neg_control_quantile = 0.75,
                                iterate_controls = FALSE,
                                max_iterations = 3L,
                                scoring = c("ratio", "composite_z"),
                                epsilon = 0.01,
                                ebayes = TRUE,
                                return_corrected = FALSE,
                                assay = NULL) {
  .check_ruv_available()
  n_features <- .validate_positive_integer(n_features, "n_features")
  neg_control_method <- match.arg(neg_control_method)
  scoring <- match.arg(scoring)

  if (!is.null(k)) {
    k <- .validate_positive_integer(k, "k")
  }
  .validate_probability(neg_control_quantile, "neg_control_quantile",
                        bounds = "closed")
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive numeric scalar.", call. = FALSE)
  }
  if (neg_control_method == "list" && is.null(neg_control_genes)) {
    stop('`neg_control_genes` must be provided when neg_control_method = "list".',
      call. = FALSE
    )
  }

  prepared <- .prepare_selection_inputs(x_list, y_list, assay = assay)
  pooled <- .pool_cohort_data(prepared$matrices, prepared$responses,
    prepared$cohort_names
  )

  ctl <- .select_negative_controls(
    Y = pooled$Y,
    X_outcome = pooled$X_outcome,
    feature_names = pooled$feature_names,
    method = neg_control_method,
    neg_control_genes = neg_control_genes,
    neg_control_quantile = neg_control_quantile
  )

  if (is.null(k)) {
    k <- .auto_select_k(pooled$Y[, ctl, drop = FALSE], nrow(pooled$Y))
  }
  k <- .cap_k(k, sum(ctl), nrow(pooled$Y))

  fit <- .fit_ruv4(pooled$Y, pooled$X_outcome, ctl, k, ebayes)

  if (iterate_controls) {
    result <- .iterate_ruv_controls(
      pooled$Y, pooled$X_outcome, ctl, k, ebayes,
      neg_control_quantile, max_iterations
    )
    fit <- result$fit
    ctl <- result$ctl
  }

  scored <- .score_ruv_features(fit, scoring, epsilon, pooled$feature_names)
  diagnostics <- .ruv_diagnostics(fit, pooled$cohort, pooled$feature_names)

  top_n <- min(n_features, nrow(scored))
  if (top_n == 0L) {
    selected <- scored[integer(0), , drop = FALSE]
  } else {
    selected <- scored[seq_len(top_n), , drop = FALSE]
  }

  corrected_out <- NULL
  if (return_corrected) {
    corrected_out <- pooled$Y - fit$W %*% fit$alpha
    colnames(corrected_out) <- pooled$feature_names
  }

  list(
    features = selected$feature,
    scores = selected,
    ruv_summary = diagnostics,
    corrected = corrected_out,
    settings = list(
      k = k,
      neg_control_method = neg_control_method,
      neg_control_quantile = neg_control_quantile,
      iterate_controls = iterate_controls,
      scoring = scoring,
      epsilon = epsilon,
      ebayes = ebayes
    )
  )
}
