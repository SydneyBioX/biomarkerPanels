#' Redundancy-Aware Feature Selection
#'
#' Greedy minimum-redundancy maximum-relevance (mRMR) selection, for pools
#' where every panel slot costs an assay and should not be spent on
#' co-regulated features.
#'
#' @name nonredundant_features
NULL

#' Select Non-Redundant Features by Greedy mRMR
#'
#' Builds a ranked list one feature at a time. Each feature's relevance is its
#' moderated limma |t| for the case/control contrast on the pooled cohorts
#' (cohort as a covariate), scaled so the most relevant feature scores 1. At
#' each step the candidate with the highest relevance minus its redundancy
#' with the features already selected is added. Redundancy is the absolute
#' correlation between two features measured *within* each cohort and averaged
#' over cohorts, so shared batch structure is not mistaken for biological
#' co-regulation.
#'
#' The first feature is always the most relevant one, so the ranking diverges
#' from [select_de_features()] only from the second feature onwards.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Provide a list when
#'   `x_list` is a list.
#' @param n_features Number of features to return (default `50`).
#' @param penalty How a candidate's redundancy with the already-selected
#'   features is summarised: `"max"` (default) uses its largest |correlation|
#'   with any selected feature, so a near-duplicate of any member is always
#'   penalised in full; `"mean"` uses the average, which lets a near-duplicate
#'   of one member through once enough unrelated features have been selected.
#' @param cor_method Correlation used for redundancy: `"spearman"` (default)
#'   or `"pearson"`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to analysis.
#' @return Character vector of feature identifiers in selection order.
#'   Features whose relevance is undefined are dropped.
#' @details Only features shared by every cohort are considered. Each cohort
#'   needs at least three samples to estimate correlations. With missing
#'   values, Spearman ranks are taken per feature over its observed samples and
#'   correlations use pairwise-complete observations; an undefined correlation
#'   (e.g. a feature constant within a cohort) counts as zero redundancy.
#'
#'   Each step correlates every candidate with the newly added feature only,
#'   so no features x features matrix is formed and the cost grows linearly
#'   in `n_features`.
#' @seealso [select_de_features()], [select_ruleout_features()],
#'   [compute_diversity()]
#' @export
select_nonredundant_features <- function(x_list,
                                         y_list,
                                         n_features = 50L,
                                         penalty = c("max", "mean"),
                                         cor_method = c("spearman", "pearson"),
                                         assay = NULL) {
  n_features <- .validate_positive_integer(n_features, "n_features")
  penalty <- match.arg(penalty)
  cor_method <- match.arg(cor_method)

  prepared <- .prepare_selection_inputs(
    x_list, y_list,
    assay = assay,
    response = "factor"
  )
  matrices <- prepared$matrices
  cohort_names <- prepared$cohort_names

  too_small <- vapply(matrices, nrow, integer(1)) < 3L
  if (any(too_small)) {
    stop(
      "Cohort(s) ", paste0("`", cohort_names[too_small], "`", collapse = ", "),
      " need at least three samples to estimate redundancy.",
      call. = FALSE
    )
  }

  relevance <- .pooled_de_relevance(matrices, prepared$responses)
  relevance <- relevance[!is.na(relevance)]
  candidates <- names(relevance)

  # Ranking once up front makes each step a plain Pearson correlation.
  cor_inputs <- lapply(matrices, function(mat) {
    mat <- mat[, candidates, drop = FALSE]
    if (cor_method == "spearman") .rank_columns(mat) else mat
  })
  redundancy_with <- function(j) {
    abs_cor <- vapply(cor_inputs, .abs_cor_with_column, numeric(length(candidates)), j = j)
    if (is.null(dim(abs_cor))) {
      abs_cor <- matrix(abs_cor, ncol = length(cor_inputs))
    }
    redundancy <- rowMeans(abs_cor, na.rm = TRUE)
    redundancy[is.nan(redundancy)] <- 0
    redundancy
  }

  n_select <- min(n_features, length(candidates))
  selected <- integer(n_select)
  selected[1L] <- which.max(relevance)
  redundancy_sum <- redundancy_max <- redundancy_with(selected[1L])

  for (i in seq_len(n_select)[-1L]) {
    redundancy <- if (penalty == "mean") redundancy_sum / (i - 1L) else redundancy_max
    score <- relevance - redundancy
    score[selected[seq_len(i - 1L)]] <- -Inf
    selected[i] <- which.max(score)
    added <- redundancy_with(selected[i])
    redundancy_sum <- redundancy_sum + added
    redundancy_max <- pmax(redundancy_max, added)
  }

  candidates[selected]
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Pooled Differential Expression Relevance
#'
#' Fits `~ cohort + label` with limma on the row-bound cohorts (`~ label` for a
#' single cohort) and returns |moderated t| for the label coefficient, scaled
#' so the maximum is 1.
#'
#' @param matrices List of samples x features matrices with identical columns.
#' @param responses List of factors with levels `c("No", "Yes")`.
#' @return Named numeric vector of relevance scores in `[0, 1]` (`NA` where
#'   limma could not estimate the contrast).
#' @noRd
.pooled_de_relevance <- function(matrices, responses) {
  x_pool <- do.call(rbind, matrices)
  label <- factor(
    unlist(lapply(responses, as.character), use.names = FALSE),
    levels = c("No", "Yes")
  )
  design <- if (length(matrices) > 1L) {
    cohort <- factor(rep(seq_along(matrices), vapply(matrices, nrow, integer(1))))
    stats::model.matrix(~ cohort + label)
  } else {
    stats::model.matrix(~label)
  }

  fit <- limma::eBayes(limma::lmFit(t(x_pool), design))
  abs_t <- abs(fit$t[, "labelYes"])
  names(abs_t) <- colnames(x_pool)

  max_t <- suppressWarnings(max(abs_t, na.rm = TRUE))
  if (!is.finite(max_t) || max_t == 0) {
    stop(
      "No feature shows any case/control difference, so relevance is ",
      "undefined. Check that cases and controls are not confounded with cohort.",
      call. = FALSE
    )
  }
  abs_t / max_t
}

#' Column-Wise Ranks Keeping Missing Values
#'
#' @param mat Numeric matrix.
#' @return Matrix of the same shape holding within-column ranks, `NA` kept.
#' @noRd
.rank_columns <- function(mat) {
  matrix(
    apply(mat, 2, rank, na.last = "keep"),
    nrow = nrow(mat),
    dimnames = dimnames(mat)
  )
}

#' Absolute Correlation of Every Column With One Column
#'
#' @param mat Numeric matrix.
#' @param j Column index.
#' @return Numeric vector of |correlations|, `NA` where undefined.
#' @noRd
.abs_cor_with_column <- function(mat, j) {
  use <- if (anyNA(mat)) "pairwise.complete.obs" else "everything"
  # A constant column has no defined correlation; it is scored as zero
  # redundancy by the caller, so the warning carries no information.
  withCallingHandlers(
    abs(stats::cor(mat, mat[, j], use = use))[, 1L],
    warning = function(w) {
      if (grepl("standard deviation is zero", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
