#' Discriminative Feature Selection with Cohort Penalty
#'
#' Functions for selecting features based on within-cohort discriminative AUC
#' penalised by cross-cohort separability and AUC variability.
#'
#' @name discriminative_features
NULL

#' Select Discriminative Features with Cohort Penalty
#'
#' Ranks features by a composite score that rewards within-cohort discrimination
#' (outcome AUC) while penalising features that separate cohorts (batch
#' confounding) or have inconsistent performance across cohorts (AUC variance).
#'
#' The composite score for feature \eqn{j} is:
#' \deqn{score_j = \bar{AUC}_j - \lambda_1 \cdot cohort\_disc_j - \lambda_2 \cdot sd(AUC_j)}
#'
#' where \eqn{\bar{AUC}_j} is the mean within-cohort AUC, \eqn{cohort\_disc_j}
#' measures how well feature \eqn{j} separates cohorts (irrespective of
#' outcome), and \eqn{sd(AUC_j)} captures cross-cohort variability.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Must be a factor with
#'   levels `c("No", "Yes")`. Provide a list when `x_list` is a list.
#' @param n_features Maximum number of features to return (default `50`).
#' @param lambda_cohort Penalty weight for cohort discrimination (default `1.0`).
#'   Higher values more aggressively penalise features that separate cohorts.
#' @param lambda_sd Penalty weight for cross-cohort AUC variability (default
#'   `1.0`). Higher values favour features with consistent AUC across cohorts.
#' @param min_auc Minimum mean within-cohort AUC required (default `0.6`).
#'   Features below this threshold are excluded before ranking.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to analysis.
#' @return A list with components:
#' \describe{
#'   \item{features}{Character vector of selected feature names.}
#'   \item{scores}{Data frame with per-feature scoring metadata: `feature`,
#'     `mean_auc`, `sd_auc`, `min_auc`, `cohort_disc`, `score`.}
#'   \item{per_cohort_auc}{Matrix of within-cohort AUC values (features x
#'     cohorts).}
#'   \item{settings}{List of parameter values used.}
#' }
#' @export
#'
#' @importFrom stats sd
select_discriminative_features <- function(x_list,
                                           y_list,
                                           n_features = 50L,
                                           lambda_cohort = 1.0,
                                           lambda_sd = 1.0,
                                           min_auc = 0.6,
                                           assay = NULL) {
  n_features <- .validate_positive_integer(n_features, "n_features")
  .validate_numeric_scalar(lambda_cohort, "lambda_cohort")
  .validate_numeric_scalar(lambda_sd, "lambda_sd")
  .validate_probability(min_auc, "min_auc", bounds = "closed")

  prepared <- .prepare_selection_inputs(x_list, y_list, assay = assay)
  matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  feature_names <- prepared$feature_names
  n_cohorts <- length(matrices)

  # Step 1: Within-cohort discrimination AUC (features x cohorts)
  auc_matrix <- vapply(seq_along(matrices), function(k) {
    .compute_rank_auc(matrices[[k]], responses[[k]])
  }, numeric(length(feature_names)))

  if (is.null(dim(auc_matrix))) {
    auc_matrix <- matrix(auc_matrix, ncol = n_cohorts)
  }
  rownames(auc_matrix) <- feature_names
  colnames(auc_matrix) <- cohort_names

  # Step 2: Cohort discrimination (only meaningful with 2+ cohorts)
  if (n_cohorts >= 2L) {
    x_pooled <- do.call(rbind, matrices)
    cohort_labels <- rep(seq_len(n_cohorts), vapply(matrices, nrow, integer(1)))
    cohort_disc <- .compute_cohort_discrimination(x_pooled, cohort_labels)
  } else {
    cohort_disc <- rep(0, length(feature_names))
  }
  names(cohort_disc) <- feature_names

  # Step 3 & 4: Score, filter, rank
  scored <- .score_discriminative_features(
    auc_matrix, cohort_disc,
    lambda_cohort = lambda_cohort,
    lambda_sd = lambda_sd,
    min_auc = min_auc
  )

  if (!nrow(scored)) {
    return(list(
      features = character(),
      scores = scored,
      per_cohort_auc = auc_matrix[integer(0), , drop = FALSE],
      settings = list(
        lambda_cohort = lambda_cohort,
        lambda_sd = lambda_sd,
        min_auc = min_auc
      )
    ))
  }

  top_n <- min(n_features, nrow(scored))
  selected <- scored[seq_len(top_n), , drop = FALSE]

  list(
    features = selected$feature,
    scores = selected,
    per_cohort_auc = auc_matrix[selected$feature, , drop = FALSE],
    settings = list(
      lambda_cohort = lambda_cohort,
      lambda_sd = lambda_sd,
      min_auc = min_auc
    )
  )
}


#' Vectorised Rank-Based AUC for All Columns
#'
#' Computes Wilcoxon-Mann-Whitney AUC for each column of a matrix against a
#' binary label. Uses the rank-sum formula for speed:
#' \eqn{AUC = (R_1 - n_1(n_1+1)/2) / (n_1 \cdot n_0)} where \eqn{R_1} is the
#' sum of ranks for the positive class.
#'
#' AUC is folded to \code{[0.5, 1]} so direction doesn't matter.
#'
#' @param x Numeric matrix (samples x features).
#' @param y Integer vector of 0/1 labels.
#' @return Named numeric vector of AUC values (length = ncol(x)).
#' @keywords internal
.compute_rank_auc <- function(x, y) {
  pos <- which(y == 1L)
  neg <- which(y == 0L)
  n1 <- length(pos)
  n0 <- length(neg)

  if (n1 < 1L || n0 < 1L) {
    return(rep(NA_real_, ncol(x)))
  }

  ranks <- apply(x, 2, rank, ties.method = "average")
  r1 <- colSums(ranks[pos, , drop = FALSE])
  auc <- (r1 - n1 * (n1 + 1) / 2) / (n1 * n0)

  # Fold to [0.5, 1] — direction-agnostic
  pmax(auc, 1 - auc)
}


#' Compute Cohort Discrimination per Feature
#'
#' For each feature, measures how well it separates cohorts (ignoring outcome).
#' Uses pairwise AUC between all cohort pairs, averaged into a single score.
#' For >2 cohorts, uses one-vs-rest macro AUC.
#'
#' @param x_pooled Numeric matrix of pooled samples (samples x features).
#' @param cohort_labels Integer vector of cohort assignments.
#' @return Numeric vector of cohort discrimination scores per feature, scaled
#'   to \code{[0.5, 1]} where 0.5 = no cohort separation.
#' @keywords internal
.compute_cohort_discrimination <- function(x_pooled, cohort_labels) {
  cohorts <- sort(unique(cohort_labels))
  n_cohorts <- length(cohorts)
  p <- ncol(x_pooled)

  if (n_cohorts < 2L) {
    return(rep(0.5, p))
  }

  if (n_cohorts == 2L) {
    idx1 <- which(cohort_labels == cohorts[1])
    idx2 <- which(cohort_labels == cohorts[2])
    binary <- integer(nrow(x_pooled))
    binary[idx2] <- 1L
    return(.compute_rank_auc(x_pooled, binary))
  }

  # >2 cohorts: one-vs-rest macro AUC
  auc_sum <- numeric(p)
  for (k in cohorts) {
    binary <- as.integer(cohort_labels == k)
    auc_sum <- auc_sum + .compute_rank_auc(x_pooled, binary)
  }
  auc_sum / n_cohorts
}


#' Score and Rank Features by Discriminative-Penalised Composite
#'
#' @param auc_matrix Features x cohorts matrix of within-cohort AUC values.
#' @param cohort_disc Numeric vector of cohort discrimination scores.
#' @param lambda_cohort Penalty weight for cohort discrimination.
#' @param lambda_sd Penalty weight for AUC variability.
#' @param min_auc Minimum mean AUC threshold.
#' @return Data frame sorted by decreasing score, filtered to `mean_auc > min_auc`.
#' @keywords internal
.score_discriminative_features <- function(auc_matrix,
                                           cohort_disc,
                                           lambda_cohort,
                                           lambda_sd,
                                           min_auc) {
  feature_names <- rownames(auc_matrix)
  n_cohorts <- ncol(auc_matrix)

  mean_auc <- rowMeans(auc_matrix)

  if (n_cohorts > 1L) {
    sd_auc <- apply(auc_matrix, 1, stats::sd)
  } else {
    sd_auc <- rep(0, length(mean_auc))
  }

  min_auc_vec <- apply(auc_matrix, 1, min)
  score <- mean_auc - lambda_cohort * (cohort_disc - 0.5) - lambda_sd * sd_auc

  df <- data.frame(
    feature = feature_names,
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    min_auc = min_auc_vec,
    cohort_disc = cohort_disc,
    score = score,
    stringsAsFactors = FALSE
  )

  df <- df[!is.na(df$score) & df$mean_auc >= min_auc, , drop = FALSE]
  df <- df[order(df$score, decreasing = TRUE), , drop = FALSE]
  rownames(df) <- NULL
  df
}
