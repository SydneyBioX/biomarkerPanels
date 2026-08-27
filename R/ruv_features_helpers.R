#' RUV Feature Selection Helpers
#'
#' Internal helpers backing [select_ruv_features()]: control-gene selection,
#' RUV-4 fitting, factor-count selection, feature scoring, and diagnostics.
#' Split out of `ruv_features.R` to separate the public entry point from its
#' implementation detail.
#'
#' @name ruv_features_helpers
NULL


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Pool Cohort Data into Single Matrix
#'
#' Combines per-cohort matrices and responses from `.prepare_selection_inputs()`
#' into a single pooled dataset suitable for RUV-4.
#'
#' @param matrices List of per-cohort feature matrices (samples x features).
#' @param responses List of per-cohort 0/1 integer response vectors.
#' @param cohort_names Character vector of cohort identifiers.
#' @return List with `Y` (pooled matrix), `X_outcome` (design matrix),
#'   `cohort` (factor), `feature_names`.
#' @noRd
.pool_cohort_data <- function(matrices, responses, cohort_names) {
  Y <- do.call(rbind, matrices)
  outcome <- unlist(responses, use.names = FALSE)
  X_outcome <- cbind(outcome = outcome)
  cohort <- factor(
    rep(cohort_names, vapply(matrices, nrow, integer(1))),
    levels = cohort_names
  )
  list(
    Y = Y,
    X_outcome = X_outcome,
    cohort = cohort,
    feature_names = colnames(Y)
  )
}

#' Empirical Negative Controls from a Pooled Limma Fit
#'
#' Fits a naive pooled limma model (expression ~ outcome, ignoring cohort) and
#' flags genes whose p-value falls above the `neg_control_quantile` threshold
#' as negative controls. Shared by the initial control selection and the
#' iterative refinement so the two passes cannot drift apart.
#'
#' @param Y Pooled expression matrix (samples x genes).
#' @param X_outcome Design matrix for outcome (samples x 1).
#' @param neg_control_quantile Quantile threshold for empirical selection.
#' @return Logical vector of length `ncol(Y)`.
#' @noRd
.pooled_pvalue_controls <- function(Y, X_outcome, neg_control_quantile) {
  expr <- t(Y) # genes x samples for limma
  design <- stats::model.matrix(~X_outcome[, 1])
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  pvals <- fit$p.value[, 2]
  threshold <- stats::quantile(pvals, probs = neg_control_quantile, na.rm = TRUE)
  pvals > threshold
}

#' Select Negative Control Genes
#'
#' Identifies genes to serve as negative controls for RUV-4. Empirical mode
#' fits a naive pooled limma model (expression ~ outcome, ignoring cohort) and
#' selects non-significant genes.
#'
#' @param Y Pooled expression matrix (samples x genes).
#' @param X_outcome Design matrix for outcome (samples x 1).
#' @param feature_names Character vector of gene names.
#' @param method `"empirical"` or `"list"`.
#' @param neg_control_genes User-supplied gene names (for `method = "list"`).
#' @param neg_control_quantile Quantile threshold for empirical selection.
#' @return Logical vector of length `ncol(Y)`.
#' @noRd
.select_negative_controls <- function(Y, X_outcome, feature_names, method,
                                      neg_control_genes, neg_control_quantile) {
  if (method == "list") {
    missing <- setdiff(neg_control_genes, feature_names)
    if (length(missing)) {
      stop(
        "Negative control genes not found in feature names: ",
        paste(head(missing, 5), collapse = ", "),
        if (length(missing) > 5) paste0(" (and ", length(missing) - 5, " more)"),
        call. = FALSE
      )
    }
    return(feature_names %in% neg_control_genes)
  }

  # Empirical: pooled limma ignoring cohort
  ctl <- .pooled_pvalue_controls(Y, X_outcome, neg_control_quantile)

  if (sum(ctl) < 20L) {
    warning(
      "Only ", sum(ctl), " negative control genes selected. ",
      "Consider lowering `neg_control_quantile` for more robust RUV-4 estimation.",
      call. = FALSE
    )
  }
  ctl
}

#' Auto-Select Number of Unwanted Factors
#'
#' Determines k by PCA on the negative control gene submatrix, finding the
#' number of components needed to explain 90% of variance.
#'
#' @param Y_ctl Control gene submatrix (samples x control genes).
#' @param n_samples Total number of samples.
#' @return Integer k >= 1.
#' @noRd
.auto_select_k <- function(Y_ctl, n_samples) {
  if (ncol(Y_ctl) < 2L) {
    return(1L)
  }
  pca <- stats::prcomp(Y_ctl, center = TRUE, scale. = FALSE)
  var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  k <- which(var_explained >= 0.90)[1]
  if (is.na(k)) k <- length(pca$sdev)
  k <- min(k, 10L, floor(n_samples / 5))
  max(k, 1L)
}

#' Cap k to Prevent Singular Matrices
#' @param k Requested k.
#' @param n_ctl Number of control genes.
#' @param n_samples Number of samples.
#' @return Capped integer k.
#' @noRd
.cap_k <- function(k, n_ctl, n_samples) {
  max_k <- min(n_ctl - 1L, floor(n_samples / 5))
  if (k > max_k) {
    message("Reducing k from ", k, " to ", max_k,
      " to avoid singular matrix (n_controls=", n_ctl,
      ", n_samples=", n_samples, ")."
    )
    k <- max_k
  }
  max(as.integer(k), 1L)
}

#' Fit RUV-4 Model
#'
#' @param Y Pooled expression matrix (samples x genes).
#' @param X_outcome Design matrix (samples x 1).
#' @param ctl Logical vector of negative controls.
#' @param k Number of unwanted factors.
#' @param ebayes Whether to apply empirical Bayes variance adjustment.
#' @return RUV-4 fit object.
#' @noRd
.fit_ruv4 <- function(Y, X_outcome, ctl, k, ebayes) {
  fit <- ruv::RUV4(Y = Y, X = X_outcome, ctl = ctl, k = k)
  if (ebayes) {
    tryCatch(
      {
        fit <- ruv::variance_adjust(fit, ebayes = TRUE)
      },
      error = function(e) {
        warning(
          "ruv::variance_adjust() failed: ", conditionMessage(e),
          ". Using unadjusted t-statistics.",
          call. = FALSE
        )
      }
    )
  }
  fit
}

#' Iteratively Refine Negative Controls
#'
#' Corrects expression using current RUV-4 fit, re-runs pooled limma to update
#' the control set, and re-fits RUV-4. Repeats up to `max_iterations` times or
#' until the control set stabilises.
#'
#' @param Y Pooled expression matrix.
#' @param X_outcome Design matrix.
#' @param ctl Initial logical control vector.
#' @param k Number of unwanted factors.
#' @param ebayes Whether to apply EB adjustment.
#' @param neg_control_quantile Quantile threshold.
#' @param max_iterations Maximum iterations.
#' @return List with `fit` and `ctl`.
#' @noRd
.iterate_ruv_controls <- function(Y, X_outcome, ctl, k, ebayes,
                                  neg_control_quantile, max_iterations) {
  fit <- .fit_ruv4(Y, X_outcome, ctl, k, ebayes)

  for (iter in seq_len(max_iterations)) {
    Y_corrected <- Y - fit$W %*% fit$alpha
    new_ctl <- .pooled_pvalue_controls(Y_corrected, X_outcome,
      neg_control_quantile
    )

    if (identical(which(new_ctl), which(ctl))) break
    ctl <- new_ctl
    k_safe <- .cap_k(k, sum(ctl), nrow(Y))
    fit <- .fit_ruv4(Y, X_outcome, ctl, k_safe, ebayes)
  }

  list(fit = fit, ctl = ctl)
}

#' Score Features by RUV-4 Dual Objective
#'
#' @param fit RUV-4 fit object.
#' @param scoring `"ratio"` or `"composite_z"`.
#' @param epsilon Small constant for ratio denominator.
#' @param feature_names Character vector of gene names.
#' @return Data frame sorted by decreasing score.
#' @noRd
.score_ruv_features <- function(fit, scoring, epsilon, feature_names) {
  t_stat <- fit$t
  if (is.matrix(t_stat)) {
    discrimination <- abs(t_stat[1, ])
  } else {
    discrimination <- abs(t_stat)
  }
  names(discrimination) <- feature_names

  alpha <- fit$alpha
  if (is.matrix(alpha) && nrow(alpha) > 0L) {
    alpha_norm <- sqrt(colSums(alpha^2))
  } else {
    alpha_norm <- rep(0, length(feature_names))
  }
  names(alpha_norm) <- feature_names

  if (scoring == "ratio") {
    score <- discrimination / (alpha_norm + epsilon)
  } else {
    z_disc <- as.numeric(scale(discrimination))
    z_alpha <- as.numeric(scale(-alpha_norm))
    z_disc[is.na(z_disc)] <- 0
    z_alpha[is.na(z_alpha)] <- 0
    score <- z_disc + z_alpha
  }

  df <- data.frame(
    feature = feature_names,
    discrimination = discrimination,
    alpha_norm = alpha_norm,
    score = score,
    stringsAsFactors = FALSE
  )

  df <- df[order(df$score, decreasing = TRUE), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  rownames(df) <- NULL
  df
}

#' Compute RUV-4 Diagnostics
#'
#' @param fit RUV-4 fit object.
#' @param cohort Factor of cohort labels.
#' @param feature_names Character vector of gene names, ordered as the columns
#'   of the pooled expression matrix passed to RUV-4.
#' @return List of diagnostic values.
#' @noRd
.ruv_diagnostics <- function(fit, cohort, feature_names) {
  k <- ncol(fit$W)
  W <- fit$W
  outcome <- fit$X[, 1]

  W_outcome_cor <- vapply(seq_len(k), function(j) {
    stats::cor(W[, j], outcome)
  }, numeric(1))

  if (nlevels(cohort) >= 2L) {
    cohort_int <- as.integer(cohort)
    W_cohort_cor <- vapply(seq_len(k), function(j) {
      stats::cor(W[, j], cohort_int)
    }, numeric(1))
  } else {
    W_cohort_cor <- rep(NA_real_, k)
  }

  # Warn if unwanted factors capture biological signal
  if (any(abs(W_outcome_cor) > 0.3)) {
    flagged <- which(abs(W_outcome_cor) > 0.3)
    warning(
      "Unwanted factor(s) ", paste(flagged, collapse = ", "),
      " correlate with outcome (|cor| > 0.3). ",
      "Consider reducing k to avoid removing biological signal.",
      call. = FALSE
    )
  }

  # ruv::RUV4 does not return the expression matrix, so control gene names are
  # recovered from `feature_names` (columns of the pooled matrix) via fit$ctl.
  ctl_idx <- if (is.logical(fit$ctl)) which(fit$ctl) else fit$ctl
  control_genes <- if (length(feature_names) && length(ctl_idx)) {
    feature_names[ctl_idx]
  } else {
    character()
  }

  list(
    k = k,
    n_controls = length(ctl_idx),
    control_genes = control_genes,
    W_outcome_cor = W_outcome_cor,
    W_cohort_cor = W_cohort_cor
  )
}
