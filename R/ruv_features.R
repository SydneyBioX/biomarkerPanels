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
#'   [get_top_de_features()]
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
    if (!is.numeric(k) || length(k) != 1L || k < 1L) {
      stop("`k` must be a positive integer or NULL.", call. = FALSE)
    }
    k <- as.integer(k)
  }
  if (!is.numeric(neg_control_quantile) || length(neg_control_quantile) != 1L ||
    neg_control_quantile < 0 || neg_control_quantile > 1) {
    stop("`neg_control_quantile` must be a numeric scalar in [0, 1].",
      call. = FALSE
    )
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive numeric scalar.", call. = FALSE)
  }
  if (neg_control_method == "list" && is.null(neg_control_genes)) {
    stop('`neg_control_genes` must be provided when neg_control_method = "list".',
      call. = FALSE
    )
  }

  prepared <- .prepare_ridge_inputs(x_list, y_list, assay = assay)
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
  diagnostics <- .ruv_diagnostics(fit, pooled$cohort)

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


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Check ruv Package Availability
#' @keywords internal
.check_ruv_available <- function() {
  if (!requireNamespace("ruv", quietly = TRUE)) {
    stop(
      "Package 'ruv' is required for RUV-4 feature selection but not installed.\n",
      "Install it via: install.packages('ruv')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Pool Cohort Data into Single Matrix
#'
#' Combines per-cohort matrices and responses from [.prepare_ridge_inputs()]
#' into a single pooled dataset suitable for RUV-4.
#'
#' @param matrices List of per-cohort feature matrices (samples x features).
#' @param responses List of per-cohort 0/1 integer response vectors.
#' @param cohort_names Character vector of cohort identifiers.
#' @return List with `Y` (pooled matrix), `X_outcome` (design matrix),
#'   `cohort` (factor), `feature_names`.
#' @keywords internal
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
#' @keywords internal
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
  expr <- t(Y) # genes x samples for limma
  design <- stats::model.matrix(~X_outcome[, 1])
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  pvals <- fit$p.value[, 2]

  threshold <- stats::quantile(pvals, probs = neg_control_quantile, na.rm = TRUE)
  ctl <- pvals > threshold

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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
.iterate_ruv_controls <- function(Y, X_outcome, ctl, k, ebayes,
                                  neg_control_quantile, max_iterations) {
  fit <- .fit_ruv4(Y, X_outcome, ctl, k, ebayes)

  for (iter in seq_len(max_iterations)) {
    Y_corrected <- Y - fit$W %*% fit$alpha
    expr_corr <- t(Y_corrected)
    design <- stats::model.matrix(~X_outcome[, 1])
    lm_fit <- limma::lmFit(expr_corr, design)
    lm_fit <- limma::eBayes(lm_fit)
    pvals <- lm_fit$p.value[, 2]

    threshold <- stats::quantile(pvals, probs = neg_control_quantile,
      na.rm = TRUE
    )
    new_ctl <- pvals > threshold

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
#' @keywords internal
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
#' @return List of diagnostic values.
#' @keywords internal
.ruv_diagnostics <- function(fit, cohort) {
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

  ctl_genes <- colnames(fit$Y)
  if (is.null(ctl_genes)) ctl_genes <- character()
  ctl_idx <- if (is.logical(fit$ctl)) which(fit$ctl) else fit$ctl
  control_genes <- if (length(ctl_genes) && length(ctl_idx)) {
    ctl_genes[ctl_idx]
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
