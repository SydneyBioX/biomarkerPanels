#' Cross-Platform Omics Prediction (CPOP) Feature Selection
#'
#' Faithful implementation of the CPOP algorithm (Wang et al., 2022). CPOP
#' selects pairwise log-ratio features that are *jointly* informative across
#' independent cohorts, and prunes those whose effects disagree in direction.
#' This sits alongside [select_transferable_features()] — that function is a
#' single-pass ridge surrogate, whereas [select_cpop_features()] is the full
#' iterated CPOP procedure.
#'
#' @name cpop_features
NULL

#' Select Features via the CPOP Algorithm
#'
#' Runs CPOP step 1 (iterated elastic-net selection with cohort intersection)
#' followed by CPOP step 2 (sign-agreement pruning using ridge refits) on the
#' pairwise log-ratio representation of the input. Features are returned as
#' `"A--B"` pair names produced by [pairwise_col_diff()]; for log-normalised
#' expression data these are mathematically equivalent to log-ratios.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows are samples, columns are features. Original CPOP is defined
#'   for two cohorts; this implementation generalises by intersecting selections
#'   across all cohorts at every CPOP step.
#' @param y_list A binary response factor aligned with `x_list` (levels
#'   `c("No", "Yes")`). Provide a list when `x_list` is a list.
#' @param n_features Maximum number of pair features to retain (default `50`).
#' @param n_iter Maximum number of iterations for both steps (default `20`).
#' @param alpha Elastic-net mixing parameter for step 1. Default `1` (lasso).
#'   A numeric vector triggers CPOP's alpha sweep: features carried over from
#'   earlier alphas are exempted from the penalty in later passes.
#' @param lambda_choice Either `"lambda.min"` (default) or `"lambda.1se"`.
#' @param penalty Optional numeric vector of `glmnet` `penalty.factor` values
#'   named by pair feature. When `NULL`, [colmeans_penalty()] is used, matching
#'   the CPOP default.
#' @param sign_consistency_threshold Minimum fraction of cohorts required to
#'   agree on coefficient sign in step 2. Default `1.0` is the canonical CPOP
#'   rule (unanimous agreement).
#' @param prune_collinear Logical; if `TRUE` (default), removes redundant
#'   ratios by retaining a minimum spanning tree of the implied feature graph
#'   (matching CPOP's `mst_lratio`).
#' @param family `glmnet` family. Default `"binomial"`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index.
#' @param cv_control Optional named list of arguments forwarded to
#'   [glmnet::cv.glmnet()].
#' @return A list with the selected pair features, the union of base genes
#'   they reference, the step-1 intermediate set, per-iteration logs, and the
#'   ridge coefficient matrix from the final step-2 pass.
#' @seealso [fit_cpop_panel()], [select_transferable_features()]
#' @export
select_cpop_features <- function(x_list, y_list,
                                 n_features = 50L,
                                 n_iter = 20L,
                                 alpha = 1,
                                 lambda_choice = c("lambda.min", "lambda.1se"),
                                 penalty = NULL,
                                 sign_consistency_threshold = 1.0,
                                 prune_collinear = TRUE,
                                 family = "binomial",
                                 assay = NULL,
                                 cv_control = list()) {
  lambda_choice <- match.arg(lambda_choice)
  n_features <- .validate_positive_integer(n_features, "n_features")
  n_iter <- .validate_positive_integer(n_iter, "n_iter")
  if (!is.numeric(alpha) || any(!is.finite(alpha)) || any(alpha < 0) || any(alpha > 1)) {
    stop("`alpha` must be numeric in [0, 1].", call. = FALSE)
  }
  .validate_probability(sign_consistency_threshold, "sign_consistency_threshold",
                        bounds = "closed")
  alpha <- sort(unique(alpha), decreasing = TRUE)

  prepared <- .cpop_prepare_inputs(x_list, y_list, assay = assay)
  pair_matrices <- prepared$pair_matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  pair_features <- prepared$pair_features

  if (is.null(penalty)) {
    penalty_vec <- colmeans_penalty(pair_matrices)
    names(penalty_vec) <- pair_features
  } else {
    if (!is.numeric(penalty) || length(penalty) != length(pair_features) ||
        is.null(names(penalty)) ||
        !setequal(names(penalty), pair_features)) {
      stop("`penalty` must be a numeric vector named with pair feature names.",
           call. = FALSE)
    }
    penalty_vec <- penalty[pair_features]
  }

  cpop1_features <- character(0)
  step1_logs <- list()
  for (a in alpha) {
    w <- penalty_vec
    w[cpop1_features] <- 0
    s1 <- .cpop_step1(pair_matrices, responses, w, n_iter, a, n_features,
                      lambda_choice, family, cv_control)
    step1_logs[[paste0("alpha_", a)]] <- s1$log
    cpop1_features <- unique(c(cpop1_features, s1$features))
    if (length(cpop1_features) >= n_features) {
      cpop1_features <- cpop1_features[seq_len(n_features)]
      break
    }
  }

  empty_coef <- matrix(numeric(0), nrow = 0L, ncol = length(pair_matrices),
                       dimnames = list(NULL, cohort_names))
  settings <- list(
    alpha = alpha, n_iter = n_iter, n_features = n_features,
    lambda_choice = lambda_choice,
    sign_consistency_threshold = sign_consistency_threshold,
    prune_collinear = prune_collinear, family = family
  )

  if (!length(cpop1_features)) {
    return(list(features = character(0), base_features = character(0),
                cpop1_features = character(0), step1 = step1_logs,
                step2_coefficients = empty_coef, settings = settings))
  }

  s2 <- .cpop_step2(pair_matrices, responses, cpop1_features, n_iter,
                    lambda_choice, sign_consistency_threshold, family,
                    cv_control)
  cpop2_features <- s2$features
  if (!length(cpop2_features)) {
    warning("CPOP step 2 removed all features; using step 1 features.",
            call. = FALSE)
    cpop2_features <- cpop1_features
  }

  final_features <- if (prune_collinear && length(cpop2_features) >= 2L) {
    .mst_lratio_pruning(cpop2_features)
  } else {
    cpop2_features
  }

  coef_out <- if (!is.null(s2$coefficients) && nrow(s2$coefficients) > 0L) {
    colnames(s2$coefficients) <- cohort_names
    s2$coefficients
  } else {
    empty_coef
  }

  list(
    features = final_features,
    base_features = sort(unique(unlist(
      strsplit(final_features, "--", fixed = TRUE)))),
    cpop1_features = cpop1_features,
    step1 = step1_logs,
    step2_coefficients = coef_out,
    settings = settings
  )
}

#' CPOP Column-Mean Difference Penalty
#'
#' Computes the CPOP default `penalty.factor` for [glmnet::cv.glmnet()]: a
#' p-transformed vector of cross-cohort absolute mean differences. Pair
#' features with similar means across cohorts receive smaller penalties (and
#' are therefore preferred), encouraging selection of platform-stable ratios.
#'
#' Generalises to N cohorts by using the per-feature mean range (max minus
#' min across cohort means). For two cohorts this reduces to the absolute
#' difference used in the original CPOP implementation.
#'
#' @param x_list A list of pair-feature matrices (one per cohort), or a single
#'   matrix (returns an unit-weight vector).
#' @return A numeric vector of penalty factors of length `ncol(x_list[[1]])`,
#'   normalised so the values sum to the vector length.
#' @export
colmeans_penalty <- function(x_list) {
  x_list <- .as_cohort_list(x_list)
  if (length(x_list) < 2L) {
    return(rep(1, ncol(x_list[[1]])))
  }
  ncols <- vapply(x_list, ncol, integer(1))
  if (length(unique(ncols)) != 1L) {
    stop("All cohort matrices must have the same number of columns.",
         call. = FALSE)
  }
  means <- do.call(rbind, lapply(x_list, colMeans, na.rm = TRUE))
  raw <- if (nrow(means) == 2L) {
    abs(means[1L, ] - means[2L, ])
  } else {
    apply(means, 2L, function(v) max(v) - min(v))
  }
  .p_transform(raw)
}

.cpop_step1 <- function(pair_matrices, responses, w, n_iter, alpha, n_features,
                        lambda_choice, family, cv_control) {
  pair_features <- colnames(pair_matrices[[1]])
  p <- length(pair_features)
  remaining <- pair_features
  selected <- character(0)
  step_log <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    if (length(selected) >= n_features || !length(remaining)) break

    per_cohort <- vector("list", length(pair_matrices))
    for (k in seq_along(pair_matrices)) {
      cv_args <- utils::modifyList(list(
        x = pair_matrices[[k]][, remaining, drop = FALSE],
        y = if (family == "binomial") {
          as.integer(responses[[k]]) - 1L
        } else responses[[k]],
        alpha = alpha,
        family = family,
        penalty.factor = w[remaining]
      ), cv_control)
      cv_fit <- do.call(glmnet::cv.glmnet, cv_args)
      coef_mat <- as.matrix(stats::coef(cv_fit, s = lambda_choice))
      nz <- rownames(coef_mat)[coef_mat[, 1] != 0]
      per_cohort[[k]] <- setdiff(nz, "(Intercept)")
    }
    added <- Reduce(intersect, per_cohort)
    step_log[[i]] <- list(iteration = i, per_cohort = per_cohort,
                          added = added)
    selected <- unique(c(selected, added))
    remaining <- setdiff(pair_features, selected)
  }
  list(features = selected, log = step_log[!vapply(step_log, is.null, logical(1))])
}
