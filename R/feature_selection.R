#' Get Top Differentially Expressed Features
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
#' @export
#'
#' @importFrom stats model.matrix pnorm qnorm weighted.mean sd
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
get_top_de_features <- function(x_list,
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

  if (!length(limma_stats$feature_names)) {
    return(character())
  }

  ordered_p <- .aggregate_de_pvalues(limma_stats$t_matrix, combination_method)
  if (!length(ordered_p)) {
    return(character())
  }

  head(names(ordered_p), n = min(n_features, length(ordered_p)))
}

#' Select Features for Ratio Construction
#'
#' Identifies both highly informative (differentially expressed) features and
#' features with stable expression across cohorts, supporting ratio-based panel
#' design. Stability can be measured via precision-weighted inverse t-statistics,
#' coefficient-of-variation penalties, or inverse t-statistic standard errors.
#'
#' @inheritParams get_top_de_features
#' @param n_stable Number of approximately invariant features to return
#'   (default `50`).
#' @param n_informative Number of informative features to return (default `50`).
#' @param stability_method Method used to score stability. One of
#'   `"precision_weighted"`, `"cv_t_stats"`, or `"inverse_t_se"`.
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

#' Select Ridge-Stable Transferable Features
#'
#' Fits ridge-penalised logistic regression models for each cohort and selects
#' features whose coefficients are both large in magnitude and consistent across
#' cohorts. This emphasises biomarkers that deliver transferable signal.
#'
#' @inheritParams get_top_de_features
#' @param n_features Maximum number of features to return (default `50`).
#' @param lambda Optional penalty parameter supplied to [glmnet::glmnet()]. When
#'   `NULL`, cohort-specific values are chosen via cross-validation.
#' @param lambda_choice When `lambda` is `NULL`, choose either `"lambda_1se"`
#'   (default) or `"lambda_min"` from [glmnet::cv.glmnet()] as the operating
#'   point for coefficient extraction.
#' @param min_coefficient Minimum absolute coefficient required (applied to the
#'   cohort-wise minimum magnitude) before a feature is considered transferable.
#' @param require_sign_consistency If `TRUE`, only retain features whose ridge
#'   coefficients share the same sign (ignoring near-zero coefficients) across
#'   cohorts.
#' @param standardize Logical; passed to [glmnet::glmnet()] to control feature
#'   standardisation prior to fitting. Defaults to `TRUE`.
#' @param cv_control Optional named list of additional arguments forwarded to
#'   [glmnet::cv.glmnet()] when `lambda` is `NULL`.
#' @return A list containing the selected feature identifiers, per-feature
#'   scoring metadata, coefficient matrix, and the lambda value used for each
#'   cohort.
#' @details Features are currently aligned across cohorts via a simple
#'   intersection of shared column names. Future versions will add more flexible
#'   harmonisation strategies (e.g., imputation or reference-based mapping). 
#'   I promise I'll get to this lol #TODO
#' @export
#'
#' @importFrom glmnet cv.glmnet glmnet
select_transferable_features <- function(x_list,
                                         y_list,
                                         n_features = 50L,
                                         lambda = NULL,
                                         lambda_choice = c("lambda_1se", "lambda_min"),
                                         min_coefficient = 0,
                                         require_sign_consistency = TRUE,
                                         standardize = TRUE,
                                         assay = NULL,
                                         cv_control = list()) {
  n_features <- .validate_positive_integer(n_features, "n_features")
  lambda_choice <- match.arg(lambda_choice)
  if (!is.numeric(min_coefficient) || length(min_coefficient) != 1L) {
    stop("`min_coefficient` must be a numeric scalar.", call. = FALSE)
  }

  prepared <- .prepare_ridge_inputs(x_list, y_list, assay = assay)
  matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  feature_names <- prepared$feature_names

  if (!length(feature_names)) {
    return(list(
      features = character(),
      scores = data.frame(),
      coefficients = matrix(numeric(), nrow = 0, ncol = length(cohort_names),
                            dimnames = list(character(), cohort_names)),
      lambda = numeric(),
      settings = list(
        lambda_choice = lambda_choice,
        min_coefficient = min_coefficient,
        require_sign_consistency = require_sign_consistency,
        standardize = standardize
      )
    ))
  }

  lambda_vec <- .resolve_lambda_vector(lambda, length(matrices))

  ridge_results <- lapply(seq_along(matrices), function(i) {
    x_mat <- matrices[[i]]
    y_vec <- responses[[i]]

    if (is.null(lambda_vec)) {
      cv_args <- utils::modifyList(list(
        x = x_mat,
        y = y_vec,
        alpha = 0,
        family = "binomial",
        standardize = standardize
      ), cv_control)

      cv_fit <- do.call(glmnet::cv.glmnet, cv_args)
      chosen_lambda <- if (lambda_choice == "lambda_min") {
        cv_fit$lambda.min
      } else {
        cv_fit$lambda.1se
      }

      coef_vec <- .extract_glmnet_coefficients(cv_fit, chosen_lambda, feature_names)
      list(coefficients = coef_vec, lambda = chosen_lambda)
    } else {
      fit <- glmnet::glmnet(
        x = x_mat,
        y = y_vec,
        alpha = 0,
        family = "binomial",
        standardize = standardize,
        lambda = lambda_vec[[i]]
      )
      coef_vec <- .extract_glmnet_coefficients(fit, lambda_vec[[i]], feature_names)
      list(coefficients = coef_vec, lambda = lambda_vec[[i]])
    }
  })

  coefficient_matrix <- do.call(cbind, lapply(ridge_results, `[[`, "coefficients"))
  colnames(coefficient_matrix) <- cohort_names

  lambda_used <- vapply(ridge_results, `[[`, numeric(1), "lambda")
  names(lambda_used) <- cohort_names

  combined_scores <- .score_transferable_features(
    coefficient_matrix,
    min_coefficient = min_coefficient,
    require_sign_consistency = require_sign_consistency
  )

  if (!nrow(combined_scores)) {
    return(list(
      features = character(),
      scores = combined_scores,
      coefficients = coefficient_matrix[NULL, , drop = FALSE],
      lambda = lambda_used,
      settings = list(
        lambda_choice = if (is.null(lambda_vec)) lambda_choice else "fixed",
        min_coefficient = min_coefficient,
        require_sign_consistency = require_sign_consistency,
        standardize = standardize
      )
    ))
  }

  ranked <- order(combined_scores$score, decreasing = TRUE)
  top_indices <- head(ranked, n = min(n_features, length(ranked)))
  selected_features <- rownames(combined_scores)[top_indices]

  list(
    features = selected_features,
    scores = combined_scores[selected_features, , drop = FALSE],
    coefficients = coefficient_matrix[selected_features, , drop = FALSE],
    lambda = lambda_used,
    settings = list(
      lambda_choice = if (is.null(lambda_vec)) lambda_choice else "fixed",
      min_coefficient = min_coefficient,
      require_sign_consistency = require_sign_consistency,
      standardize = standardize
    )
  )
}

.compute_limma_statistics <- function(x_list,
                                      y_list,
                                      contrast = NULL,
                                      assay = NULL) {
  x_list <- .as_cohort_list(x_list)
  y_list <- .as_cohort_list(y_list)

  if (length(x_list) != length(y_list)) {
    stop("`x_list` and `y_list` must have the same length.", call. = FALSE)
  }

  cohort_names <- names(x_list)
  if (is.null(cohort_names) || any(!nzchar(cohort_names))) {
    cohort_names <- sprintf("cohort_%02d", seq_along(x_list))
  }

  t_list <- vector("list", length(x_list))
  se_list <- vector("list", length(x_list))

  for (i in seq_along(x_list)) {
    x_mat <- .extract_feature_matrix(x_list[[i]], assay = assay)
    y_vec <- ensure_binary_response(y_list[[i]])

    if (nrow(x_mat) != length(y_vec)) {
      stop("Number of rows in `x_list[[", i, "]]` must match the length of ",
           "`y_list[[", i, "]]`.", call. = FALSE)
    }

    expr <- t(x_mat)
    design <- stats::model.matrix(~0 + y_vec)
    colnames(design) <- levels(y_vec)

    contrast_str <- .resolve_contrast(
      contrast = contrast,
      cohort_index = i,
      total_cohorts = length(x_list),
      level_names = levels(y_vec)
    )

    fit <- limma::lmFit(expr, design = design)
    cm <- limma::makeContrasts(contrasts = contrast_str, levels = design)
    fit2 <- limma::contrasts.fit(fit, cm)
    efit <- limma::eBayes(fit2, robust = TRUE)

    tt <- limma::topTable(
      efit,
      coef = 1,
      number = Inf,
      sort.by = "none"
    )

    if (!"t" %in% colnames(tt)) {
      stop("`limma::topTable()` output is missing the 't' column.", call. = FALSE)
    }

    t_vals <- tt[["t"]]
    names(t_vals) <- rownames(tt)

    se_vals <- sqrt(efit$s2.post) * efit$stdev.unscaled[, 1]
    se_vals <- se_vals[rownames(tt)]
    names(se_vals) <- rownames(tt)

    t_list[[i]] <- t_vals
    se_list[[i]] <- se_vals
  }

  feature_names <- unique(unlist(lapply(t_list, names), use.names = FALSE))

  if (!length(feature_names)) {
    empty <- matrix(numeric(0), nrow = 0, ncol = length(x_list))
    return(list(
      t_matrix = empty,
      se_matrix = empty,
      feature_names = character(),
      cohort_names = cohort_names
    ))
  }

  t_matrix <- matrix(
    NA_real_,
    nrow = length(feature_names),
    ncol = length(t_list),
    dimnames = list(feature_names, cohort_names)
  )
  se_matrix <- matrix(
    NA_real_,
    nrow = length(feature_names),
    ncol = length(se_list),
    dimnames = list(feature_names, cohort_names)
  )

  for (i in seq_along(t_list)) {
    idx <- match(names(t_list[[i]]), feature_names)
    t_matrix[idx, i] <- t_list[[i]]
    se_matrix[idx, i] <- se_list[[i]]
  }

  list(
    t_matrix = t_matrix,
    se_matrix = se_matrix,
    feature_names = feature_names,
    cohort_names = cohort_names
  )
}

.aggregate_de_pvalues <- function(t_matrix, combination_method) {
  if (length(t_matrix) == 0L) {
    return(numeric())
  }

  z_scores <- apply(t_matrix, 2, .rank_to_z)
  if (is.vector(z_scores)) {
    z_scores <- matrix(
      z_scores,
      ncol = 1,
      dimnames = list(rownames(t_matrix), colnames(t_matrix))
    )
  }

  combined_p <- apply(z_scores, 1, function(z_vec) {
    valid <- !is.na(z_vec) & is.finite(z_vec)
    if (!any(valid)) {
      return(NA_real_)
    }
    .combine_pvalues(z_vec[valid], combination_method)
  })

  if (is.null(names(combined_p)) && !is.null(rownames(z_scores))) {
    names(combined_p) <- rownames(z_scores)
  }

  combined_p <- pmax(combined_p, .Machine$double.xmin, na.rm = FALSE)
  combined_z <- stats::qnorm(combined_p, lower.tail = FALSE)
  two_sided <- 2 * stats::pnorm(-abs(combined_z))
  two_sided <- two_sided[!is.na(two_sided)]

  sort(two_sided, decreasing = FALSE)
}

.select_stable_genes <- function(t_matrix, se_matrix, method, top_n) {
  if (length(t_matrix) == 0L || nrow(t_matrix) == 0L) {
    return(character())
  }

  abs_t <- abs(t_matrix)

  scores <- vapply(seq_len(nrow(abs_t)), function(idx) {
    t_row <- abs_t[idx, ]
    se_row <- se_matrix[idx, ]

    if (method == "precision_weighted") {
      precision_row <- 1 / (se_row + 1e-8)
      valid <- !is.na(t_row) & !is.na(precision_row) & is.finite(precision_row) & precision_row > 0
      if (!any(valid)) {
        return(NA_real_)
      }
      stats::weighted.mean(1 / (t_row[valid] + 0.01), precision_row[valid])
    } else if (method == "cv_t_stats") {
      valid <- !is.na(t_row) & is.finite(t_row)
      if (!any(valid)) {
        return(NA_real_)
      }
      vals <- t_row[valid]
      mean_abs <- mean(vals)
      sd_abs <- stats::sd(vals)
      if (is.na(sd_abs)) {
        sd_abs <- 0
      }
      1 / ((sd_abs + 1e-6) * (mean_abs + 0.01))
    } else if (method == "inverse_t_se") {
      valid <- !is.na(t_row) & !is.na(se_row) & is.finite(t_row) & is.finite(se_row)
      if (!any(valid)) {
        return(NA_real_)
      }
      mean((t_row[valid] + 0.01) * (se_row[valid] + 1e-8))
    } else {
      NA_real_
    }
  }, numeric(1))

  if (method == "inverse_t_se") {
    scores <- 1 / (scores + 1e-6)
  }

  scores[!is.finite(scores)] <- NA_real_
  valid_idx <- which(!is.na(scores))
  if (!length(valid_idx)) {
    return(character())
  }

  ordered <- valid_idx[order(scores[valid_idx], decreasing = TRUE)]
  top_k <- head(ordered, n = min(top_n, length(ordered)))
  rownames(abs_t)[top_k]
}

.rank_to_z <- function(values) { # nocov start
  if (all(is.na(values))) {
    out <- rep(NA_real_, length(values))
    names(out) <- names(values)
    return(out)
  }
  valid <- !is.na(values) & is.finite(values)
  out <- rep(NA_real_, length(values))
  if (!any(valid)) {
    names(out) <- names(values)
    return(out)
  }
  ranks <- rank(values[valid], ties.method = "average")
  probs <- ranks / (length(ranks) + 1)
  out[valid] <- stats::qnorm(probs)
  names(out) <- names(values)
  out
} # nocov end

.combine_pvalues <- function(z_scores, method) { # nocov start
  if (!length(z_scores)) {
    return(NA_real_)
  }
  if (method == "Stouffer") {
    return(stats::pnorm(sum(z_scores), sd = sqrt(length(z_scores)), lower.tail = FALSE))
  }
  if (method == "OSP") {
    p_vals <- stats::pnorm(z_scores, lower.tail = TRUE)
    return(stats::pchisq(-2 * sum(log(p_vals)), df = 2 * length(p_vals), lower.tail = TRUE))
  }
  if (method == "Fisher") {
    p_vals <- stats::pnorm(z_scores, lower.tail = FALSE)
    return(stats::pchisq(-2 * sum(log(p_vals)), df = 2 * length(p_vals), lower.tail = FALSE))
  }
  if (method == "maxP") {
    return(stats::pnorm(max(z_scores), lower.tail = FALSE))
  }
  stop("Unsupported combination method: ", method, call. = FALSE)
} # nocov end

.resolve_contrast <- function(contrast, cohort_index, total_cohorts, level_names) {
  if (length(level_names) < 2L) {
    stop("`y_list` must contain exactly two levels.", call. = FALSE)
  }
  if (is.null(contrast)) {
    return(paste(level_names[2], level_names[1], sep = "-"))
  }
  if (length(contrast) == 1L) {
    return(contrast)
  }
  if (length(contrast) == total_cohorts) {
    return(contrast[[cohort_index]])
  }
  stop("`contrast` must be length 1 or match the number of cohorts.",
       call. = FALSE)
}

.as_cohort_list <- function(x) { # nocov start
  if (is.null(x)) {
    stop("Input cannot be NULL.", call. = FALSE)
  }
  if (is.list(x)) {
    return(x)
  }
  list(x)
} # nocov end

.validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1L) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.prepare_ridge_inputs <- function(x_list, y_list, assay = NULL) {
  x_list <- .as_cohort_list(x_list)
  y_list <- .as_cohort_list(y_list)

  if (length(x_list) != length(y_list)) {
    stop("`x_list` and `y_list` must have the same length.", call. = FALSE)
  }

  cohort_names <- names(x_list)
  if (is.null(cohort_names) || any(!nzchar(cohort_names))) {
    cohort_names <- sprintf("cohort_%02d", seq_along(x_list))
  }

  matrices <- vector("list", length(x_list))
  responses <- vector("list", length(x_list))

  for (i in seq_along(x_list)) {
    x_mat <- .extract_feature_matrix(x_list[[i]], assay = assay)
    y_vec <- ensure_binary_response(y_list[[i]])

    if (nrow(x_mat) != length(y_vec)) {
      stop("Number of rows in `x_list[[", i, "]]` must match the length of ",
           "`y_list[[", i, "]]`.", call. = FALSE)
    }

    if (is.null(colnames(x_mat))) {
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }

    matrices[[i]] <- x_mat
    responses[[i]] <- as.integer(y_vec) - 1L
  }

  feature_sets <- lapply(matrices, colnames)
  # TODO(#transferability): replace simple intersection alignment with robust harmonisation.
  common_features <- Reduce(intersect, feature_sets)
  if (is.null(common_features) || !length(common_features)) {
    stop(
      "No shared features were found across cohorts. Provide data with ",
      "overlapping feature identifiers (future releases will support more ",
      "flexible alignment).",
      call. = FALSE
    )
  }
  ordered_features <- feature_sets[[1]][feature_sets[[1]] %in% common_features]
  matrices <- lapply(matrices, function(mat) {
    mat[, ordered_features, drop = FALSE]
  })

  list(
    matrices = matrices,
    responses = responses,
    cohort_names = cohort_names,
    feature_names = ordered_features
  )
}

.resolve_lambda_vector <- function(lambda, num_cohorts) {
  if (is.null(lambda)) {
    return(NULL)
  }
  if (!is.numeric(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("`lambda` must be positive numeric.", call. = FALSE)
  }
  if (length(lambda) == 1L) {
    rep(lambda, num_cohorts)
  } else if (length(lambda) == num_cohorts) {
    lambda
  } else {
    stop("`lambda` must be length 1 or match the number of cohorts.", call. = FALSE)
  }
}

.extract_glmnet_coefficients <- function(model, lambda, feature_names) {
  coef_mat <- stats::coef(model, s = lambda)
  coef_vec <- as.numeric(coef_mat)[-1]
  names(coef_vec) <- rownames(coef_mat)[-1]
  coef_vec[feature_names]
}

.score_transferable_features <- function(coefficient_matrix,
                                         min_coefficient,
                                         require_sign_consistency) {
  if (!nrow(coefficient_matrix)) {
    return(data.frame())
  }

  abs_coeff <- abs(coefficient_matrix)
  mean_abs <- rowMeans(abs_coeff)
  sd_coeff <- apply(coefficient_matrix, 1, stats::sd)
  min_abs <- apply(abs_coeff, 1, min)

  sign_consistent <- apply(coefficient_matrix, 1, function(vals) {
    nz <- vals[abs(vals) > 1e-6]
    if (!length(nz)) {
      return(TRUE)
    }
    all(nz >= 0) || all(nz <= 0)
  })

  score <- mean_abs / (sd_coeff + 1e-6)

  keep <- is.finite(score) & !is.na(score) & (min_abs >= min_coefficient)
  if (require_sign_consistency) {
    keep <- keep & sign_consistent
  }

  score_df <- data.frame(
    mean_abs = mean_abs,
    sd = sd_coeff,
    min_abs = min_abs,
    sign_consistent = sign_consistent,
    score = score,
    stringsAsFactors = FALSE
  )

  score_df <- score_df[keep, , drop = FALSE]
  score_df[order(score_df$score, decreasing = TRUE), , drop = FALSE]
}
