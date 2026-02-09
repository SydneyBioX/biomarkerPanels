#' Differential Expression Analysis
#'
#' Functions for computing differential expression statistics via limma and
#' aggregating p-values across cohorts using meta-analysis methods.
#'
#' @name de_analysis
NULL

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
    stop(
      "No features available for differential expression analysis. ",
      "Check that input matrices have column names and contain numeric data.",
      call. = FALSE
    )
  }

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

#' Compute Limma Statistics Across Cohorts
#'
#' Runs differential expression analysis via limma on each cohort and returns
#' the t-statistics and standard errors for all features.
#'
#' @param x_list A matrix, `SummarizedExperiment`, or list of such objects.
#' @param y_list A binary response aligned with `x_list`.
#' @param contrast Optional contrast string or list of strings.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index.
#' @return A list with `t_matrix`, `se_matrix`, `feature_names`, and
#'   `cohort_names`.
#' @keywords internal
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
        "`y_list[[", i, "]]`.",
        call. = FALSE
      )
    }

    # Ensure feature names exist (columns of x_mat become rows of expr)
    if (is.null(colnames(x_mat))) {
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }

    expr <- t(x_mat)
    design <- stats::model.matrix(~ 0 + y_vec)
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

  # Throw an error if no features are detected.
  if (!length(feature_names)) {
    stop("No features detected in input data after computing test statistics.", call. = FALSE)
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

#' Aggregate Differential Expression P-values Across Cohorts
#'
#' Converts t-statistics to z-scores and combines them across cohorts using the
#' specified meta-analysis method. Uses C++ for performance.
#'
#' @param t_matrix A matrix of t-statistics (features x cohorts).
#' @param combination_method One of `"OSP"`, `"Stouffer"`, `"Fisher"`, or
#'   `"maxP"`.
#' @return A named numeric vector of two-sided p-values, sorted ascending.
#' @keywords internal
.aggregate_de_pvalues <- function(t_matrix, combination_method) {
  if (length(t_matrix) == 0L) {
    return(numeric())
  }

  # Map method names to C++ expected format
  method_map <- c(
    "stouffer" = "Stouffer",
    "Stouffer" = "Stouffer",
    "fisher" = "Fisher",
    "Fisher" = "Fisher",
    "osp" = "OSP",
    "OSP" = "OSP",
    "maxp" = "maxP",
    "maxP" = "maxP"
  )

  cpp_method <- method_map[combination_method]
  if (is.na(cpp_method)) {
    stop("Unsupported combination method: ", combination_method, call. = FALSE)
  }

  .aggregate_de_pvalues_cpp(t_matrix, cpp_method)
}

#' Resolve Contrast String for Limma
#'
#' Returns the appropriate contrast string for a given cohort, handling both
#' single contrasts and cohort-specific contrast lists.
#'
#' @param contrast User-supplied contrast or `NULL` for default.
#' @param cohort_index Index of the current cohort.
#' @param total_cohorts Total number of cohorts.
#' @param level_names Factor levels from the response.
#' @return A contrast string for [limma::makeContrasts()].
#' @keywords internal
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
    call. = FALSE
  )
}

#' Convert Values to Z-scores via Rank Transformation
#'
#' Transforms values to z-scores using rank-based normal quantile
#' transformation. Used internally for p-value combination methods.
#'
#' @param values A numeric vector to transform.
#' @return A numeric vector of z-scores with the same names as input.
#' @keywords internal
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

#' Combine P-values using aggregation
#'
#' Combines z-scores into a single p-value using the specified method.
#'
#' @param z_scores A numeric vector of z-scores.
#' @param method One of `"Stouffer"`, `"OSP"`, `"Fisher"`, or `"maxP"`.
#' @return A single combined p-value.
#' @keywords internal
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

# ==============================================================================
# PURE R REFERENCE IMPLEMENTATIONS (for testing)
# ==============================================================================

#' Pure R Reference Implementation of DE P-value Aggregation
#'
#' Kept for regression testing against the C++ implementation.
#' TODO: Remove once Rcpp equivalents are fully debugged and validated.
#'
#' @inheritParams .aggregate_de_pvalues
#' @keywords internal
.aggregate_de_pvalues_pure_r <- function(t_matrix, combination_method) {
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
