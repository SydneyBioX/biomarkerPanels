#' Differential Expression Analysis
#'
#' Internal helpers for computing differential expression statistics via limma
#' and aggregating p-values across cohorts using meta-analysis methods. The
#' public entry point is [select_de_features()].
#'
#' @name de_features_helpers
NULL

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
  # Union alignment: a feature present in only some cohorts is still modelled
  # where it exists and padded with NA elsewhere, so partially-overlapping
  # cohorts still contribute evidence to the meta-analysis.
  prepared <- .prepare_selection_inputs(
    x_list, y_list,
    assay = assay,
    response = "factor",
    align = "union"
  )
  matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  feature_names <- prepared$feature_names

  t_list <- vector("list", length(matrices))
  se_list <- vector("list", length(matrices))

  for (i in seq_along(matrices)) {
    x_mat <- matrices[[i]]
    y_vec <- responses[[i]]

    expr <- t(x_mat)
    design <- stats::model.matrix(~ 0 + y_vec)
    colnames(design) <- levels(y_vec)

    contrast_str <- .resolve_contrast(
      contrast = contrast,
      cohort_index = i,
      total_cohorts = length(matrices),
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
