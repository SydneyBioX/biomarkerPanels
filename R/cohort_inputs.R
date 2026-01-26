#' Cohort Input Preparation
#'
#' Internal functions for preparing and standardizing cohort data before
#' optimization. Handles multi-cohort inputs, feature extraction, and
#' response standardization.
#'
#' @name cohort_inputs
#' @keywords internal
NULL

#' Prepare Cohort Inputs for Optimization
#'
#' Main internal function that processes input data (single or multi-cohort)
#' and prepares it for NSGA-II optimization. Handles feature extraction,
#' alignment, transformation, and response standardization.
#'
#' @param x Matrix-like object, SummarizedExperiment, or list of such objects.
#' @param y Binary response vector or list of response vectors.
#' @param assay For SummarizedExperiment inputs, assay name or index.
#' @param transform Name of the feature transform to apply.
#' @param feature_subset Optional subset of features to keep.
#' @param feature_alignment Alignment strategy for multi-cohort data.
#' @return List with x (combined matrix), truth (factor), cohort (factor),
#'   cohort_names, and cohort_counts.
#' @keywords internal
.prepare_cohort_inputs <- function(x, y, assay = NULL,
                                   transform = "none",
                                   feature_subset = NULL,
                                   feature_alignment = "intersection") {

  # Handle list input (matrices, data.frames, or SummarizedExperiments)
  if (is.list(x) && !inherits(x, "SummarizedExperiment")) {
    if (!is.list(y)) {
      stop("When `x` is a list, `y` must also be a list.", call. = FALSE)
    }
    if (length(x) != length(y)) {
      stop("`x` and `y` lists must have the same length.", call. = FALSE)
    }

    cohort_names <- names(x)
    if (is.null(cohort_names) || any(cohort_names == "")) {
      cohort_names <- sprintf("cohort_%02d", seq_along(x))
    }

    matrices <- lapply(seq_along(x), function(i) {
      .extract_feature_matrix(x[[i]], assay = assay)
    })

    for (i in seq_along(matrices)) {
      if (is.null(colnames(matrices[[i]]))) {
        colnames(matrices[[i]]) <- sprintf("feature_%04d", seq_len(ncol(matrices[[i]])))
      }
    }

    feature_sets <- lapply(matrices, colnames)

    # Apply feature alignment strategy
    aligned <- .align_features(matrices, feature_sets, feature_alignment)
    matrices <- aligned$matrices
    ordered_features <- aligned$features

    if (!is.null(feature_subset)) {
      missing <- setdiff(feature_subset, ordered_features)
      if (length(missing)) {
        stop(
          "Feature(s) not found across all cohorts: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      matrices <- lapply(matrices, function(mat) {
        # Preserve reference_feature attribute through subsetting
        ref_attr <- attr(mat, "reference_feature")
        result <- mat[, feature_subset, drop = FALSE]
        if (!is.null(ref_attr)) {
          attr(result, "reference_feature") <- ref_attr
        }
        result
      })
    }

    matrices <- .apply_feature_transform(matrices, transform)

    feature_sets <- lapply(matrices, colnames)
    if (any(vapply(feature_sets, is.null, logical(1)))) {
      stop("Feature transform produced matrices without column names.", call. = FALSE)
    }
    common_features <- Reduce(intersect, feature_sets)
    if (is.null(common_features) || !length(common_features)) {
      # Provide detailed diagnostics to help users understand the issue
      feature_counts <- vapply(feature_sets, length, integer(1))
      stop(
        "No shared features were found across cohorts after aggregation.\n",
        "Feature counts per cohort after aggregation: ",
        paste(sprintf("cohort %d: %d", seq_along(feature_counts), feature_counts),
              collapse = ", "), ".\n",
        "This typically occurs when:\n",
        "  1. Cohorts have non-overlapping base features before aggregation\n",
        "  2. feature_pool specifies aggregated names (e.g., 'A--B') that don't exist ",
        "in all cohorts\n",
        "Solution: Specify base feature names (not aggregated names like 'A--B') in ",
        "feature_pool, or ensure cohorts share common features.",
        call. = FALSE
      )
    }
    ordered_features <- feature_sets[[1]][feature_sets[[1]] %in% common_features]
    matrices <- lapply(matrices, function(mat) {
      mat[, ordered_features, drop = FALSE]
    })

    counts <- vapply(matrices, nrow, integer(1))
    response_list <- lapply(seq_along(y), function(i) {
      yi <- ensure_binary_response(y[[i]])
      if (length(yi) != counts[[i]]) {
        stop("Length of `y[[", i, "]]` must match number of samples in `x[[", i, "]]`.",
             call. = FALSE)
      }
      yi
    })

    # Validate that all cohorts have consistent factor levels after standardization
    # ensure_binary_response() should produce c("No", "Yes") levels, but verify
    reference_levels <- levels(response_list[[1]])
    for (i in seq_along(response_list)[-1]) {
      current_levels <- levels(response_list[[i]])
      if (!identical(reference_levels, current_levels)) {
        stop(
          "Inconsistent response factor levels across cohorts. ",
          "Cohort 1 has levels: ", paste(reference_levels, collapse = ", "), "; ",
          "Cohort ", i, " has levels: ", paste(current_levels, collapse = ", "), ". ",
          "Use ensure_binary_response() with consistent positive/negative arguments.",
          call. = FALSE
        )
      }
    }

    combined_x <- do.call(rbind, matrices)
    # Safe concatenation: all responses are now guaranteed to have identical levels
    truth <- factor(
      unlist(lapply(response_list, as.character), use.names = FALSE),
      levels = reference_levels
    )
    cohort <- factor(rep(cohort_names, counts), levels = cohort_names)

    list(
      x = combined_x,
      truth = truth,
      cohort = cohort,
      cohort_names = cohort_names,
      cohort_counts = setNames(as.list(counts), cohort_names)
    )
  }

  else if (inherits(x, "SummarizedExperiment")) {
    # Handle single SummarizedExperiment input (same as matrix case)
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    if (is.null(colnames(x_mat))) {
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }

    if (!is.null(feature_subset)) {
      missing <- setdiff(feature_subset, colnames(x_mat))
      if (length(missing)) {
        stop(
          "Feature(s) not found in `x`: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      ref_attr <- attr(x_mat, "reference_feature")
      x_mat <- x_mat[, feature_subset, drop = FALSE]
      if (!is.null(ref_attr)) {
        attr(x_mat, "reference_feature") <- ref_attr
      }
    }

    x_mat <- .apply_feature_transform(list(x_mat), transform)[[1]]
    if (is.null(colnames(x_mat))) {
      stop("`x` must have column names in order to align with panel features.",
           call. = FALSE)
    }
    cohort_names <- "cohort_01"
    list(
      x = x_mat,
      truth = truth,
      cohort = factor(rep(cohort_names, nrow(x_mat)), levels = cohort_names),
      cohort_names = cohort_names,
      cohort_counts = setNames(list(nrow(x_mat)), cohort_names)
    )
  }

  # Handle single matrix/data.frame input
  else {
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    if (is.null(colnames(x_mat))) {
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }

    if (!is.null(feature_subset)) {
      missing <- setdiff(feature_subset, colnames(x_mat))
      if (length(missing)) {
        stop(
          "Feature(s) not found in `x`: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      # Preserve reference_feature attribute through subsetting
      ref_attr <- attr(x_mat, "reference_feature")
      x_mat <- x_mat[, feature_subset, drop = FALSE]
      if (!is.null(ref_attr)) {
        attr(x_mat, "reference_feature") <- ref_attr
      }
    }

    x_mat <- .apply_feature_transform(list(x_mat), transform)[[1]]
    if (is.null(colnames(x_mat))) {
      stop("`x` must have column names in order to align with panel features.",
           call. = FALSE)
    }
    cohort_names <- "cohort_01"
    list(
      x = x_mat,
      truth = truth,
      cohort = factor(rep(cohort_names, nrow(x_mat)), levels = cohort_names),
      cohort_names = cohort_names,
      cohort_counts = setNames(list(nrow(x_mat)), cohort_names)
    )
  }
}

#' Extract Feature Matrix from Various Input Types
#'
#' Converts input data (matrix, data.frame, or SummarizedExperiment) to a
#' numeric matrix.
#'
#' @param x Input object.
#' @param assay For SummarizedExperiment, assay name or index.
#' @return Numeric matrix.
#' @keywords internal
.extract_feature_matrix <- function(x, assay = NULL) {
  if (inherits(x, "SummarizedExperiment")) {
    assays <- SummarizedExperiment::assayNames(x)
    if (is.null(assay)) {
      assay <- assays[1]
    }
    return(SummarizedExperiment::assay(x, assay))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`x` must be a matrix-like object.", call. = FALSE)
  }
  mode(x) <- "numeric"
  x
}
