#' Cohort Input Preparation
#'
#' Internal functions for preparing and standardizing cohort data before
#' optimization. Handles multi-cohort inputs, feature extraction, and
#' response standardization.
#'
#' @name cohort_inputs
#' @keywords internal
NULL

#' Default Cohort Names
#'
#' Returns the names of a cohort list, substituting `cohort_01`, `cohort_02`,
#' ... wherever a name is missing or empty, so every cohort has a stable
#' identifier for per-cohort bookkeeping and messages.
#'
#' @param x A list (or vector) of cohorts.
#' @return A character vector of names, one per element of `x`.
#' @keywords internal
.default_cohort_names <- function(x) {
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- sprintf("cohort_%02d", seq_along(x))
  }
  nm
}

#' Ensure a Feature Matrix Has Column Names
#'
#' Assigns generated `feature_0001`, `feature_0002`, ... column names when a
#' matrix has none, so downstream feature bookkeeping always has stable
#' identifiers. Matrices that already have column names are returned unchanged.
#'
#' @param mat A matrix-like object.
#' @return `mat` with non-`NULL` column names.
#' @keywords internal
.ensure_feature_colnames <- function(mat) {
  if (is.null(colnames(mat))) {
    colnames(mat) <- sprintf("feature_%04d", seq_len(ncol(mat)))
  }
  mat
}

#' Prepare Cohort Inputs for Optimization
#'
#' Main internal function that processes input data (single or multi-cohort)
#' and prepares it for NSGA optimization. Handles feature extraction,
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

    cohort_names <- .default_cohort_names(x)

    matrices <- lapply(seq_along(x), function(i) {
      .extract_feature_matrix(x[[i]], assay = assay)
    })

    matrices <- lapply(matrices, .ensure_feature_colnames)

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
    ref_attr <- attr(matrices[[1]], "reference_feature")
    if (!is.null(ref_attr)) {
      attr(combined_x, "reference_feature") <- ref_attr
    }

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

  # Handle single-cohort input; .extract_feature_matrix() dispatches on
  # matrix, data.frame, and SummarizedExperiment alike.
  else {
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    x_mat <- .ensure_feature_colnames(x_mat)

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
    return(t(SummarizedExperiment::assay(x, assay)))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`x` must be a matrix-like object.", call. = FALSE)
  }
  na_before <- sum(is.na(x))
  mode(x) <- "numeric"
  na_after <- sum(is.na(x))
  if (na_after > na_before) {
    warning(
      "Coercing matrix to numeric introduced ",
      na_after - na_before, " NA value(s). ",
      "Check that input data contains only numeric values.",
      call. = FALSE
    )
  }
  x
}

#' Report an Empty Cross-Cohort Feature Intersection
#'
#' Shared error for every path that aligns cohorts by intersecting feature
#' names. Reports the per-cohort feature counts so the caller can tell an
#' identifier-convention mismatch (all cohorts populated, zero overlap) from an
#' empty input.
#'
#' @param feature_sets List of per-cohort feature name vectors.
#' @param cohort_names Character vector of cohort identifiers.
#' @return Never returns; always signals an error.
#' @keywords internal
.stop_no_shared_features <- function(feature_sets, cohort_names) {
  counts <- vapply(feature_sets, length, integer(1))
  stop(
    "No shared features were found across cohorts.\n",
    "Feature counts per cohort: ",
    paste(sprintf("%s: %d", cohort_names, counts), collapse = ", "), ".\n",
    "Cohorts must share at least one feature identifier. Check that every ",
    "cohort uses the same naming convention (e.g. all gene symbols, or all ",
    "Ensembl IDs).",
    call. = FALSE
  )
}

#' Prepare Cohort Inputs for Feature Selection
#'
#' Shared input preparation for the feature-selection entry points. Coerces
#' single objects to cohort lists, extracts numeric feature matrices from
#' matrices / data frames / `SummarizedExperiment`s, standardises responses,
#' checks sample-count agreement, and aligns features across cohorts.
#'
#' This is deliberately separate from [.prepare_cohort_inputs()], which serves
#' the optimisation path and additionally handles feature transforms,
#' `feature_subset`, and the richer [.align_features()] strategies.
#'
#' @param x_list A matrix, `SummarizedExperiment`, or list of such objects.
#' @param y_list A binary response aligned with `x_list`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index.
#' @param response Response encoding to return: `"integer"` for 0/1 (the form
#'   `glmnet` and RUV want) or `"factor"` for the standard `No`/`Yes` factor.
#' @param align `"intersection"` keeps only features shared by every cohort and
#'   subsets each matrix to them. `"union"` keeps every feature seen in any
#'   cohort and leaves the matrices ragged — the caller is then responsible for
#'   handling features absent from a given cohort (limma fits per cohort and
#'   pads the resulting statistic matrix with `NA`).
#' @param feature_order `"input"` preserves the column order of the first
#'   cohort; `"sorted"` sorts feature names alphabetically.
#' @return A list with `matrices`, `responses`, `cohort_names`, and
#'   `feature_names`. Under `align = "union"` the matrices are *not* subset to
#'   `feature_names`.
#' @keywords internal
.prepare_selection_inputs <- function(x_list,
                                      y_list,
                                      assay = NULL,
                                      response = c("integer", "factor"),
                                      align = c("intersection", "union"),
                                      feature_order = c("input", "sorted")) {
  response <- match.arg(response)
  align <- match.arg(align)
  feature_order <- match.arg(feature_order)

  x_list <- .as_cohort_list(x_list)
  y_list <- .as_cohort_list(y_list)

  if (length(x_list) != length(y_list)) {
    stop("`x_list` and `y_list` must have the same length.", call. = FALSE)
  }

  cohort_names <- .default_cohort_names(x_list)

  matrices <- vector("list", length(x_list))
  responses <- vector("list", length(x_list))

  for (i in seq_along(x_list)) {
    x_mat <- .ensure_feature_colnames(
      .extract_feature_matrix(x_list[[i]], assay = assay)
    )
    y_vec <- ensure_binary_response(y_list[[i]])

    if (nrow(x_mat) != length(y_vec)) {
      stop("Number of rows in `x_list[[", i, "]]` must match the length of ",
        "`y_list[[", i, "]]`.",
        call. = FALSE
      )
    }

    matrices[[i]] <- x_mat
    responses[[i]] <- if (response == "integer") {
      as.integer(y_vec) - 1L
    } else {
      y_vec
    }
  }

  feature_sets <- lapply(matrices, colnames)

  if (align == "union") {
    ordered_features <- unique(unlist(feature_sets, use.names = FALSE))
    if (feature_order == "sorted") {
      ordered_features <- sort(ordered_features)
    }
    if (!length(ordered_features)) {
      stop(
        "No features were found in any cohort. Check that the input matrices ",
        "have at least one column of numeric data.",
        call. = FALSE
      )
    }
    # Matrices stay ragged: callers under "union" model each cohort separately.
    return(list(
      matrices = matrices,
      responses = responses,
      cohort_names = cohort_names,
      feature_names = ordered_features
    ))
  }

  common_features <- Reduce(intersect, feature_sets)
  if (is.null(common_features) || !length(common_features)) {
    .stop_no_shared_features(feature_sets, cohort_names)
  }

  ordered_features <- if (feature_order == "sorted") {
    sort(common_features)
  } else {
    feature_sets[[1]][feature_sets[[1]] %in% common_features]
  }

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
