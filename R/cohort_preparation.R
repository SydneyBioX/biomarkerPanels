#' Cohort Input Preparation Functions
#'
#' Internal functions for preparing and aligning cohort data before
#' optimization. Handles multi-cohort inputs, feature alignment strategies,
#' and cohort aggregation.
#'
#' @name cohort_preparation
#' @keywords internal
NULL

#' Prepare Cohort Inputs for Optimization
#'
#' Main internal function that processes input data (single or multi-cohort)
#' and prepares it for NSGA-II optimization. Handles feature extraction,
#' alignment, aggregation, and response standardization.
#'
#' @param x Matrix-like object, SummarizedExperiment, or list of such objects.
#' @param y Binary response vector or list of response vectors.
#' @param assay For SummarizedExperiment inputs, assay name or index.
#' @param aggregator Name of the aggregator to apply.
#' @param feature_subset Optional subset of features to keep.
#' @param feature_alignment Alignment strategy for multi-cohort data.
#' @return List with x (combined matrix), truth (factor), cohort (factor),
#'   cohort_names, and cohort_counts.
#' @keywords internal
.prepare_cohort_inputs <- function(x, y, assay = NULL,
                                   aggregator = "none",
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

    matrices <- .apply_cohort_aggregator(matrices, aggregator)

    feature_sets <- lapply(matrices, colnames)
    if (any(vapply(feature_sets, is.null, logical(1)))) {
      stop("Aggregator produced matrices without column names.", call. = FALSE)
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

    x_mat <- .apply_cohort_aggregator(list(x_mat), aggregator)[[1]]
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

    x_mat <- .apply_cohort_aggregator(list(x_mat), aggregator)[[1]]
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

#' Resolve Feature Pool to Feature Names
#'
#' Converts integer indices to feature names and validates that all
#' requested features exist.
#'
#' @param pool Feature identifiers (names or indices).
#' @param feature_names Available feature names.
#' @return Character vector of validated feature names.
#' @keywords internal
.resolve_feature_pool <- function(pool, feature_names) {
  if (is.numeric(pool)) {
    pool <- feature_names[pool]
  }
  if (!all(pool %in% feature_names)) {
    missing <- setdiff(pool, feature_names)
    stop("Feature(s) not found in `x`: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(pool) != length(unique(pool))) {
    duplicates <- pool[duplicated(pool)]
    warning(
      "Feature pool contained duplicates which were removed: ",
      paste(unique(duplicates), collapse = ", "),
      call. = FALSE
    )
  }
  unique(pool)
}

#' Apply Cohort Aggregator to Matrices
#'
#' Applies the specified aggregation function to each cohort matrix.
#'
#' @param matrices List of feature matrices.
#' @param aggregator Name of the aggregator.
#' @return List of aggregated matrices.
#' @keywords internal
.apply_cohort_aggregator <- function(matrices, aggregator) {
  if (!length(matrices)) {
    return(matrices)
  }
  agg_spec <- .get_aggregator(aggregator)
  lapply(matrices, agg_spec$fun)
}

#' Normalize Constraints to Standard Format
#'
#' Processes constraint specifications into a standardized list format.
#'
#' @param constraints List of constraint functions or descriptors.
#' @return List of standardized constraint objects with fun and label.
#' @keywords internal
.normalize_constraints <- function(constraints) {
  if (!length(constraints)) {
    return(list())
  }
  if (!is.list(constraints)) {
    stop("`constraints` must be supplied as a list.", call. = FALSE)
  }
  lapply(seq_along(constraints), function(i) {
    entry <- constraints[[i]]
    label <- names(constraints)[i]
    if (is.function(entry)) {
      fun <- entry
    } else if (is.list(entry) && !is.null(entry$fun) && is.function(entry$fun)) {
      fun <- entry$fun
      if (!is.null(entry$label) && nzchar(entry$label)) {
        label <- entry$label
      }
    } else {
      stop(
        "`constraints` entries must be functions or lists containing a `fun` element.",
        call. = FALSE
      )
    }
    if (is.null(label) || !nzchar(label)) {
      label <- sprintf("constraint_%02d", i)
    }
    list(
      fun = fun,
      label = label
    )
  })
}

#' Get Adaptive NSGA-II Defaults Based on Feature Pool Size
#'
#' Computes appropriate NSGA-II hyperparameters based on the dimensionality
#' of the optimization problem. Larger feature pools require more generations
#' and larger populations to adequately explore the search space.
#'
#' @param n_features Number of features in the decision space.
#' @return Named list of NSGA-II parameters.
#' @keywords internal
.get_adaptive_nsga_defaults <- function(n_features) {
  # Base parameters that don't change with problem size
  # Note: rmoo uses SBX crossover (nc=20) and polynomial mutation (nm=0.2) by default
  base_params <- list(
    pcrossover = 0.7,
    pmutation = 0.2
  )

  # Scale population and generations based on feature pool size
  if (n_features <= 30) {
    adaptive_params <- list(popSize = 64, maxiter = 60)
  } else if (n_features <= 100) {
    adaptive_params <- list(popSize = 128, maxiter = 150)
  } else {
    adaptive_params <- list(popSize = 200, maxiter = 300)
  }

  c(adaptive_params, base_params)
}

#' Align Features Across Cohorts Using Specified Strategy
#'
#' Internal helper to align feature matrices across cohorts using different
#' strategies for handling missing features.
#'
#' @param matrices List of feature matrices (one per cohort).
#' @param feature_sets List of feature name vectors (one per cohort).
#' @param strategy Alignment strategy: "intersection", "majority", or "impute_median".
#' @return List with `matrices` (aligned) and `features` (ordered feature names).
#' @keywords internal
.align_features <- function(matrices, feature_sets, strategy = "intersection") {
  n_cohorts <- length(matrices)

  if (strategy == "intersection") {
    # Original behavior: only features in ALL cohorts
    common_features <- Reduce(intersect, feature_sets)
    if (is.null(common_features) || !length(common_features)) {
      stop(
        "No shared features were found across cohorts. Provide data with ",
        "overlapping feature identifiers, or use feature_alignment = 'majority' ",
        "or 'impute_median' to handle partial overlap.",
        call. = FALSE
      )
    }
    ordered_features <- feature_sets[[1]][feature_sets[[1]] %in% common_features]
    aligned_matrices <- lapply(matrices, function(mat) {
      mat[, ordered_features, drop = FALSE]
    })
    return(list(matrices = aligned_matrices, features = ordered_features))
  }

  # For majority and impute_median, we need to count feature occurrences
  all_features <- unique(unlist(feature_sets))
  feature_counts <- vapply(all_features, function(feat) {
    sum(vapply(feature_sets, function(fs) feat %in% fs, logical(1)))
  }, integer(1))

  if (strategy == "majority") {
    # Keep features in >= 50% of cohorts
    threshold <- ceiling(n_cohorts / 2)
    selected_features <- all_features[feature_counts >= threshold]
    if (!length(selected_features)) {
      stop(
        "No features found in >= 50% of cohorts. Consider using ",
        "feature_alignment = 'intersection' or 'impute_median'.",
        call. = FALSE
      )
    }
  } 
  else if (strategy == "impute_median") {
    # Keep all features from any cohort
    selected_features <- all_features
  } 
  else {
    stop("Unknown feature_alignment strategy: ", strategy, call. = FALSE)
  }

  # Order by frequency (most common first), then alphabetically for ties
  feature_order <- order(-feature_counts[selected_features], selected_features)
  ordered_features <- selected_features[feature_order]

  # Create aligned matrices with imputation for missing features
  aligned_matrices <- lapply(matrices, function(mat) {
    mat_features <- colnames(mat)
    n_samples <- nrow(mat)

    # Initialize result matrix
    result <- matrix(
      NA_real_,
      nrow = n_samples,
      ncol = length(ordered_features),
      dimnames = list(rownames(mat), ordered_features)
    )

    # Copy existing features
    present <- ordered_features[ordered_features %in% mat_features]
    if (length(present)) {
      result[, present] <- mat[, present, drop = FALSE]
    }

    # Impute missing features with cohort-specific median of available features
    # This is more reasonable than 0 since it centers missing features around
    # the cohort's typical expression level
    missing <- ordered_features[!ordered_features %in% mat_features]
    if (length(missing)) {
      warning(
        "Imputing ", length(missing), " missing feature(s) using cohort median. ",
        "Missing: ", paste(head(missing, 5), collapse = ", "),
        if (length(missing) > 5) paste0("... and ", length(missing) - 5, " more") else "",
        ". Consider using feature_alignment = 'intersection' for strict matching.",
        call. = FALSE
      )
      present_features <- ordered_features[ordered_features %in% mat_features]
      if (length(present_features) > 0) {
        cohort_median <- median(result[, present_features, drop = FALSE], na.rm = TRUE)
        if (is.na(cohort_median)) cohort_median <- 0
      } else {
        cohort_median <- 0
      }
      result[, missing] <- cohort_median
    }

    result
  })

  list(matrices = aligned_matrices, features = ordered_features)
}
