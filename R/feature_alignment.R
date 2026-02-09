#' Feature Alignment and Transformation
#'
#' Internal functions for aligning features across cohorts and applying
#' feature transformations.
#'
#' @name feature_alignment
#' @keywords internal
NULL

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

#' Apply Feature Transform to Matrices
#'
#' Applies the specified transform function to each cohort matrix.
#'
#' @param matrices List of feature matrices.
#' @param transform Name of the feature transform.
#' @return List of transformed matrices.
#' @keywords internal
.apply_feature_transform <- function(matrices, transform) {
  if (!length(matrices)) {
    return(matrices)
  }
  transform_spec <- .get_feature_transform(transform)
  lapply(matrices, transform_spec$fun)
}

#' Apply Feature Transform to Single Matrix
#'
#' Applies the specified transform function to a single matrix.
#' Used during optimization to transform selected base features on-the-fly.
#'
#' @param x A numeric matrix with column names.
#' @param transform_name Name of the feature transform.
#' @return Transformed matrix.
#' @keywords internal
.apply_feature_transform_single <- function(x, transform_name) {
  if (transform_name == "none" || ncol(x) < 2L) {
    return(x)
  }
  transform_spec <- .get_feature_transform(transform_name)
  transform_spec$fun(x)
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
