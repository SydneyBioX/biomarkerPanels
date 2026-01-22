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

#' Get Adaptive NSGA Defaults Based on Feature Pool Size and Algorithm
#'
#' Computes appropriate NSGA-II or NSGA-III hyperparameters based on the
#' dimensionality of the optimization problem and the selected algorithm.
#' Larger feature pools require more generations and larger populations to
#' adequately explore the search space. NSGA-III typically benefits from
#' larger populations for good reference point coverage.
#'
#' @param n_features Number of features in the decision space.
#' @param algorithm Algorithm to use: `"NSGA-III"` or `"NSGA-II"`.
#' @return Named list of NSGA parameters.
#' @keywords internal
.get_adaptive_nsga_defaults <- function(n_features, algorithm = "NSGA-III") {
  # Base parameters that don't change with problem size
  # Note: rmoo uses SBX crossover (nc=20) and polynomial mutation (nm=0.2) by default
  base_params <- list(
    pcrossover = 0.7,
    pmutation = 0.2
  )

  # NSGA-III typically needs larger populations for good reference point coverage
  if (algorithm == "NSGA-III") {
    if (n_features <= 30) {
      adaptive_params <- list(popSize = 92, maxiter = 80)
    } else if (n_features <= 100) {
      adaptive_params <- list(popSize = 156, maxiter = 180)
    } else {
      adaptive_params <- list(popSize = 252, maxiter = 350)
    }
  } else {
    # NSGA-II defaults (existing)
    if (n_features <= 30) {
      adaptive_params <- list(popSize = 64, maxiter = 60)
    } else if (n_features <= 100) {
      adaptive_params <- list(popSize = 128, maxiter = 150)
    } else {
      adaptive_params <- list(popSize = 200, maxiter = 300)
    }
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

#' Compute NSGA-III Partition Count Based on Number of Objectives
#'
#' Determines the appropriate number of partitions for generating reference
#' points in NSGA-III based on the number of objectives. More objectives
#' require fewer partitions to avoid exponential growth in reference points.
#'
#' @param n_objectives Number of objectives in the optimization problem.
#' @return Integer number of partitions.
#' @keywords internal
.compute_nsga3_partitions <- function(n_objectives) {

  # Rule of thumb: partitions scale inversely with objectives
  # For 2-3 objectives: 12 partitions
  # For 4-5 objectives: 6 partitions
  # For 6+ objectives: 4 partitions
  if (n_objectives <= 3) 12L
  else if (n_objectives <= 5) 6L
  else 4L
}

#' Generate Sparse Initial Population Suggestions
#'
#' Creates a diverse set of initial weight vectors that span different panel
#' sizes, from `min_features` to `max_features`. This helps NSGA explore
#' the full range of panel sizes rather than converging only toward max_features.
#'
#' @param n_features Number of features in the decision space.
#' @param n_suggestions Number of suggestion vectors to generate.
#' @param min_features Minimum number of features to include (active weights).
#' @param max_features Maximum number of features to include.
#' @param seed Optional integer seed for reproducibility.
#' @return Matrix with `n_suggestions` rows and `n_features` columns, where
#'   each row is a weight vector with varying numbers of high/low values.
#' @keywords internal
.generate_sparse_suggestions <- function(n_features,
                                          n_suggestions = 20L,
                                          min_features = 2L,
                                          max_features = 10L,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  target_sizes <- unique(round(seq(min_features, max_features,
                                    length.out = n_suggestions)))
  n_suggestions <- length(target_sizes)
  suggestions <- matrix(0, nrow = n_suggestions, ncol = n_features)

  for (i in seq_len(n_suggestions)) {
    k <- target_sizes[i]
    on_idx <- sample(n_features, size = min(k, n_features))
    suggestions[i, on_idx] <- runif(length(on_idx), 0.7, 0.95)
    off_idx <- setdiff(seq_len(n_features), on_idx)
    suggestions[i, off_idx] <- runif(length(off_idx), 0.05, 0.3)
  }
  suggestions
}

#' Validate Partition Ratios for Train/Validation/Held-out Split
#'
#' Ensures that train, validation, and held-out ratios satisfy minimum
#' requirements for meaningful model fitting and evaluation.
#'
#' @param train_ratio Proportion of data for training (must be >= 0.5).
#' @param val_ratio Proportion of data for validation (must be >= 0.1).
#' @return Invisibly returns TRUE if valid; otherwise throws an error.
#' @keywords internal
.validate_partition_ratios <- function(train_ratio, val_ratio) {
  if (!is.numeric(train_ratio) || length(train_ratio) != 1L ||
      is.na(train_ratio) || train_ratio < 0 || train_ratio > 1) {
    stop("`train_ratio` must be a single numeric value between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(val_ratio) || length(val_ratio) != 1L ||
      is.na(val_ratio) || val_ratio < 0 || val_ratio > 1) {
    stop("`val_ratio` must be a single numeric value between 0 and 1.",
         call. = FALSE)
  }

  heldout_ratio <- 1 - train_ratio - val_ratio

  if (train_ratio < 0.5) {
    stop("`train_ratio` must be at least 0.5 to ensure sufficient training data.",
         call. = FALSE)
  }
  if (val_ratio < 0.1) {
    stop("`val_ratio` must be at least 0.1 to ensure meaningful validation.",
         call. = FALSE)
  }
  if (heldout_ratio < 0.05) {
    stop("Held-out ratio (1 - train_ratio - val_ratio) must be at least 0.05. ",
         "Current: ", round(heldout_ratio, 3), call. = FALSE)
  }

  invisible(TRUE)
}

#' Stratified Partitioning of Multi-Cohort Data
#'
#' Partitions each cohort's data into train, validation, and held-out sets
#' while maintaining class balance within each partition.
#'
#' @param x_list List of feature matrices (one per cohort).
#' @param y_list List of binary response factors (one per cohort).
#' @param train_ratio Proportion of data for training.
#' @param val_ratio Proportion of data for validation.
#' @return List containing:
#'   \describe{
#'     \item{train_x}{List of training matrices (one per cohort)}
#'     \item{train_y}{List of training response factors}
#'     \item{val_x}{List of validation matrices}
#'     \item{val_y}{List of validation response factors}
#'     \item{heldout_x}{List of held-out matrices}
#'     \item{heldout_y}{List of held-out response factors}
#'     \item{cohort_names}{Character vector of cohort names}
#'     \item{partition_info}{Data frame with partition sizes per cohort}
#'   }
#' @keywords internal
.stratified_partition_cohorts <- function(x_list, y_list, train_ratio, val_ratio) {
  if (!is.list(x_list) || !is.list(y_list)) {
    stop("`x_list` and `y_list` must be lists.", call. = FALSE)
  }
  if (length(x_list) != length(y_list)) {
    stop("`x_list` and `y_list` must have the same length.", call. = FALSE)
  }

  n_cohorts <- length(x_list)
  cohort_names <- names(x_list)
  if (is.null(cohort_names) || any(cohort_names == "")) {
    cohort_names <- sprintf("cohort_%02d", seq_len(n_cohorts))
  }

  train_x <- vector("list", n_cohorts)
  train_y <- vector("list", n_cohorts)
  val_x <- vector("list", n_cohorts)
  val_y <- vector("list", n_cohorts)
  heldout_x <- vector("list", n_cohorts)
  heldout_y <- vector("list", n_cohorts)

  partition_info <- data.frame(
    cohort = cohort_names,
    n_total = integer(n_cohorts),
    n_train = integer(n_cohorts),
    n_val = integer(n_cohorts),
    n_heldout = integer(n_cohorts),
    train_yes = integer(n_cohorts),
    train_no = integer(n_cohorts),
    val_yes = integer(n_cohorts),
    val_no = integer(n_cohorts),
    heldout_yes = integer(n_cohorts),
    heldout_no = integer(n_cohorts),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_cohorts)) {
    x_i <- x_list[[i]]
    y_i <- ensure_binary_response(y_list[[i]])

    n <- nrow(x_i)
    if (n != length(y_i)) {
      stop("Cohort ", i, ": `x` and `y` have different sample sizes.", call. = FALSE)
    }

    # Stratified sampling within each class
    yes_idx <- which(y_i == "Yes")
    no_idx <- which(y_i == "No")

    # Partition each class
    partition_class <- function(idx) {
      n_class <- length(idx)
      if (n_class == 0L) {
        return(list(train = integer(0), val = integer(0), heldout = integer(0)))
      }

      shuffled <- sample(idx)
      n_train <- max(1L, round(n_class * train_ratio))
      n_val <- max(1L, round(n_class * val_ratio))
      n_heldout <- n_class - n_train - n_val

      # Ensure at least 1 sample in held-out if possible

      if (n_heldout < 1L && n_class >= 3L) {
        n_heldout <- 1L
        if (n_train > n_val) {
          n_train <- n_train - 1L
        } else {
          n_val <- n_val - 1L
        }
      }

      list(
        train = shuffled[seq_len(n_train)],
        val = if (n_val > 0L) shuffled[(n_train + 1L):min(n_train + n_val, n_class)] else integer(0),
        heldout = if (n_heldout > 0L) shuffled[(n_train + n_val + 1L):n_class] else integer(0)
      )
    }

    yes_parts <- partition_class(yes_idx)
    no_parts <- partition_class(no_idx)

    train_idx <- c(yes_parts$train, no_parts$train)
    val_idx <- c(yes_parts$val, no_parts$val)
    heldout_idx <- c(yes_parts$heldout, no_parts$heldout)

    # Check for small partitions
    min_per_class <- 2L
    warn_threshold <- 5L

    check_partition_size <- function(part_name, yes_n, no_n) {
      if (yes_n < min_per_class || no_n < min_per_class) {
        stop(
          "Cohort '", cohort_names[i], "' has insufficient samples in ", part_name, " set. ",
          "Yes: ", yes_n, ", No: ", no_n, ". Minimum required: ", min_per_class, " per class.",
          call. = FALSE
        )
      }
      if (yes_n < warn_threshold || no_n < warn_threshold) {
        warning(
          "Cohort '", cohort_names[i], "' has few samples in ", part_name, " set. ",
          "Yes: ", yes_n, ", No: ", no_n, ". Results may be unreliable.",
          call. = FALSE
        )
      }
    }

    check_partition_size("training", length(yes_parts$train), length(no_parts$train))
    check_partition_size("validation", length(yes_parts$val), length(no_parts$val))
    check_partition_size("held-out", length(yes_parts$heldout), length(no_parts$heldout))

    # Store partitions
    train_x[[i]] <- x_i[train_idx, , drop = FALSE]
    train_y[[i]] <- y_i[train_idx]
    val_x[[i]] <- x_i[val_idx, , drop = FALSE]
    val_y[[i]] <- y_i[val_idx]
    heldout_x[[i]] <- x_i[heldout_idx, , drop = FALSE]
    heldout_y[[i]] <- y_i[heldout_idx]

    # Record partition info
    partition_info$n_total[i] <- n
    partition_info$n_train[i] <- length(train_idx)
    partition_info$n_val[i] <- length(val_idx)
    partition_info$n_heldout[i] <- length(heldout_idx)
    partition_info$train_yes[i] <- length(yes_parts$train)
    partition_info$train_no[i] <- length(no_parts$train)
    partition_info$val_yes[i] <- length(yes_parts$val)
    partition_info$val_no[i] <- length(no_parts$val)
    partition_info$heldout_yes[i] <- length(yes_parts$heldout)
    partition_info$heldout_no[i] <- length(no_parts$heldout)
  }

  names(train_x) <- cohort_names
  names(train_y) <- cohort_names
  names(val_x) <- cohort_names
  names(val_y) <- cohort_names
  names(heldout_x) <- cohort_names
  names(heldout_y) <- cohort_names

  list(
    train_x = train_x,
    train_y = train_y,
    val_x = val_x,
    val_y = val_y,
    heldout_x = heldout_x,
    heldout_y = heldout_y,
    cohort_names = cohort_names,
    partition_info = partition_info
  )
}
