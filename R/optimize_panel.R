#' Optimize Biomarker Panels with NSGA-II
#'
#' Wrapper around [rmoo::nsga2()] (or [rmoo::nsga3()]) that composes registered
#' loss functions into a multi-objective search for compact biomarker panels.
#' Candidate solutions are represented as weights over a feature pool; features
#' with weight > 0.5 are included in the panel (up to `max_features`). This
#' threshold-based selection allows variable panel sizes, making the
#' `num_features` objective meaningful. Inputs may be a single cohort (`matrix`,
#' `data.frame`, or `SummarizedExperiment`) or multiple cohorts supplied as
#' lists of such objects. Optimisation should be run on training data only---use
#' [evaluate_panel()] for held-out validation.
#'
#' @param x Matrix-like object, `SummarizedExperiment`, or list of matrices /
#'   experiments representing one or more cohorts.
#' @param y Binary response (`factor`, `character`, or `logical`) aligned with
#'   `x`. When `x` is a list, `y` must be a list of the same length.
#' @param objectives Named list of objective descriptors as returned by
#'   [define_objectives()].
#' @param max_features Maximum number of biomarkers permitted in a panel. Acts
#'   as an upper bound; actual panel size varies based on how many features have
#'   weight > 0.5 in the decision vector. Minimum is 2 features when
#'   `regularized = TRUE` (glmnet requirement), or 1 when `regularized = FALSE`.
#' @param feature_pool Optional subset of base feature identifiers (names or
#'   integer indices) to consider during optimization. These are the original
#'   feature names before any transformation. When using pairwise transforms,
#'   transformed labels like `"A--B"` are also accepted and will be mapped
#'   back to their constituent base features. Defaults to all features.
#' @param feature_transform Transformation applied to selected features
#'   after base feature selection. The NSGA algorithm selects base features,
#'   then this transform is applied on-the-fly. Defaults to `"pairwise_ratios"`,
#'   which generates pairwise contrasts via [`pairwise_col_diff()`] to dampen
#'   batch effects. Other options include `"pairwise_log_ratios"` (log-ratios
#'   for compositional data), `"reference_norm"` (normalize relative to a
#'   reference feature), and `"none"` (no transformation). Use
#'   [`register_feature_transform()`] to add custom strategies. See
#'   [`feature_transform_registry()`] for available options.
#' @param feature_alignment Strategy for aligning features across cohorts:
#'   \describe{
#'     \item{"intersection"}{(Default) Keep only features present in ALL cohorts.
#'       Features missing in any cohort are dropped.}
#'     \item{"majority"}{Keep features present in >= 50% of cohorts.
#'       Missing values are imputed with cohort medians.}
#'     \item{"impute_median"}{Keep all features from any cohort.
#'       Missing values are imputed with cohort medians.}
#'   }
#' @param constraints Optional list of constraint descriptors (e.g.,
#'   from [min_metric_constraint()]) that must evaluate to `TRUE` for a candidate
#'   solution to be considered feasible.
#' @param scoring_fn Function producing per-sample scores from the selected
#'   features. Signature:
#'   `function(x_selected, selected_features, truth, cohort = NULL, ...)`.
#' @param algorithm Multi-objective optimization algorithm. `"NSGA-II"` (default)
#'   uses crowding distance for diversity preservation. `"NSGA-III"` uses
#'   reference-point-based selection for many-objective problems but has a bug
#'   in rmoo 0.3.0 (use NSGA-II as workaround).
#' @param nsga_control Named list of arguments passed to [rmoo::nsga2()] (or
#'   [rmoo::nsga3()]). For NSGA-III, `n_partitions` controls reference point
#'   density (default computed from number of objectives). Defaults to adaptive
#'   values based on feature pool size. Note: `parallel = TRUE` is not recommended
#'   as it is slower than sequential execution due to overhead.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @param seed Optional integer seed for reproducibility. When provided, sets
#'   the random seed before running NSGA-II optimization. This ensures
#'   reproducible results across runs. If `NULL` (default), no seed is set and
#'   results may vary between runs.
#' @param fitness_cv Logical; if `TRUE` (default), use cross-validation when
#'   evaluating candidate solutions during NSGA-II optimization. This prevents
#'   overfitting by computing objective metrics on held-out fold predictions
#'   rather than in-sample predictions. Recommended for fair comparison with
#'   regularized methods like Lasso.
#' @param fitness_cv_folds Number of cross-validation folds for fitness
#'   evaluation when `fitness_cv = TRUE`. Default is 5.
#' @param regularized Logical; if `TRUE` (default), use regularized regression
#'   (elastic net via glmnet) for scoring candidates during optimization. This
#'   typically produces more stable fitness estimates and better out-of-sample
#'   performance compared to unregularized logistic regression.
#' @param regularized_alpha Elastic net mixing parameter when `regularized = TRUE`.
#'   `alpha = 1` is lasso, `alpha = 0` is ridge, and values in between give
#'   elastic net. Default is 0.5 (elastic net).
#' @param selection_threshold Either a fixed numeric threshold in (0,1) for
#'   feature selection, or `"adaptive"` (default) to use gap-based selection
#'   that allows variable panel sizes in the Pareto front. Fixed thresholds
#'   (e.g., 0.5) tend to produce panels at `max_features` because more features

#'   improve accuracy. Adaptive selection finds natural breakpoints in the
#'   weight distribution, enabling the `num_features` objective to drive panel
#'   size diversity.
#' @return An `OptimizationResult` containing the Pareto-optimal solutions.
#'   Use [summarize_solutions()] to inspect solutions and [fit_panel()] to
#'   fit a model on a selected solution.
#' @export
#' @seealso [fit_panel()], [summarize_solutions()], [evaluate_panel()]
optimize_panel <- function(x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "num_features")
                           ),
                           max_features = 5L,
                           feature_pool = NULL,
                           feature_transform = "pairwise_ratios",
                           feature_alignment = "intersection",
                           constraints = list(),
                           scoring_fn = NULL,
                           algorithm = c("NSGA-II", "NSGA-III"),
                           nsga_control = list(),
                           assay = NULL,
                           seed = NULL,
                           fitness_cv = TRUE,
                           fitness_cv_folds = 5L,
                           regularized = TRUE,
                           regularized_alpha = 0.5,
                           selection_threshold = "adaptive") {
  algorithm <- match.arg(algorithm)
  feature_alignment <- match.arg(feature_alignment)
  if (!requireNamespace("rmoo", quietly = TRUE)) {
    stop("The 'rmoo' package is required. Install it via install.packages('rmoo').",
         call. = FALSE)
  }

  if (!is.character(feature_transform) || length(feature_transform) != 1L ||
      !nzchar(feature_transform)) {
    stop("`feature_transform` must be a single non-empty character string.", call. = FALSE)
  }
  if (!exists(feature_transform, envir = .transform_registry, inherits = FALSE)) {
    available <- ls(.transform_registry)
    stop(
      "Unknown feature transform '", feature_transform, "'. ",
      "Available: ", paste(available, collapse = ", "), ". ",
      "Use register_feature_transform() to add custom transforms.",
      call. = FALSE
    )
  }

  inputs_raw <- .prepare_cohort_inputs(x, y, assay = assay, transform = "none",
                                       feature_alignment = feature_alignment)
  raw_x_mat <- inputs_raw$x
  raw_feature_names <- colnames(raw_x_mat)
  if (is.null(raw_feature_names)) {
    raw_feature_names <- sprintf("feature_%04d", seq_len(ncol(raw_x_mat)))
    colnames(raw_x_mat) <- raw_feature_names
  }

  feature_pool_arg <- feature_pool
  feature_pool_base <- raw_feature_names
  feature_pool_final <- NULL

  if (!is.null(feature_pool_arg)) {
    if (is.numeric(feature_pool_arg)) {
      feature_pool_base <- .resolve_feature_pool(feature_pool_arg, raw_feature_names)
    } else {
      feature_pool_arg <- unique(feature_pool_arg)
      missing_raw <- setdiff(feature_pool_arg, raw_feature_names)
      if (!length(missing_raw)) {
        feature_pool_base <- feature_pool_arg
      } else if (
        feature_transform != "none" &&
        all(grepl("--", feature_pool_arg, fixed = TRUE))
      ) {
        components <- unique(unlist(strsplit(feature_pool_arg, "--", fixed = TRUE)))
        missing_components <- setdiff(components, raw_feature_names)
        if (length(missing_components)) {
          stop(
            "Feature(s) referenced in `feature_pool` not found in training data: ",
            paste(missing_components, collapse = ", "),
            call. = FALSE
          )
        }
        feature_pool_base <- components
        feature_pool_final <- feature_pool_arg
      } else {
        missing <- setdiff(feature_pool_arg, raw_feature_names)
        stop(
          "Feature(s) not found in `x`: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  # Note about pairwise transforms (no longer causes O(n^2) search space)
  # With base-feature-first selection, search space is O(n), not O(n^2)
  if (feature_transform %in% c("pairwise_ratios", "pairwise_log_ratios")) {
    n_base_features <- length(feature_pool_base)
    if (n_base_features > 200) {
      message(
        sprintf(
          "Note: Searching %d base features with pairwise transform.",
          n_base_features
        )
      )
    }
  }

  # Prepare raw (untransformed) data - transform is applied on-the-fly during fitness evaluation
  inputs <- .prepare_cohort_inputs(
    x,
    y,
    assay = assay,
    transform = "none",  # Always raw - transform applied after feature selection
    feature_subset = feature_pool_base,
    feature_alignment = feature_alignment
  )
  x_raw <- inputs$x
  truth <- inputs$truth
  cohort <- inputs$cohort

  # Feature pool is now BASE features (not aggregated)
  feature_pool <- feature_pool_base

  if (!length(feature_pool)) {
    stop("`feature_pool` produced zero features.", call. = FALSE)
  }

  # Minimum features required depends on transform and regularization

  # With pairwise transforms: n base features -> n*(n-1)/2 ratios
  # For regularized fitting: need >= 2 transformed features
  # 2 base features -> 1 ratio (not enough for glmnet)
  # 3 base features -> 3 ratios (sufficient for glmnet)
  if (regularized && feature_transform %in% c("pairwise_ratios", "pairwise_log_ratios")) {
    min_features_for_regularized <- 3L
  } else if (regularized) {
    min_features_for_regularized <- 2L
  } else {
    min_features_for_regularized <- 1L
  }

  if (max_features < min_features_for_regularized) {
    stop("`max_features` must be at least ", min_features_for_regularized,
         " when `regularized = TRUE`",
         if (feature_transform %in% c("pairwise_ratios", "pairwise_log_ratios"))
           paste0(" with feature_transform = '", feature_transform, "'") else "",
         ". Either increase `max_features` or set `regularized = FALSE`.",
         call. = FALSE)
  }

  if (max_features < 1L) {
    stop("`max_features` must be at least 1.", call. = FALSE)
  }

  if (max_features > length(feature_pool)) {
    max_features <- length(feature_pool)
  }

  # x_pool contains RAW (base) features - transform applied on-the-fly
  x_pool <- x_raw[, feature_pool, drop = FALSE]
  decision_dim <- ncol(x_pool)  # Decision dimension = number of BASE features

  if (decision_dim > 200) {
    warning("Optimizing over more than 200 features may be slow; consider ",
            "reducing `feature_pool` for exploration.")
  }

  if (is.null(scoring_fn)) {
    scoring_fn <- .default_scoring_fn
  }

  if (!is.function(scoring_fn)) {
    stop("`scoring_fn` must be a function.", call. = FALSE)
  }

  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  names(objective_directions) <- names(objectives)

  constraint_specs <- .normalize_constraints(constraints)

  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim, algorithm)
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)


  # Create CV folds for fitness evaluation (if enabled)
  cv_fold_ids <- NULL
  if (fitness_cv) {
    n_samples <- nrow(x_pool)
    if (n_samples < fitness_cv_folds * 2L) {
      warning("Too few samples for ", fitness_cv_folds, "-fold CV in fitness evaluation. ",
              "Falling back to in-sample scoring.", call. = FALSE)
      fitness_cv <- FALSE
    } else {
      cv_fold_ids <- .create_stratified_folds(truth, fitness_cv_folds)
    }
  }

  # Minimum base features required depends on transform and regularization
  # min_features_for_regularized is already computed above
  min_features_required <- min_features_for_regularized

  evaluate_candidate <- function(decision_vec) {
    # feature_names_pool contains BASE feature names
    feature_names_pool <- colnames(x_pool)
    n_pool <- length(feature_names_pool)

    # Step 1: Select BASE features using decision weights
    if (identical(selection_threshold, "adaptive")) {
      # Adaptive: select based on largest gap in sorted weights
      # This allows the GA to "encode" panel size in the weight distribution
      ord <- order(decision_vec, feature_names_pool,
                   decreasing = c(TRUE, FALSE), method = "radix")
      sorted_weights <- decision_vec[ord]

      if (n_pool > max_features) {
        top_n <- min(max_features + 5L, n_pool)
        weight_diffs <- -diff(sorted_weights[seq_len(top_n)])

        # Find largest gap after min_features_required
        search_range <- seq(min_features_required, top_n - 1L)
        if (length(search_range) > 0L) {
          gap_idx <- which.max(weight_diffs[search_range])
          natural_cutoff <- min_features_required + gap_idx
        } else {
          natural_cutoff <- min_features_required
        }
        n_selected <- min(natural_cutoff, max_features)
        n_selected <- max(n_selected, min_features_required)
      } else {
        above_median <- sum(decision_vec > 0.5)
        n_selected <- max(min_features_required, min(above_median, max_features))
      }
      selected_idx <- ord[seq_len(n_selected)]

    } else {
      # Fixed threshold (backward compatible)
      threshold <- as.numeric(selection_threshold)
      above_threshold <- which(decision_vec > threshold)

      if (length(above_threshold) < min_features_required) {
        # Fallback: take top min_features_required by weight
        ord <- order(decision_vec, feature_names_pool,
                     decreasing = c(TRUE, FALSE), method = "radix")
        selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
      } else if (length(above_threshold) > max_features) {
        # Cap at max_features, using feature names for reproducible tie-breaking
        weights_above <- decision_vec[above_threshold]
        names_above <- feature_names_pool[above_threshold]
        ord <- order(weights_above, names_above,
                     decreasing = c(TRUE, FALSE), method = "radix")
        selected_idx <- above_threshold[ord[seq_len(max_features)]]
      } else {
        selected_idx <- above_threshold
      }
    }

    # Sort by descending weight for consistent output ordering
    selected_idx <- selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
    selected_base_features <- feature_names_pool[selected_idx]
    x_base_selected <- x_pool[, selected_base_features, drop = FALSE]

    # Step 2: Apply feature transform ON-THE-FLY to selected base features
    if (feature_transform != "none" && ncol(x_base_selected) >= 2L) {
      x_transformed <- .apply_feature_transform_single(x_base_selected, feature_transform)
      selected_features <- colnames(x_transformed)  # e.g., "A--B" for pairwise
      x_selected <- x_transformed
    } else {
      # No transform or single feature - use base features directly
      x_selected <- x_base_selected
      selected_features <- selected_base_features
    }

    # Step 3: Compute scores using CV (out-of-fold predictions) or in-sample
    if (fitness_cv && !is.null(cv_fold_ids)) {
      # Cross-validation based scoring: prevents overfitting
      scores <- .compute_cv_scores(
        x_selected = x_selected,
        truth = truth,
        fold_ids = cv_fold_ids,
        cohort = cohort,
        regularized = regularized,
        alpha = regularized_alpha
      )
    } else {
      # In-sample scoring (for backward compatibility or small datasets)
      if (regularized || identical(scoring_fn, .default_scoring_fn)) {
        scores <- .default_scoring_fn(
          x_selected = x_selected,
          selected_features = selected_features,
          truth = truth,
          cohort = cohort,
          regularized = regularized,
          alpha = regularized_alpha
        )
      } else {
        # Custom scoring function provided by user (only used when regularized = FALSE)
        score_args <- list(
          x_selected = x_selected,
          selected_features = selected_features,
          truth = truth
        )
        if (!is.null(cohort)) {
          score_args$cohort <- cohort
        }
        scores <- do.call(scoring_fn, score_args)
      }
    }

    if (!is.numeric(scores) || length(scores) != nrow(x_pool)) {
      stop("Scoring must return a numeric vector matching the number of samples.",
           call. = FALSE)
    }

    # Step 4: Evaluate constraints on TRANSFORMED features
    constraint_results <- if (length(constraint_specs)) {
      setNames(vapply(seq_along(constraint_specs), function(j) {
        res <- constraint_specs[[j]]$fun(
          truth = truth,
          scores = scores,
          selected = selected_features,
          cohort = cohort,
          x = x_selected
        )
        isTRUE(res)
      }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
    } else {
      logical(0)
    }
    feasible <- if (length(constraint_results)) all(constraint_results) else TRUE

    # Step 5: Evaluate objectives on TRANSFORMED features
    metrics <- vapply(objectives, function(obj) {
      obj$fun(
        truth,
        scores,
        selected = selected_features,
        cohort = cohort,
        x = x_selected
      )
    }, numeric(1))

    # Return BOTH base features and transformed features
    list(
      base_features = selected_base_features,  # Original genes selected
      features = selected_features,             # Transformed features (for model)
      scores = scores,
      metrics = metrics,
      constraint_results = constraint_results,
      feasible = feasible
    )
  }

  # Large finite penalty instead of Inf — NSGA-III normalization produces NaN
  # from Inf values, causing "missing value where TRUE/FALSE needed" errors
  .PENALTY <- 1e6

  # Evaluate a single decision vector and return converted objective values
  evaluate_single <- function(decision_vec) {
    evaluated <- evaluate_candidate(decision_vec)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(.PENALTY, length(objectives)))
    }
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      # Handle NA, Inf, and -Inf explicitly
      if (is.na(val)) {
        return(.PENALTY)  # Treat NA as infeasible
      }
      if (is.infinite(val)) {
        return(.PENALTY)
      }
      if (dir == "maximize") {
        return(-val)
      }
      val
    }, val = metrics, dir = objective_directions, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    as.numeric(converted)
  }

  # rmoo fitness function: can receive matrix (rows=individuals) or vector (single individual)
  # Must return matrix when receiving matrix input
  # Note: Accept ... to handle rmoo bug that passes reference_dirs to fitness
  objective_wrapper <- function(x, ...) {
    if (is.null(dim(x))) {
      # Single individual as vector
      return(evaluate_single(x))
    }
    # Matrix input: evaluate each row and return matrix of objectives
    t(apply(x, 1, evaluate_single))
  }

  nsga_params <- c(
    list(
      type = "real-valued",
      fitness = objective_wrapper,
      nObj = length(objectives),
      lower = rep(0, decision_dim),
      upper = rep(1, decision_dim),
      monitor = FALSE,
      summary = FALSE
    ),
    nsga_args
  )

  # Add reference points for NSGA-III if not user-specified
  if (algorithm == "NSGA-III" && is.null(nsga_args$n_partitions)) {
    nsga_params$n_partitions <- .compute_nsga3_partitions(length(objectives))
  }

  # Set random seed for reproducibility if provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("`seed` must be a single integer value.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  # Generate sparse initialization suggestions to encourage diverse panel sizes
  # This seeds the initial population with solutions spanning min to max features
  n_suggestions <- min(20L, nsga_params$popSize %/% 4L)
  if (n_suggestions >= 2L) {
    sparse_suggestions <- .generate_sparse_suggestions(
      n_features = decision_dim,
      n_suggestions = n_suggestions,
      min_features = min_features_required,
      max_features = max_features,
      seed = seed
    )
    nsga_params$suggestions <- sparse_suggestions
  }

  if (algorithm == "NSGA-II") {
    nsga_result <- do.call(rmoo::nsga2, nsga_params)
  } else {
    nsga_result <- do.call(rmoo::nsga3, nsga_params)
  }

  # Filter to Pareto-optimal solutions (front rank == 1)
  # rmoo returns full population; mco returned only Pareto-optimal
  optimal_idx <- which(nsga_result@front == 1)
  pareto_pop <- nsga_result@population[optimal_idx, , drop = FALSE]

  # Handle edge case: single solution returns vector instead of matrix
  if (is.null(dim(pareto_pop))) {
    pareto_pop <- matrix(pareto_pop, nrow = 1)
  }

  solutions <- lapply(seq_len(nrow(pareto_pop)), function(i) {
    decision_vec <- pareto_pop[i, ]
    evaluate_candidate(decision_vec)
  })

  feasible_vec <- vapply(solutions, `[[`, logical(1), "feasible")
  if (!any(feasible_vec)) {
    stop("No solutions satisfied the supplied constraints. Consider relaxing them.",
         call. = FALSE)
  }

  solutions <- solutions[feasible_vec]
  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  # Build solutions data frame in wide format (one row per solution)
  solutions_df <- data.frame(
    solution_id = seq_along(solutions),
    stringsAsFactors = FALSE
  )

  # Add base_features and features as list columns
  solutions_df$base_features <- I(lapply(solutions, `[[`, "base_features"))
  solutions_df$features <- I(lapply(solutions, `[[`, "features"))

  # Add objective values as separate columns
  for (obj_name in names(objectives)) {
    solutions_df[[obj_name]] <- metric_matrix[, obj_name]
  }

  # Build control parameters to store
  control <- list(
    max_features = max_features,
    feature_pool = feature_pool,
    base_feature_pool = feature_pool_base,
    algorithm = algorithm,
    nsga2 = nsga_args,  # Keep nsga2 name for backward compatibility
    scoring_function = deparse(substitute(scoring_fn)),
    feature_transform = feature_transform,
    feature_alignment = feature_alignment,
    constraints = if (length(constraint_specs)) {
      vapply(constraint_specs, `[[`, character(1), "label")
    } else {
      character()
    },
    # Store positive class for consistent evaluation
    # ensure_binary_response() standardizes to levels c("No", "Yes")
    positive_class = "Yes",
    response_levels = levels(truth),
    # Store seed for reproducibility documentation
    seed = seed,
    selection_threshold = selection_threshold,
    regularized = regularized,
    regularized_alpha = if (regularized) regularized_alpha else NULL,
    objective_directions = objective_directions
  )

  # Build training signature
  training_signature <- list(
    n = nrow(x_raw),
    p = ncol(x_raw),
    class_balance = table(truth),
    feature_pool_size = length(feature_pool),
    base_feature_pool_size = length(feature_pool_base),
    num_cohorts = length(inputs$cohort_names),
    cohort_labels = inputs$cohort_names,
    cohort_counts = inputs$cohort_counts
  )

  # Return OptimizationResult (no model fitting)
  # Note: aggregated_x now stores RAW (untransformed) base features
  # Transform is applied on-the-fly during fit_panel() and evaluate_panel()
  new(
    "OptimizationResult",
    solutions = solutions_df,
    feature_pool = feature_pool,
    control = control,
    training_signature = training_signature,
    aggregated_x = x_pool,  # Raw base features (transform applied later)
    aggregated_y = truth,
    aggregated_cohort = cohort
  )
}
