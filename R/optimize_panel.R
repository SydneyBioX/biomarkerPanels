#' Optimize Biomarker Panels with NSGA-II
#'
#' Wrapper around [rmoo::nsga2()] that composes registered loss functions into a
#' multi-objective search for compact biomarker panels. Candidate solutions are
#' represented as weights over a feature pool; features with weight > 0.5 are
#' included in the panel (up to `max_features`). This threshold-based selection
#' allows variable panel sizes, making the `num_features` objective meaningful.
#' Inputs may be a single cohort (`matrix`, `data.frame`, or
#' `SummarizedExperiment`) or multiple cohorts supplied as lists of such objects.
#' Optimisation should be run on training data only---use [evaluate_panel()] for
#' held-out validation.
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
#' @param feature_pool Optional subset of feature identifiers (names or integer
#'   indices) considered during optimization. When a cohort aggregator is used,
#'   specify the underlying (pre-aggregation) feature names; aggregated labels
#'   such as `"A--B"` are also accepted and will be mapped back to their
#'   constituents. Defaults to all features.
#' @param cohort_aggregator Transformation applied to cohort feature matrices
#'   prior to alignment. Defaults to `"pairwise_ratios"`, which generates
#'   pairwise within-cohort contrasts via [`pairwise_col_diff()`] to dampen
#'   distributional shifts across sites. Other built-in options include
#'   `"pairwise_log_ratios"` (log-ratios for compositional data),
#'   `"reference_norm"` (normalize relative to a reference feature), and
#'   `"none"` (no transformation). Use [`register_aggregator()`] to add custom
#'   strategies. See [`aggregator_registry()`] for available options.
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
#' @param nsga_control Named list of arguments passed to [rmoo::nsga2()]. Defaults
#'   to `list(popSize = 64, maxiter = 60, pcrossover = 0.7, pmutation = 0.2)`.
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
#' @param final_model_cv Logical; if `TRUE`, use nested cross-validation when
#'   fitting the final model to reduce training data contamination. The final
#'   model coefficients are averaged across CV folds. Default is `FALSE` for
#'   backward compatibility.
#' @param cv_folds Number of cross-validation folds when `final_model_cv = TRUE`.
#'   Default is 5.
#' @param regularized Logical; if `TRUE` (default), use regularized regression
#'   (elastic net via glmnet) for scoring candidates during optimization. This
#'   typically produces more stable fitness estimates and better out-of-sample
#'   performance compared to unregularized logistic regression.
#' @param regularized_alpha Elastic net mixing parameter when `regularized = TRUE`.
#'   `alpha = 1` is lasso, `alpha = 0` is ridge, and values in between give
#'   elastic net. Default is 0.5 (elastic net).
#' @return A `BiomarkerPanelResult` with Pareto-optimal solutions summarised in
#'   the `objectives` slot.
#' @export
optimize_panel <- function(x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "num_features")
                           ),
                           max_features = 5L,
                           feature_pool = NULL,
                           cohort_aggregator = "pairwise_ratios",
                           feature_alignment = "intersection",
                           constraints = list(),
                           scoring_fn = NULL,
                           nsga_control = list(),
                           assay = NULL,
                           seed = NULL,
                           fitness_cv = TRUE,
                           fitness_cv_folds = 5L,
                           final_model_cv = FALSE,
                           cv_folds = 5L,
                           regularized = TRUE,
                           regularized_alpha = 0.5) {
  feature_alignment <- match.arg(feature_alignment)
  if (!requireNamespace("rmoo", quietly = TRUE)) {
    stop("The 'rmoo' package is required. Install it via install.packages('rmoo').",
         call. = FALSE)
  }

  if (!is.character(cohort_aggregator) || length(cohort_aggregator) != 1L ||
      !nzchar(cohort_aggregator)) {
    stop("`cohort_aggregator` must be a single non-empty character string.", call. = FALSE)
  }
  if (!exists(cohort_aggregator, envir = .aggregator_registry, inherits = FALSE)) {
    available <- ls(.aggregator_registry)
    stop(
      "Unknown aggregator '", cohort_aggregator, "'. ",
      "Available: ", paste(available, collapse = ", "), ". ",
      "Use register_aggregator() to add custom aggregators.",
      call. = FALSE
    )
  }

  inputs_raw <- .prepare_cohort_inputs(x, y, assay = assay, aggregator = "none",
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
        cohort_aggregator != "none" &&
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

  # Warn about memory usage BEFORE aggregation for pairwise aggregators
  if (cohort_aggregator %in% c("pairwise_ratios", "pairwise_log_ratios")) {
    n_base_features <- length(feature_pool_base)
    if (n_base_features > 100) {
      n_pairs <- (as.numeric(n_base_features) * (n_base_features - 1)) / 2
      message(
        sprintf(
          "Note: Pairwise aggregation of %d base features will create %.0f pairs.",
          n_base_features, n_pairs
        )
      )
      if (n_pairs > 10000) {
        warning(
          sprintf(
            "Large feature pool (%d features -> %.0f pairs) may cause slow ",
            n_base_features, n_pairs
          ),
          "optimization. Consider reducing feature_pool via get_top_de_features() ",
          "or select_transferable_features().",
          call. = FALSE
        )
      }
    }
  }

  inputs <- .prepare_cohort_inputs(
    x,
    y,
    assay = assay,
    aggregator = cohort_aggregator,
    feature_subset = feature_pool_base,
    feature_alignment = feature_alignment
  )
  x_mat <- inputs$x
  truth <- inputs$truth
  cohort <- inputs$cohort

  feature_names <- colnames(x_mat)
  if (is.null(feature_names)) {
    feature_names <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    colnames(x_mat) <- feature_names
  }

  if (is.null(feature_pool_final)) {
    feature_pool <- feature_names
  } else {
    feature_pool <- .resolve_feature_pool(feature_pool_final, feature_names)
  }

  if (!length(feature_pool)) {
    stop("`feature_pool` produced zero features.", call. = FALSE)
  }

  min_features_for_regularized <- 2L
  if (regularized && max_features < min_features_for_regularized) {
    stop("`max_features` must be at least ", min_features_for_regularized,
         " when `regularized = TRUE` (glmnet requirement). ",
         "Either increase `max_features` or set `regularized = FALSE`.",
         call. = FALSE)
  }

  if (max_features < 1L) {
    stop("`max_features` must be at least 1.", call. = FALSE)
  }

  if (max_features > length(feature_pool)) {
    max_features <- length(feature_pool)
  }

  x_pool <- x_mat[, feature_pool, drop = FALSE]
  decision_dim <- ncol(x_pool)

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

  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim)
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

  # Minimum features needed: glmnet requires >= 2 features when regularized

  min_features_required <- if (regularized) 2L else 1L

  evaluate_candidate <- function(decision_vec) {
    # Threshold-based selection: include features with weight > 0.5
    # This allows variable panel sizes, making num_features objective meaningful
    feature_names_pool <- colnames(x_pool)
    above_threshold <- which(decision_vec > 0.5)

    if (length(above_threshold) < min_features_required) {
      # Fallback: take top min_features_required by weight
      ord <- order(decision_vec, feature_names_pool, decreasing = c(TRUE, FALSE),
                   method = "radix")
      selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
    } else if (length(above_threshold) > max_features) {
      # Cap at max_features, using feature names for reproducible tie-breaking
      weights_above <- decision_vec[above_threshold]
      names_above <- feature_names_pool[above_threshold]
      ord <- order(weights_above, names_above, decreasing = c(TRUE, FALSE),
                   method = "radix")
      selected_idx <- above_threshold[ord[seq_len(max_features)]]
    } else {
      selected_idx <- above_threshold
    }

    # Sort by descending weight for consistent output ordering
    selected_idx <- selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
    selected_features <- feature_names_pool[selected_idx]
    x_selected <- x_pool[, selected_features, drop = FALSE]

    # Compute scores using CV (out-of-fold predictions) or in-sample
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
    metrics <- vapply(objectives, function(obj) {
      obj$fun(
        truth,
        scores,
        selected = selected_features,
        cohort = cohort,
        x = x_selected
      )
    }, numeric(1))
    list(
      features = selected_features,
      scores = scores,
      metrics = metrics,
      constraint_results = constraint_results,
      feasible = feasible
    )
  }

  # Evaluate a single decision vector and return converted objective values
  evaluate_single <- function(decision_vec) {
    evaluated <- evaluate_candidate(decision_vec)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(Inf, length(objectives)))
    }
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      # Handle NA, Inf, and -Inf explicitly
      if (is.na(val)) {
        return(Inf)  # Treat NA as infeasible
      }
      if (is.infinite(val)) {
        # For infinite values, always return +Inf to indicate infeasibility
        # This prevents issues with NSGA-II when mixing +Inf and -Inf
        return(Inf)
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
  objective_wrapper <- function(x) {
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

  # Set random seed for reproducibility if provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("`seed` must be a single integer value.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  nsga_result <- do.call(rmoo::nsga2, nsga_params)

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

  select_primary <- function(idx_vec, direction) {
    if (direction == "maximize") {
      which.max(idx_vec)
    } else {
      which.min(idx_vec)
    }
  }

  primary_obj <- names(objectives)[1]
  primary_dir <- objective_directions[[primary_obj]]
  primary_idx <- select_primary(metric_matrix[, primary_obj], primary_dir)

  primary_solution <- solutions[[primary_idx]]
  primary_metrics <- primary_solution$metrics

  objective_df <- do.call(rbind, lapply(seq_along(solutions), function(i) {
    data.frame(
      solution_id = i,
      objective = names(objectives),
      value = solutions[[i]]$metrics,
      direction = objective_directions,
      features = I(rep(list(solutions[[i]]$features), length(objectives))),
      constraints = I(rep(list(solutions[[i]]$constraint_results), length(objectives))),
      stringsAsFactors = FALSE
    )
  }))

  # Fit the final model on the selected features for storage
  selected_features <- primary_solution$features
  x_selected <- x_pool[, selected_features, drop = FALSE]
  if (regularized) {
    final_model <- .fit_final_model_regularized(x_selected, truth, cohort,
                                                 alpha = regularized_alpha)
  } else if (final_model_cv) {
    final_model <- .fit_final_model_cv(x_selected, truth, cohort, cv_folds)
  } else {
    final_model <- .fit_final_model(x_selected, truth, cohort)
  }

  panel <- new(
    "BiomarkerPanelResult",
    features = primary_solution$features,
    metrics = setNames(as.numeric(primary_metrics), names(objectives)),
    objectives = objective_df,
    control = list(
      max_features = max_features,
      feature_pool = feature_pool,
      base_feature_pool = feature_pool_base,
      nsga2 = nsga_args,
      scoring_function = deparse(substitute(scoring_fn)),
      cohort_aggregator = cohort_aggregator,
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
      final_model_cv = final_model_cv,
      cv_folds = if (final_model_cv) cv_folds else NULL,
      regularized = regularized,
      regularized_alpha = if (regularized) regularized_alpha else NULL
    ),
    training_data = list(
      n = nrow(x_mat),
      p = ncol(x_mat),
      class_balance = table(truth),
      feature_pool_size = length(feature_pool),
      base_feature_pool_size = length(feature_pool_base),
      num_cohorts = length(inputs$cohort_names),
      cohort_labels = inputs$cohort_names,
      cohort_counts = inputs$cohort_counts
    ),
    model = final_model
  )

  panel
}
