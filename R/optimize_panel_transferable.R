#' Optimize Biomarker Panels with Transferability Focus
#'
#' Wrapper implementing a train/validation/held-out split strategy for
#' cross-cohort generalizability. Partitions each cohort into training
#' (for NSGA optimization), validation (for solution selection), and
#' held-out (stored for downstream NP threshold calibration via
#' [calibrate_panel()]).
#'
#' Uses base-feature-first selection: NSGA searches over base features
#' (O(n) search space), then feature transforms (e.g., pairwise ratios) are
#' applied on-the-fly during fitness evaluation.
#'
#' @param x List of matrix-like objects representing cohorts.
#' @param y List of binary response vectors aligned with `x`.
#' @param objectives Named list of objective descriptors as returned by
#'   [define_objectives()].
#' @param max_features Maximum number of base biomarkers permitted in a panel.
#' @param feature_pool Optional subset of feature identifiers. If `NULL`,
#'   `select_transferable_features()` is run on training data only.
#' @param train_ratio Proportion of each cohort for training (default 0.7).
#' @param val_ratio Proportion of each cohort for validation (default 0.2).
#' @param np_alpha Type I error rate for Neyman-Pearson threshold (default 0.15).
#'   Stored in the result for use by [calibrate_panel()].
#' @param np_delta Tolerance for NP threshold selection (default 0.05).
#'   Stored in the result for use by [calibrate_panel()].
#' @param feature_transform Transformation applied to selected base features
#'   on-the-fly during optimization (default `"pairwise_ratios"`). See
#'   [feature_transform_registry()] for options.
#' @param feature_alignment Strategy for aligning features across cohorts
#'   (default `"intersection"`).
#' @param constraints Optional list of constraint descriptors.
#' @param algorithm Multi-objective algorithm: `"NSGA-III"` (default) or
#'   `"NSGA-II"`.
#' @param nsga_control Named list of arguments passed to rmoo.
#' @param assay For `SummarizedExperiment` inputs, assay name or index.
#' @param seed Optional integer seed for reproducibility.
#' @param regularized Logical; if `TRUE` (default), use regularized regression.
#' @param regularized_alpha Elastic net mixing parameter (default 0.5).
#' @param n_top_features Number of top features to select via
#'   `select_transferable_features()` when `feature_pool` is NULL (default 50).
#' @return An `OptimizationResult` containing Pareto-optimal solutions and
#'   held-out data for downstream calibration. Use [summarize_solutions()] to
#'   inspect solutions, [fit_panel()] to fit a model, [calibrate_panel()] for
#'   NP threshold calibration, and [evaluate_panel()] for validation.
#' @export
#' @seealso [fit_panel()], [calibrate_panel()], [evaluate_panel()],
#'   [summarize_solutions()]
optimize_panel_transferable <- function(
  x, y,
  objectives = define_objectives(
    metrics = c("sensitivity", "specificity", "num_features")
  ),
  max_features = 5L,
  feature_pool = NULL,
  train_ratio = 0.7,
  val_ratio = 0.2,
  np_alpha = 0.15,
  np_delta = 0.05,
  feature_transform = "pairwise_ratios",
  feature_alignment = "intersection",
  constraints = list(),
  algorithm = c("NSGA-III", "NSGA-II"),
  nsga_control = list(),
  assay = NULL,
  seed = NULL,
  regularized = TRUE,
  regularized_alpha = 0.5,
  n_top_features = 50L
) {
  algorithm <- match.arg(algorithm)

  # Validate inputs
  if (!is.list(x) || length(x) < 1L) {
    stop("`x` must be a list of cohort matrices.", call. = FALSE)
  }
  if (!is.list(y) || length(y) != length(x)) {
    stop("`y` must be a list of response vectors matching `x`.", call. = FALSE)
  }

  # 1. Validate partition ratios
  .validate_partition_ratios(train_ratio, val_ratio)

  # 2. Stratified partition each cohort
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  partitions <- .stratified_partition_cohorts(x, y, train_ratio, val_ratio)

  # 3. Feature selection on TRAINING data only (if feature_pool is NULL)
  if (is.null(feature_pool)) {
    message("Selecting transferable features from training data...")
    fs_result <- select_transferable_features(
      partitions$train_x,
      partitions$train_y,
      n_features = n_top_features
    )
    feature_pool <- fs_result$features
    message("Selected ", length(feature_pool), " features for optimization.")
  }

  # 4. Prepare combined datasets with NO transform (base features only)
  train_inputs <- .prepare_cohort_inputs(
    partitions$train_x, partitions$train_y,
    assay = assay,
    transform = "none",
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )
  val_inputs <- .prepare_cohort_inputs(
    partitions$val_x, partitions$val_y,
    assay = assay,
    transform = "none",
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )
  heldout_inputs <- .prepare_cohort_inputs(
    partitions$heldout_x, partitions$heldout_y,
    assay = assay,
    transform = "none",
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )

  # Base feature pool (O(n) search space, not O(n^2))
  base_feature_pool <- colnames(train_inputs$x)

  # Minimum base features required depends on transform and regularization
  if (regularized && feature_transform %in% c("pairwise_ratios", "pairwise_log_ratios")) {
    min_features_required <- 3L
  } else if (regularized) {
    min_features_required <- 2L
  } else {
    min_features_required <- 1L
  }

  # 5. Create validation-based fitness function (base-feature-first)
  fitness_fn <- .make_validation_fitness(
    train_x = train_inputs$x,
    train_y = train_inputs$truth,
    train_cohort = train_inputs$cohort,
    val_x = val_inputs$x,
    val_y = val_inputs$truth,
    val_cohort = val_inputs$cohort,
    feature_pool = base_feature_pool,
    max_features = max_features,
    objectives = objectives,
    constraints = constraints,
    regularized = regularized,
    alpha = regularized_alpha,
    feature_transform = feature_transform,
    min_features_required = min_features_required
  )

  # 6. Get NSGA defaults and run optimization
  decision_dim <- length(base_feature_pool)
  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim, algorithm)
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)

  nsga_params <- c(
    list(
      type = "real-valued",
      fitness = fitness_fn$wrapper,
      nObj = length(objectives),
      lower = rep(0, decision_dim),
      upper = rep(1, decision_dim),
      monitor = FALSE,
      summary = FALSE
    ),
    nsga_args
  )

  if (algorithm == "NSGA-III" && is.null(nsga_args$n_partitions)) {
    nsga_params$n_partitions <- .compute_nsga3_partitions(length(objectives))
  }

  message("Running ", algorithm, " optimization...")
  if (algorithm == "NSGA-II") {
    nsga_result <- do.call(rmoo::nsga2, nsga_params)
  } else {
    nsga_result <- do.call(rmoo::nsga3, nsga_params)
  }

  # 7. Extract Pareto front and build solutions
  optimal_idx <- which(nsga_result@front == 1)
  pareto_pop <- nsga_result@population[optimal_idx, , drop = FALSE]
  if (is.null(dim(pareto_pop))) {
    pareto_pop <- matrix(pareto_pop, nrow = 1)
  }

  solutions <- lapply(seq_len(nrow(pareto_pop)), function(i) {
    decision_vec <- pareto_pop[i, ]
    fitness_fn$evaluate(decision_vec)
  })

  feasible_vec <- vapply(solutions, `[[`, logical(1), "feasible")
  if (!any(feasible_vec)) {
    stop("No solutions satisfied the supplied constraints.", call. = FALSE)
  }
  solutions <- solutions[feasible_vec]

  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  # Build solutions data frame (all Pareto solutions)
  solutions_df <- data.frame(
    solution_id = seq_along(solutions),
    stringsAsFactors = FALSE
  )
  solutions_df$base_features <- I(lapply(solutions, `[[`, "base_features"))
  solutions_df$features <- I(lapply(solutions, `[[`, "features"))
  for (obj_name in names(objectives)) {
    solutions_df[[obj_name]] <- metric_matrix[, obj_name]
  }

  # 8. Combine train + val raw base features for fit_panel()
  combined_x <- rbind(train_inputs$x, val_inputs$x)
  combined_y <- factor(
    c(as.character(train_inputs$truth), as.character(val_inputs$truth)),
    levels = levels(train_inputs$truth)
  )
  combined_cohort <- factor(
    c(as.character(train_inputs$cohort), as.character(val_inputs$cohort)),
    levels = unique(c(levels(train_inputs$cohort), levels(val_inputs$cohort)))
  )

  # Build control parameters
  control <- list(
    max_features = max_features,
    feature_pool = base_feature_pool,
    algorithm = algorithm,
    nsga2 = nsga_args,
    feature_transform = feature_transform,
    feature_alignment = feature_alignment,
    positive_class = "Yes",
    response_levels = levels(train_inputs$truth),
    seed = seed,
    selection_threshold = 0.5,
    regularized = regularized,
    regularized_alpha = if (regularized) regularized_alpha else NULL,
    objective_directions = objective_directions,
    # Transferable-specific: held-out data for calibrate_panel()
    heldout_x = heldout_inputs$x,
    heldout_y = heldout_inputs$truth,
    heldout_cohort = heldout_inputs$cohort,
    np_alpha = np_alpha,
    np_delta = np_delta,
    partition_info = list(
      train_ratio = train_ratio,
      val_ratio = val_ratio,
      heldout_ratio = 1 - train_ratio - val_ratio,
      partition_sizes = partitions$partition_info
    )
  )

  # Build training signature
  training_signature <- list(
    n = nrow(combined_x),
    p = ncol(combined_x),
    class_balance = table(combined_y),
    feature_pool_size = length(base_feature_pool),
    num_cohorts = length(partitions$cohort_names),
    cohort_labels = partitions$cohort_names
  )

  new("OptimizationResult",
      solutions = solutions_df,
      feature_pool = base_feature_pool,
      control = control,
      training_signature = training_signature,
      aggregated_x = combined_x,
      aggregated_y = combined_y,
      aggregated_cohort = combined_cohort)
}


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Create Validation-Based Fitness Function (Base-Feature-First)
#'
#' Factory that returns a fitness function for rmoo that selects base features,
#' applies feature transforms on-the-fly, trains on training data, and
#' evaluates objectives on validation data.
#'
#' @param train_x Training feature matrix (base features).
#' @param train_y Training response factor.
#' @param train_cohort Training cohort factor.
#' @param val_x Validation feature matrix (base features).
#' @param val_y Validation response factor.
#' @param val_cohort Validation cohort factor.
#' @param feature_pool Character vector of base feature names.
#' @param max_features Maximum base features to include.
#' @param objectives Objective list from define_objectives().
#' @param constraints Constraint list.
#' @param regularized Whether to use regularized regression.
#' @param alpha Elastic net mixing parameter.
#' @param feature_transform Name of the feature transform to apply.
#' @param min_features_required Minimum number of base features.
#' @return List with `wrapper` (fitness for rmoo) and `evaluate` (full eval).
#' @keywords internal
.make_validation_fitness <- function(
  train_x, train_y, train_cohort,
  val_x, val_y, val_cohort,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL
) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)

  if (is.null(min_features_required)) {
    min_features_required <- if (regularized) 2L else 1L
  }

  evaluate_candidate <- function(decision_vec) {
    feature_names_pool <- feature_pool
    above_threshold <- which(decision_vec > 0.5)

    if (length(above_threshold) < min_features_required) {
      ord <- order(decision_vec, feature_names_pool,
        decreasing = c(TRUE, FALSE),
        method = "radix"
      )
      selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
    } else if (length(above_threshold) > max_features) {
      weights_above <- decision_vec[above_threshold]
      names_above <- feature_names_pool[above_threshold]
      ord <- order(weights_above, names_above,
        decreasing = c(TRUE, FALSE),
        method = "radix"
      )
      selected_idx <- above_threshold[ord[seq_len(max_features)]]
    } else {
      selected_idx <- above_threshold
    }

    selected_idx <- selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
    selected_base_features <- feature_names_pool[selected_idx]

    # Extract base features from both sets
    x_train_base <- train_x[, selected_base_features, drop = FALSE]
    x_val_base <- val_x[, selected_base_features, drop = FALSE]

    # Apply transform on-the-fly (mirrors optimize_panel.R)
    if (feature_transform != "none" && ncol(x_train_base) >= 2L) {
      x_train_sel <- .apply_feature_transform_single(x_train_base, feature_transform)
      x_val_sel <- .apply_feature_transform_single(x_val_base, feature_transform)
      selected_features <- colnames(x_train_sel)
    } else {
      x_train_sel <- x_train_base
      x_val_sel <- x_val_base
      selected_features <- selected_base_features
    }

    tryCatch(
      {
        if (regularized) {
          fit <- .fit_final_model_regularized(
            x_train_sel, train_y, train_cohort,
            alpha = alpha
          )
        } else {
          fit <- .fit_final_model(x_train_sel, train_y, train_cohort)
        }

        # Predict on validation data
        val_scores <- .predict_from_model(fit, x_val_sel, cohort = val_cohort)

        # Compute constraints on validation data
        constraint_results <- if (length(constraint_specs)) {
          setNames(vapply(seq_along(constraint_specs), function(j) {
            res <- constraint_specs[[j]]$fun(
              truth = val_y,
              scores = val_scores,
              selected = selected_features,
              cohort = val_cohort,
              x = x_val_sel
            )
            isTRUE(res)
          }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
        } else {
          logical(0)
        }

        feasible <- if (length(constraint_results)) all(constraint_results) else TRUE

        metrics <- vapply(objectives, function(obj) {
          obj$fun(
            val_y,
            val_scores,
            selected = selected_features,
            cohort = val_cohort,
            x = x_val_sel
          )
        }, numeric(1))

        list(
          base_features = selected_base_features,
          features = selected_features,
          scores = val_scores,
          metrics = metrics,
          constraint_results = constraint_results,
          feasible = feasible
        )
      },
      error = function(e) {
        list(
          base_features = selected_base_features,
          features = selected_features,
          scores = rep(0.5, nrow(val_x)),
          metrics = setNames(rep(Inf, length(objectives)), names(objectives)),
          constraint_results = logical(0),
          feasible = length(constraint_specs) == 0L
        )
      }
    )
  }

  # Large finite penalty instead of Inf -- NSGA-III normalization produces NaN
  # from Inf values, causing "missing value where TRUE/FALSE needed" errors
  .PENALTY <- 1e6

  evaluate_single <- function(decision_vec) {
    evaluated <- evaluate_candidate(decision_vec)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(.PENALTY, length(objectives)))
    }
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      if (is.na(val) || is.infinite(val)) {
        return(.PENALTY)
      }
      if (dir == "maximize") {
        return(-val)
      }
      val
    }, val = metrics, dir = objective_directions, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    as.numeric(converted)
  }

  objective_wrapper <- function(x, ...) {
    if (is.null(dim(x))) {
      return(evaluate_single(x))
    }
    t(apply(x, 1, evaluate_single))
  }

  list(
    wrapper = objective_wrapper,
    evaluate = evaluate_candidate
  )
}
