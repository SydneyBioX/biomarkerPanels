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
#' @param selection_threshold Either a fixed numeric threshold in (0,1) for
#'   feature selection, or `"adaptive"` (default) to use gap-based selection
#'   that allows variable panel sizes in the Pareto front. See
#'   [optimize_panel()] for details.
#' @param n_top_features Number of top features to select via
#'   `select_transferable_features()` when `feature_pool` is NULL (default 50).
#' @param fitness_mode Fitness-evaluation strategy:
#'   `"within_cohort_val"` (default) trains on the combined train partition and
#'   evaluates on the combined validation partition — both partitions contain
#'   samples from every cohort, so a candidate can win by fitting
#'   within-cohort structure. `"within_cohort_rotating"` pools train and
#'   validation samples, pre-computes `n_val_splits` stratified train/val
#'   splits, and rotates to the next split every NSGA generation. This is
#'   the MOO analogue of mini-batch SGD: no single split can be overfit to,
#'   and final Pareto-solution metrics are averaged across all K splits for
#'   an honest cross-split number. Feature selection still runs on the
#'   original training partition only, so held-out labels never influence
#'   the feature pool. `"loco"` (leave-one-cohort-out) ignores the
#'   validation partition and pools train+val per cohort; each candidate is
#'   scored by iteratively training on all-but-one cohort and predicting on the
#'   held-out cohort. Metrics are computed on the concatenated out-of-cohort
#'   predictions, making cross-cohort transferability the objective NSGA
#'   actually optimizes for. Requires `>=2` cohorts.
#' @param n_val_splits Number of rotating train/val splits to pre-compute
#'   when `fitness_mode = "within_cohort_rotating"` (default 10).
#' @param cache_fitness Logical; if `TRUE` (default), cache candidate fitness
#'   by selected base-feature panel within each validation context. Set
#'   `FALSE` for intentionally stochastic custom objectives.
#' @param cache_max_entries Maximum number of selected-panel entries retained
#'   per fitness cache. Defaults to `Inf`.
#' @param record_history Logical; if `TRUE`, capture the population, fitness,
#'   and NSGA rank at every generation and attach the result to the returned
#'   `OptimizationResult`. Retrieve via [nsga_history()]. Default `FALSE`.
#'   Adds modest overhead (one matrix copy per generation) and is intended for
#'   diagnostic / visualization use, e.g. plotting how the Pareto front
#'   evolves over `maxiter`.
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
  selection_threshold = "adaptive",
  n_top_features = 50L,
  fitness_mode = c("within_cohort_val", "within_cohort_rotating", "loco"),
  n_val_splits = 10L,
  cache_fitness = TRUE,
  cache_max_entries = Inf,
  record_history = FALSE
) {
  algorithm <- match.arg(algorithm)
  fitness_mode <- match.arg(fitness_mode)

  if (fitness_mode == "within_cohort_rotating") {
    if (!is.numeric(n_val_splits) || length(n_val_splits) != 1L ||
        is.na(n_val_splits) || n_val_splits < 2L) {
      stop("`n_val_splits` must be an integer >= 2.", call. = FALSE)
    }
    n_val_splits <- as.integer(n_val_splits)
  }

  # Validate inputs
  if (!is.list(x) || length(x) < 1L) {
    stop("`x` must be a list of cohort matrices.", call. = FALSE)
  }
  if (!is.list(y) || length(y) != length(x)) {
    stop("`y` must be a list of response vectors matching `x`.", call. = FALSE)
  }

  # Validate numeric parameters
  if (!is.numeric(regularized_alpha) || length(regularized_alpha) != 1L ||
      is.na(regularized_alpha) || regularized_alpha < 0 || regularized_alpha > 1) {
    stop("`regularized_alpha` must be a single numeric value in [0, 1].", call. = FALSE)
  }
  if (!identical(selection_threshold, "adaptive")) {
    st <- suppressWarnings(as.numeric(selection_threshold))
    if (is.na(st) || st <= 0 || st >= 1) {
      stop(
        "`selection_threshold` must be \"adaptive\" or a numeric value in (0, 1).",
        call. = FALSE
      )
    }
  }
  .validate_cache_controls(cache_fitness, cache_max_entries)

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

  # 5. Create fitness function (base-feature-first)
  if (fitness_mode == "within_cohort_rotating") {
    # Pool train + val; pre-compute K stratified splits and rotate per
    # generation so no single split can be overfit to. Feature pool was
    # already selected on the ORIGINAL train partition only, so held-out
    # labels never reach feature selection.
    pool_x <- rbind(train_inputs$x, val_inputs$x)
    pool_y <- factor(
      c(as.character(train_inputs$truth), as.character(val_inputs$truth)),
      levels = levels(train_inputs$truth)
    )
    pool_cohort <- factor(
      c(as.character(train_inputs$cohort), as.character(val_inputs$cohort)),
      levels = unique(c(levels(train_inputs$cohort), levels(val_inputs$cohort)))
    )
    split_val_ratio <- val_ratio / (train_ratio + val_ratio)
    rotating_splits <- .generate_stratified_splits(
      y = pool_y,
      cohort = pool_cohort,
      n_splits = n_val_splits,
      val_ratio = split_val_ratio,
      seed = seed
    )
    fitness_fn <- .make_rotating_validation_fitness(
      pool_x = pool_x,
      pool_y = pool_y,
      pool_cohort = pool_cohort,
      splits = rotating_splits,
      feature_pool = base_feature_pool,
      max_features = max_features,
      objectives = objectives,
      constraints = constraints,
      regularized = regularized,
      alpha = regularized_alpha,
      feature_transform = feature_transform,
      min_features_required = min_features_required,
      selection_threshold = selection_threshold,
      cache_fitness = cache_fitness,
      cache_max_entries = cache_max_entries
    )
  } else if (fitness_mode == "loco") {
    # Pool train+val per cohort; LOCO loop provides the held-out evaluation.
    loco_x <- rbind(train_inputs$x, val_inputs$x)
    loco_y <- factor(
      c(as.character(train_inputs$truth), as.character(val_inputs$truth)),
      levels = levels(train_inputs$truth)
    )
    loco_cohort <- factor(
      c(as.character(train_inputs$cohort), as.character(val_inputs$cohort)),
      levels = unique(c(levels(train_inputs$cohort), levels(val_inputs$cohort)))
    )
    if (nlevels(droplevels(loco_cohort)) < 2L) {
      stop(
        "fitness_mode = 'loco' requires >=2 cohorts but only ",
        nlevels(droplevels(loco_cohort)),
        " was available after partitioning.",
        call. = FALSE
      )
    }
    fitness_fn <- .make_loco_fitness(
      x = loco_x,
      y = loco_y,
      cohort = loco_cohort,
      feature_pool = base_feature_pool,
      max_features = max_features,
      objectives = objectives,
      constraints = constraints,
      regularized = regularized,
      alpha = regularized_alpha,
      feature_transform = feature_transform,
      min_features_required = min_features_required,
      selection_threshold = selection_threshold,
      cache_fitness = cache_fitness,
      cache_max_entries = cache_max_entries
    )
  } else {
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
      min_features_required = min_features_required,
      selection_threshold = selection_threshold,
      cache_fitness = cache_fitness,
      cache_max_entries = cache_max_entries
    )
  }

  # 6. Get NSGA defaults and run optimization
  decision_dim <- length(base_feature_pool)
  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim, algorithm)
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)

  # Hoisted: needed both by the per-generation history monitor (to convert
  # minimised fitness back to user-facing direction) and by the post-run
  # dominance filter further down.
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")

  # Same shared selector used by the fitness factories, materialized here for
  # per-generation history.
  panel_selector <- .make_panel_selector(
    feature_pool = base_feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )
  select_base_features <- function(decision_vec) {
    panel_selector(decision_vec)$base_features
  }

  # Optional per-generation history capture (see optimize_panel() for rationale).
  history_buffer <- new.env(parent = emptyenv())
  history_buffer$gens <- list()
  history_monitor <- function(object, ...) {
    iter <- object@iter
    fit <- object@fitness
    pop <- object@population
    if (is.null(dim(fit))) fit <- matrix(fit, nrow = 1L)
    if (is.null(dim(pop))) pop <- matrix(pop, nrow = 1L)
    fr <- as.integer(object@front)
    for (j in seq_along(objective_directions)) {
      if (objective_directions[j] == "maximize") fit[, j] <- -fit[, j]
    }
    colnames(fit) <- names(objectives)
    history_buffer$gens[[length(history_buffer$gens) + 1L]] <- list(
      iter = as.integer(iter),
      fit = fit,
      front = fr,
      pop = pop
    )
    invisible(NULL)
  }

  monitor_arg <- if (isTRUE(record_history)) history_monitor else FALSE

  nsga_params <- c(
    list(
      type = "real-valued",
      fitness = fitness_fn$wrapper,
      nObj = length(objectives),
      lower = rep(0, decision_dim),
      upper = rep(1, decision_dim),
      monitor = monitor_arg,
      summary = FALSE
    ),
    nsga_args
  )

  if (algorithm == "NSGA-III" && is.null(nsga_args$n_partitions)) {
    nsga_params$n_partitions <- .compute_nsga3_partitions(length(objectives))
  }

  # Seed initial population with diverse panel sizes
  n_suggestions <- min(20L, nsga_args$popSize %/% 4L)
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

  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  # Post-filter dominated solutions (see .filter_dominated docs)
  nondom_idx <- .filter_dominated(metric_matrix, objective_directions)
  solutions <- solutions[nondom_idx]
  metric_matrix <- metric_matrix[nondom_idx, , drop = FALSE]

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
    selection_threshold = selection_threshold,
    cache_fitness = cache_fitness,
    cache_max_entries = cache_max_entries,
    regularized = regularized,
    regularized_alpha = if (regularized) regularized_alpha else NULL,
    fitness_mode = fitness_mode,
    n_val_splits = if (fitness_mode == "within_cohort_rotating") n_val_splits else NULL,
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

  # Materialize per-generation history if it was captured.
  history_out <- if (isTRUE(record_history) && length(history_buffer$gens)) {
    gen_dfs <- lapply(history_buffer$gens, function(g) {
      n_ind <- nrow(g$pop)
      bf <- lapply(seq_len(n_ind), function(i) select_base_features(g$pop[i, ]))
      df <- data.frame(
        generation = g$iter,
        individual = seq_len(n_ind),
        rank = g$front,
        is_pareto = g$front == 1L,
        g$fit,
        n_features = lengths(bf),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      df$base_features <- bf
      df
    })
    do.call(rbind, gen_dfs)
  } else {
    list()
  }

  new("OptimizationResult",
      solutions = solutions_df,
      feature_pool = base_feature_pool,
      control = control,
      training_signature = training_signature,
      aggregated_x = combined_x,
      aggregated_y = combined_y,
      aggregated_cohort = combined_cohort,
      history = history_out)
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
#' @param selection_threshold Either `"adaptive"` or a fixed numeric threshold.
#' @param cache_fitness Logical; cache duplicate selected panels during fitness.
#' @param cache_max_entries Maximum entries retained in the fitness cache.
#' @return List with `wrapper` (fitness for rmoo) and `evaluate` (full eval).
#' @keywords internal
.make_validation_fitness <- function(
  train_x, train_y, train_cohort,
  val_x, val_y, val_cohort,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL,
  selection_threshold = "adaptive",
  cache_fitness = TRUE,
  cache_max_entries = Inf
) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)

  if (is.null(min_features_required)) {
    min_features_required <- if (regularized) 2L else 1L
  }

  panel_selector <- .make_panel_selector(
    feature_pool = feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )
  transform_panel <- .make_panel_transformer(
    matrices = list(train = train_x, val = val_x),
    feature_transform = feature_transform,
    cache_max_entries = cache_max_entries
  )
  objective_cache <- .new_fitness_cache(cache_max_entries)
  glm_design_terms <- if (!regularized) {
    .prepare_glm_design_terms(
      cohort_train = train_cohort,
      cohort_new = val_cohort,
      predict_cohort = "reference"
    )
  } else {
    NULL
  }

  evaluate_candidate <- function(decision_vec = NULL, selection = NULL) {
    if (is.null(selection)) {
      selection <- panel_selector(decision_vec)
    }
    selected_base_features <- selection$base_features

    transformed <- transform_panel(selected_base_features)
    x_train_sel <- transformed$matrices$train
    x_val_sel <- transformed$matrices$val
    selected_features <- transformed$features

    tryCatch(
      {
        if (regularized) {
          fit <- .fit_final_model_regularized(
            x_train_sel, train_y, train_cohort,
            alpha = alpha
          )
          # Predict on validation data
          val_scores <- .predict_from_model(fit, x_val_sel, cohort = val_cohort)
        } else {
          val_scores <- .fit_predict_binomial_glm(
            x_train = x_train_sel,
            truth = train_y,
            x_new = x_val_sel,
            cohort_train = train_cohort,
            cohort_new = val_cohort,
            predict_cohort = "reference",
            design_terms = glm_design_terms
          )
        }

        # Compute constraints on validation data
        constraint_results <- if (length(constraint_specs)) {
          setNames(vapply(seq_along(constraint_specs), function(j) {
            res <- constraint_specs[[j]]$fun(
              truth = val_y,
              scores = val_scores,
              selected = selected_base_features,
              cohort = val_cohort,
              x = x_val_sel
            )
            isTRUE(res)
          }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
        } else {
          logical(0)
        }

        feasible <- if (length(constraint_results)) all(constraint_results) else TRUE

        # Pass base features as `selected` so metric_num_features counts genes,
        # not pairwise ratios. Other metrics ignore `selected`.
        metrics <- vapply(objectives, function(obj) {
          obj$fun(
            val_y,
            val_scores,
            selected = selected_base_features,
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
        stop(
          "Candidate evaluation failed for base features [",
          paste(selected_base_features, collapse = ", "), "]: ",
          conditionMessage(e),
          "\nThis error occurred during fitness evaluation in the NSGA loop. ",
          "Common causes: feature alignment mismatch between training and ",
          "validation data, perfect separation in logistic regression, or ",
          "zero-variance features after transform. ",
          "Consider using regularized = TRUE or reducing max_features.",
          call. = FALSE
        )
      }
    )
  }

  # Large finite penalty instead of Inf -- NSGA-III normalization produces NaN
  # from Inf values, causing "missing value where TRUE/FALSE needed" errors
  .PENALTY <- 1e6

  evaluate_single <- function(decision_vec = NULL, selection = NULL) {
    evaluated <- evaluate_candidate(decision_vec = decision_vec,
                                    selection = selection)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(.PENALTY, length(objectives)))
    }
    .convert_metrics_to_objectives(evaluated$metrics, objective_directions,
                                   penalty = .PENALTY)
  }

  objective_wrapper <- function(x, ...) {
    .evaluate_fitness_population(
      x = x,
      selector = panel_selector,
      evaluate_selection = function(selection) {
        evaluate_single(selection = selection)
      },
      n_objectives = length(objectives),
      cache = objective_cache,
      cache_fitness = cache_fitness
    )
  }

  list(
    wrapper = objective_wrapper,
    evaluate = evaluate_candidate
  )
}


# -----------------------------------------------------------------------------
# Leave-One-Cohort-Out Fitness
# -----------------------------------------------------------------------------

#' Create Leave-One-Cohort-Out Fitness Function
#'
#' Factory that scores each candidate by iterating over cohorts: train on all
#' cohorts except one, predict on the held-out cohort, concatenate
#' out-of-cohort predictions, and compute objective metrics on the pooled
#' OOF scores. This forces the optimizer to select features that transfer
#' across cohort boundaries rather than features that simply fit within-cohort
#' validation samples.
#'
#' @param x Pooled base-feature matrix (samples x features).
#' @param y Pooled response factor.
#' @param cohort Factor aligned with rows of `x`; must have >=2 levels.
#' @param feature_pool Character vector of base feature names (cols of `x`).
#' @param max_features Maximum base features per candidate.
#' @param objectives Objective list from [define_objectives()].
#' @param constraints Constraint list.
#' @param regularized Whether the per-fold model is elastic-net.
#' @param alpha Elastic net mixing parameter.
#' @param feature_transform Name of the feature transform to apply on-the-fly.
#' @param min_features_required Minimum base features per candidate.
#' @param selection_threshold `"adaptive"` or a numeric gate in (0, 1).
#' @param cache_fitness Logical; cache duplicate selected panels during fitness.
#' @param cache_max_entries Maximum entries retained in the fitness cache.
#' @return List with `wrapper` (vectorized NSGA fitness) and `evaluate`
#'   (single-candidate evaluator returning full diagnostics).
#' @keywords internal
.make_loco_fitness <- function(
  x, y, cohort,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL,
  selection_threshold = "adaptive",
  cache_fitness = TRUE,
  cache_max_entries = Inf
) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)

  if (is.null(min_features_required)) {
    min_features_required <- if (regularized) 2L else 1L
  }

  cohort <- droplevels(factor(cohort))
  cohort_levels <- levels(cohort)
  if (length(cohort_levels) < 2L) {
    stop("LOCO fitness requires at least 2 cohorts.", call. = FALSE)
  }

  # Precompute row indices per cohort for speed in the fitness loop
  fold_rows <- lapply(cohort_levels, function(lv) which(cohort == lv))
  names(fold_rows) <- cohort_levels

  panel_selector <- .make_panel_selector(
    feature_pool = feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )
  transform_panel <- .make_panel_transformer(
    matrices = list(x = x),
    feature_transform = feature_transform,
    cache_max_entries = cache_max_entries
  )
  objective_cache <- .new_fitness_cache(cache_max_entries)
  fold_glm_design_terms <- if (!regularized) {
    lapply(cohort_levels, function(lv) {
      test_idx <- fold_rows[[lv]]
      train_idx <- setdiff(seq_len(nrow(x)), test_idx)
      .prepare_glm_design_terms(
        cohort_train = cohort[train_idx],
        cohort_new = cohort[test_idx],
        predict_cohort = "reference"
      )
    })
  } else {
    NULL
  }
  if (!is.null(fold_glm_design_terms)) {
    names(fold_glm_design_terms) <- cohort_levels
  }

  evaluate_candidate <- function(decision_vec = NULL, selection = NULL) {
    if (is.null(selection)) {
      selection <- panel_selector(decision_vec)
    }
    selected_base_features <- selection$base_features

    transformed <- transform_panel(selected_base_features)
    x_transformed <- transformed$matrices$x
    selected_features <- transformed$features

    oof_scores <- rep(NA_real_, nrow(x_transformed))

    tryCatch(
      {
        for (lv in cohort_levels) {
          test_idx <- fold_rows[[lv]]
          train_idx <- setdiff(seq_len(nrow(x_transformed)), test_idx)

          train_y <- y[train_idx]
          # Skip fold if training data collapses to a single class
          if (length(unique(train_y)) < 2L) next

          x_tr <- x_transformed[train_idx, , drop = FALSE]
          x_te <- x_transformed[test_idx, , drop = FALSE]
          coh_tr <- droplevels(cohort[train_idx])

          if (regularized) {
            fit <- .fit_final_model_regularized(x_tr, train_y, coh_tr, alpha = alpha)
            preds <- .predict_from_model(fit, x_te, cohort = cohort[test_idx])
          } else {
            preds <- .fit_predict_binomial_glm(
              x_train = x_tr,
              truth = train_y,
              x_new = x_te,
              cohort_train = coh_tr,
              cohort_new = cohort[test_idx],
              predict_cohort = "reference",
              design_terms = fold_glm_design_terms[[lv]]
            )
          }
          oof_scores[test_idx] <- as.numeric(preds)
        }

        if (!any(is.finite(oof_scores))) {
          stop("All LOCO folds failed to produce predictions.", call. = FALSE)
        }

        complete <- is.finite(oof_scores)
        eval_truth <- y[complete]
        eval_scores <- oof_scores[complete]
        eval_cohort <- cohort[complete]

        constraint_results <- if (length(constraint_specs)) {
          setNames(vapply(seq_along(constraint_specs), function(j) {
            res <- constraint_specs[[j]]$fun(
              truth = eval_truth,
              scores = eval_scores,
              selected = selected_base_features,
              cohort = eval_cohort,
              x = NULL
            )
            isTRUE(res)
          }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
        } else {
          logical(0)
        }

        feasible <- if (length(constraint_results)) all(constraint_results) else TRUE

        metrics <- vapply(objectives, function(obj) {
          obj$fun(
            eval_truth,
            eval_scores,
            selected = selected_base_features,
            cohort = eval_cohort,
            x = NULL
          )
        }, numeric(1))

        list(
          base_features = selected_base_features,
          features = selected_features,
          scores = oof_scores,
          metrics = metrics,
          constraint_results = constraint_results,
          feasible = feasible
        )
      },
      error = function(e) {
        stop(
          "LOCO candidate evaluation failed for base features [",
          paste(selected_base_features, collapse = ", "), "]: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  .PENALTY <- 1e6

  evaluate_single <- function(decision_vec = NULL, selection = NULL) {
    evaluated <- evaluate_candidate(decision_vec = decision_vec,
                                    selection = selection)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(.PENALTY, length(objectives)))
    }
    .convert_metrics_to_objectives(evaluated$metrics, objective_directions,
                                   penalty = .PENALTY)
  }

  objective_wrapper <- function(x, ...) {
    .evaluate_fitness_population(
      x = x,
      selector = panel_selector,
      evaluate_selection = function(selection) {
        evaluate_single(selection = selection)
      },
      n_objectives = length(objectives),
      cache = objective_cache,
      cache_fitness = cache_fitness
    )
  }

  list(
    wrapper = objective_wrapper,
    evaluate = evaluate_candidate
  )
}
