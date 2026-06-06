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
  # Shared machinery lives in nsga_history_utils.R.
  history_buffer <- .make_history_buffer()
  history_monitor <- .make_history_monitor(
    history_buffer, objective_directions, names(objectives)
  )

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
  history_out <- if (isTRUE(record_history)) {
    .materialize_history(history_buffer, select_base_features)
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
