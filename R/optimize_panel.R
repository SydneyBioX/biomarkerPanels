#' Optimize Biomarker Panels with NSGA-III / NSGA-II
#'
#' Wrapper around [rmoo::nsga3()] (or [rmoo::nsga2()]) that composes registered
#' metric functions into a multi-objective search for compact biomarker panels.
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
#' @param algorithm Multi-objective optimization algorithm. `"NSGA-III"`
#'   (default) uses reference-point-based selection and gives better Pareto
#'   front diversity for many-objective problems (3+ objectives). `"NSGA-II"`
#'   uses crowding distance for diversity preservation and is typically
#'   adequate for 2 objectives.
#' @param nsga_control Named list of arguments passed to [rmoo::nsga3()] (or
#'   [rmoo::nsga2()]). For NSGA-III, `n_partitions` controls reference point
#'   density (default computed from number of objectives). Defaults to adaptive
#'   values based on feature pool size. Note: `parallel = TRUE` is not recommended
#'   as it is slower than sequential execution due to overhead.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @param seed Optional integer seed for reproducibility. When provided, the
#'   seed is set once before partitioning, so both the data split and the NSGA
#'   search are reproducible. If `NULL` (default), no seed is set and results
#'   may vary between runs.
#' @param fitness_mode How candidate panels are scored during the NSGA search:
#'   \describe{
#'     \item{"cv"}{(Default) k-fold cross-validation within the training
#'       partition; objectives are computed on out-of-fold predictions.}
#'     \item{"in_sample"}{Train and score on the same samples. Fast but
#'       overfits; kept for small datasets and backward compatibility.}
#'     \item{"within_cohort_val"}{Train on the training partition, score on the
#'       validation partition. Requires `val_ratio > 0`.}
#'     \item{"within_cohort_rotating"}{Pool train and validation, pre-compute
#'       `n_val_splits` stratified splits, and rotate to the next split every
#'       generation so no single split can be overfit to. Requires
#'       `val_ratio > 0`.}
#'     \item{"loco"}{Leave-one-cohort-out: score each candidate on concatenated
#'       out-of-cohort predictions, making cross-cohort transferability the
#'       objective. Requires at least 2 cohorts.}
#'   }
#' @param fitness_cv_folds Number of cross-validation folds for fitness
#'   evaluation when `fitness_mode = "cv"`. Default is 5.
#' @param n_val_splits Number of rotating train/validation splits to
#'   pre-compute when `fitness_mode = "within_cohort_rotating"` (default 10).
#' @param train_ratio,val_ratio Per-cohort shares used for stratified
#'   partitioning; the remainder (`1 - train_ratio - val_ratio`) is held out and
#'   stored on the result for [calibrate_panel()]. The default
#'   `train_ratio = 1`, `val_ratio = 0` means no partitioning at all: the whole
#'   input is the training data.
#' @param np_alpha Type I error rate for Neyman-Pearson threshold selection.
#'   Stored on the result for [calibrate_panel()] when a held-out partition
#'   exists (default 0.15).
#' @param np_delta Tolerance for Neyman-Pearson threshold selection, stored
#'   alongside `np_alpha` (default 0.05).
#' @param regularized Logical; if `TRUE` (default), use regularized regression
#'   (elastic net via glmnet) for scoring candidates during optimization. This
#'   typically produces more stable fitness estimates and better out-of-sample
#'   performance compared to unregularized logistic regression.
#' @param regularized_alpha Elastic net mixing parameter when `regularized = TRUE`.
#'   `alpha = 1` is lasso, `alpha = 0` is ridge, and values in between give
#'   elastic net. Default is 0.5 (elastic net).
#' @param record_history Logical; if `TRUE`, capture the population, fitness,
#'   and NSGA rank at every generation and attach the result to the returned
#'   `OptimizationResult`. Retrieve via [nsga_history()]. Default `FALSE`.
#'   Adds modest overhead (one matrix copy per generation) and is intended for
#'   diagnostic / visualization use, e.g. plotting how the Pareto front
#'   evolves over `maxiter`.
#' @param selection_threshold Either a fixed numeric threshold in (0,1) for
#'   feature selection, or `"adaptive"` (default) to use gap-based selection
#'   that allows variable panel sizes in the Pareto front. Fixed thresholds
#'   (e.g., 0.5) tend to produce panels at `max_features` because more features

#'   improve accuracy. Adaptive selection finds natural breakpoints in the
#'   weight distribution, enabling the `num_features` objective to drive panel
#'   size diversity.
#' @param fitness_cv Deprecated. Logical alias for `fitness_mode`: `TRUE` maps
#'   to `"cv"` and `FALSE` to `"in_sample"`. Supplying it emits a deprecation
#'   warning and overrides `fitness_mode`.
#' @section Partitioning:
#' When `train_ratio + val_ratio < 1`, each cohort is stratified into
#' train / validation / held-out shares before the search runs; `fitness_mode`
#' determines which of these partitions
#' the NSGA search actually scores against (see above). The returned
#' `OptimizationResult` always carries the pooled train + validation data as
#' `aggregated_x`/`aggregated_y`/`aggregated_cohort`. The held-out partition,
#' `np_alpha`, `np_delta`, and `partition_info` are stored on `control` only
#' when a positive held-out share exists (`train_ratio + val_ratio < 1`); with
#' the default `train_ratio = 1, val_ratio = 0` no partitioning happens at all
#' and those fields are absent. [`calibrate_panel()`] and
#' [`evaluate_pareto_solutions()`] read `control$heldout_x` /
#' `control$heldout_y` / `control$heldout_cohort` when present.
#' @return An `OptimizationResult` containing the Pareto-optimal solutions.
#'   Use [summarize_solutions()] to inspect solutions and [fit_panel()] to
#'   fit a model on a selected solution.
#' @export
#' @seealso [fit_panel()], [summarize_solutions()], [evaluate_panel()]
optimize_panel <- function(x, y,
                           objectives = define_objectives(
                             metrics = c("sensitivity", "specificity", "num_features")
                           ),
                           max_features = 5L,
                           feature_pool = NULL,
                           feature_transform = "pairwise_ratios",
                           feature_alignment = c("intersection", "majority", "impute_median"),
                           constraints = list(),
                           algorithm = c("NSGA-III", "NSGA-II"),
                           nsga_control = list(),
                           assay = NULL,
                           seed = NULL,
                           fitness_mode = c("cv", "in_sample", "within_cohort_val",
                                            "within_cohort_rotating", "loco"),
                           fitness_cv_folds = 5L,
                           n_val_splits = 10L,
                           train_ratio = 1,
                           val_ratio = 0,
                           np_alpha = 0.15,
                           np_delta = 0.05,
                           regularized = TRUE,
                           regularized_alpha = 0.5,
                           selection_threshold = "adaptive",
                           record_history = FALSE,
                           fitness_cv = NULL) {
  algorithm <- match.arg(algorithm)
  feature_alignment <- match.arg(feature_alignment)
  if (!requireNamespace("rmoo", quietly = TRUE)) {
    stop("The 'rmoo' package is required. Install it via install.packages('rmoo').",
         call. = FALSE)
  }

  opts <- .validate_optimize_panel_args(
    fitness_mode = fitness_mode,
    fitness_cv = fitness_cv,
    fitness_cv_folds = fitness_cv_folds,
    n_val_splits = n_val_splits,
    train_ratio = train_ratio,
    val_ratio = val_ratio,
    regularized_alpha = regularized_alpha,
    selection_threshold = selection_threshold,
    feature_transform = feature_transform,
    np_alpha = np_alpha,
    np_delta = np_delta,
    seed = seed
  )
  fitness_mode <- opts$fitness_mode
  fitness_cv_folds <- opts$fitness_cv_folds
  n_val_splits <- opts$n_val_splits
  heldout_ratio <- opts$heldout_ratio

  if (fitness_mode == "in_sample") {
    warning("fitness_mode = \"in_sample\": NSGA will use in-sample scoring, which ",
            "risks overfitting during optimization. Consider fitness_mode = \"cv\" ",
            "for more reliable panel selection.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Resolve the base feature pool against the raw, aligned input. This runs on
  # the whole dataset (feature *names* only, no labels), before any split.
  # ---------------------------------------------------------------------------
  inputs_raw <- .prepare_cohort_inputs(x, y, assay = assay, transform = "none",
                                       feature_alignment = feature_alignment)
  raw_x_mat <- .ensure_feature_colnames(inputs_raw$x)
  raw_feature_names <- colnames(raw_x_mat)

  feature_pool_arg <- feature_pool
  feature_pool_base <- raw_feature_names

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

  # ---------------------------------------------------------------------------
  # Seed once, immediately before partitioning, so the split and the search are
  # both reproducible. Nothing above draws from the RNG; nothing may be
  # inserted between here and .stratified_partition_cohorts() without changing
  # the split (the deprecated optimize_panel_transferable() wrapper reproduces
  # it with the same seed).
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  partition_info <- NULL
  if (opts$partitioned) {
    x_list <- if (.is_cohort_list(x)) x else list(x)
    y_list <- if (.is_cohort_list(x)) y else list(y)
    if (!is.list(y_list) || length(y_list) != length(x_list)) {
      stop("`y` must be a list of response vectors matching `x`.", call. = FALSE)
    }
    # Partitioning subsets rows, which is meaningless for a
    # SummarizedExperiment; coerce cohorts to sample-by-feature matrices first.
    x_list <- lapply(x_list, .extract_feature_matrix, assay = assay)

    partitions <- .stratified_partition_cohorts(
      x_list, y_list, train_ratio, opts$val_ratio
    )
    partition_info <- partitions$partition_info

    prepare_partition <- function(part_x, part_y) {
      .prepare_cohort_inputs(
        part_x, part_y,
        assay = assay,
        transform = "none",
        feature_subset = feature_pool_base,
        feature_alignment = feature_alignment
      )
    }
    train_inputs <- prepare_partition(partitions$train_x, partitions$train_y)
    val_inputs <- if (opts$val_ratio > 0) {
      prepare_partition(partitions$val_x, partitions$val_y)
    } else {
      NULL
    }
    heldout_inputs <- if (heldout_ratio > 0) {
      prepare_partition(partitions$heldout_x, partitions$heldout_y)
    } else {
      NULL
    }
  } else {
    train_inputs <- .prepare_cohort_inputs(
      x, y,
      assay = assay,
      transform = "none",
      feature_subset = feature_pool_base,
      feature_alignment = feature_alignment
    )
    val_inputs <- NULL
    heldout_inputs <- NULL
  }

  # Train + validation pooled: what the returned result carries as its
  # training data, and what the pooled fitness modes (LOCO, rotating) search on.
  pooled <- .combine_partitions(train_inputs, val_inputs)
  x_raw <- pooled$x
  truth <- pooled$truth
  cohort <- pooled$cohort

  # Feature pool is now BASE features (not aggregated)
  feature_pool <- feature_pool_base

  if (!length(feature_pool)) {
    stop("`feature_pool` produced zero features.", call. = FALSE)
  }

  min_features_required <- .min_features_required(regularized, feature_transform)

  if (max_features < min_features_required) {
    stop("`max_features` must be at least ", min_features_required,
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
  ref_attr <- attr(x_raw, "reference_feature")
  x_pool <- x_raw[, feature_pool, drop = FALSE]
  if (!is.null(ref_attr)) {
    attr(x_pool, "reference_feature") <- ref_attr
  }
  decision_dim <- ncol(x_pool)  # Decision dimension = number of BASE features

  if (decision_dim > 200) {
    warning("Optimizing over more than 200 features may be slow; consider ",
            "reducing `feature_pool` for exploration.")
  }

  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  names(objective_directions) <- names(objectives)

  constraint_specs <- .normalize_constraints(constraints)

  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim, algorithm)
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)

  # Shared selector used by fitness evaluation (inside the factory) and by
  # history materialization here.
  panel_selector <- .make_panel_selector(
    feature_pool = colnames(x_pool),
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )
  select_base_features <- function(decision_vec) {
    panel_selector(decision_vec)$base_features
  }

  fitness_fn <- .make_optimize_panel_fitness(
    fitness_mode = fitness_mode,
    x_pool = x_pool,
    truth = truth,
    cohort = cohort,
    train_inputs = train_inputs,
    val_inputs = val_inputs,
    feature_pool = colnames(x_pool),
    max_features = max_features,
    objectives = objectives,
    constraints = constraints,
    regularized = regularized,
    regularized_alpha = regularized_alpha,
    feature_transform = feature_transform,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold,
    fitness_cv_folds = fitness_cv_folds,
    n_val_splits = n_val_splits,
    train_ratio = train_ratio,
    val_ratio = opts$val_ratio,
    seed = seed
  )
  objective_wrapper <- fitness_fn$wrapper

  # Optional per-generation history capture. We build a closure that pushes
  # each generation's population, fitness, and front into a parent-scope env.
  # The decision-weight matrix is needed to reconstruct base_features post-run
  # (we don't compute features inside the monitor to keep per-iter overhead low).
  history_buffer <- .make_history_buffer()
  history_monitor <- .make_history_monitor(
    history_buffer, objective_directions, names(objectives)
  )

  monitor_arg <- if (isTRUE(record_history)) history_monitor else FALSE

  nsga_result <- .run_nsga(
    algorithm = algorithm,
    fitness = objective_wrapper,
    n_objectives = length(objectives),
    decision_dim = decision_dim,
    nsga_args = nsga_args,
    monitor = monitor_arg,
    min_features = min_features_required,
    max_features = max_features,
    seed = seed,
    set_seed = FALSE
  )

  # Extract the Pareto front, re-evaluate, drop infeasible/dominated solutions,
  # and assemble the wide-format solutions data frame.
  solutions_df <- .build_pareto_solutions_df(
    nsga_result, fitness_fn$evaluate, objectives, objective_directions,
    infeasible_msg = "No solutions satisfied the supplied constraints. Consider relaxing them."
  )

  # Build control parameters to store
  control <- list(
    max_features = max_features,
    feature_pool = feature_pool,
    base_feature_pool = feature_pool_base,
    algorithm = algorithm,
    nsga2 = nsga_args,  # Keep nsga2 name for backward compatibility
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
    fitness_mode = fitness_mode,
    objective_directions = objective_directions,
    # Full objective descriptors (label/direction/fun, including any metric
    # params such as target specificity) so downstream evaluation can re-score
    # held-out data on the criteria the search actually optimised.
    objectives = objectives
  )
  if (fitness_mode == "cv") {
    control$fitness_cv_folds <- fitness_cv_folds
  }
  if (fitness_mode == "within_cohort_rotating") {
    control$n_val_splits <- n_val_splits
  }
  # Held-out data (and the NP parameters that consume it) only exist when a
  # held-out share was requested. evaluate_pareto_front() and calibrate_panel()
  # test these fields for NULL.
  if (heldout_ratio > 0) {
    control <- c(control, list(
      heldout_x = heldout_inputs$x,
      heldout_y = heldout_inputs$truth,
      heldout_cohort = heldout_inputs$cohort,
      np_alpha = np_alpha,
      np_delta = np_delta,
      partition_info = list(
        train_ratio = train_ratio,
        val_ratio = opts$val_ratio,
        heldout_ratio = heldout_ratio,
        partition_sizes = partition_info
      )
    ))
  }

  # Build training signature (counts describe the pooled train + val data)
  cohort_labels <- levels(cohort)
  cohort_counts <- setNames(as.list(as.integer(table(cohort))), cohort_labels)
  training_signature <- list(
    n = nrow(x_raw),
    p = ncol(x_raw),
    class_balance = table(truth),
    feature_pool_size = length(feature_pool),
    base_feature_pool_size = length(feature_pool_base),
    num_cohorts = length(cohort_labels),
    cohort_labels = cohort_labels,
    cohort_counts = cohort_counts
  )

  # Materialize per-generation history if it was captured. For each individual
  # we run the same selection logic the optimizer used to derive base_features,
  # so the recorded panels exactly match what NSGA was scoring.
  history_out <- if (isTRUE(record_history)) {
    .materialize_history(history_buffer, select_base_features)
  } else {
    list()
  }

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
    aggregated_cohort = cohort,
    history = history_out
  )
}
