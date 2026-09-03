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
#' @param seed Optional integer seed for reproducibility. When provided, sets
#'   the random seed before running NSGA optimization. This ensures
#'   reproducible results across runs. If `NULL` (default), no seed is set and
#'   results may vary between runs.
#' @param fitness_cv Logical; if `TRUE` (default), use cross-validation when
#'   evaluating candidate solutions during NSGA optimization. This prevents
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
                           fitness_cv = TRUE,
                           fitness_cv_folds = 5L,
                           regularized = TRUE,
                           regularized_alpha = 0.5,
                           selection_threshold = "adaptive",
                           record_history = FALSE) {
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
  .validate_positive_integer(fitness_cv_folds, "fitness_cv_folds", min = 2L)
  .validate_probability(regularized_alpha, "regularized_alpha", bounds = "closed")
  .validate_selection_threshold(selection_threshold)

  inputs_raw <- .prepare_cohort_inputs(x, y, assay = assay, transform = "none",
                                       feature_alignment = feature_alignment)
  raw_x_mat <- .ensure_feature_colnames(inputs_raw$x)
  raw_feature_names <- colnames(raw_x_mat)

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


  if (!fitness_cv) {
    warning("fitness_cv = FALSE: NSGA will use in-sample scoring, which risks ",
            "overfitting during optimization. Consider fitness_cv = TRUE for ",
            "more reliable panel selection.", call. = FALSE)
  }

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

  fitness_fn <- .make_cv_fitness(
    x = x_pool,
    y = truth,
    cohort = cohort,
    feature_pool = colnames(x_pool),
    max_features = max_features,
    objectives = objectives,
    constraints = constraints,
    regularized = regularized,
    alpha = regularized_alpha,
    feature_transform = feature_transform,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold,
    fitness_cv = fitness_cv,
    fitness_cv_folds = fitness_cv_folds
  )
  objective_wrapper <- fitness_fn$wrapper

  # Optional per-generation history capture. We build a closure that pushes
  # each generation's population, fitness, and front into a parent-scope env.
  # The decision-weight matrix is needed to reconstruct base_features post-run
  # (we don't compute features inside the monitor to keep per-iter overhead low).
  # The capture machinery is shared with optimize_panel_transferable() via
  # nsga_history_utils.R.
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
    set_seed = TRUE
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
    objective_directions = objective_directions,
    # Full objective descriptors (label/direction/fun, including any metric
    # params such as target specificity) so downstream evaluation can re-score
    # held-out data on the criteria the search actually optimised.
    objectives = objectives
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
