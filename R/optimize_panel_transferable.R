#' Optimize Biomarker Panels with Transferability Focus (Deprecated)
#'
#' @description
#' Deprecated thin wrapper. Use [optimize_panel()] directly, with
#' `fitness_mode`, `train_ratio`, and `val_ratio`: the train/validation/
#' held-out split and the `"within_cohort_val"` / `"within_cohort_rotating"` /
#' `"loco"` fitness strategies this function used to implement are now part of
#' [optimize_panel()] itself (a held-out share exists whenever
#' `train_ratio + val_ratio < 1`).
#'
#' When `feature_pool` is `NULL`, this wrapper still partitions the data and
#' runs `select_transferable_features()` on the training partition only, so
#' the old leakage-safe auto-selection behaviour survives for one release;
#' pass `feature_pool` explicitly (or call `select_transferable_features()`
#' yourself and pass its result) to skip that step and call [optimize_panel()]
#' directly instead.
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
#' @param seed Optional integer seed for reproducibility. When `NULL` and
#'   `feature_pool` is also `NULL`, a seed is drawn and reused for both the
#'   feature-selection partition and the [optimize_panel()] call, so the two
#'   partitions match and held-out samples never leak into feature selection.
#' @param regularized Logical; if `TRUE` (default), use regularized regression.
#' @param regularized_alpha Elastic net mixing parameter (default 0.5).
#' @param selection_threshold Either a fixed numeric threshold in (0,1) for
#'   feature selection, or `"adaptive"` (default) to use gap-based selection
#'   that allows variable panel sizes in the Pareto front. See
#'   [optimize_panel()] for details.
#' @param n_top_features Number of top features to select via
#'   `select_transferable_features()` when `feature_pool` is NULL (default 50).
#' @param fitness_mode Fitness-evaluation strategy: `"within_cohort_val"`
#'   (default), `"within_cohort_rotating"`, or `"loco"`. See
#'   [optimize_panel()] for the description of each.
#' @param n_val_splits Number of rotating train/val splits to pre-compute
#'   when `fitness_mode = "within_cohort_rotating"` (default 10).
#' @param record_history Logical; if `TRUE`, capture the population, fitness,
#'   and NSGA rank at every generation and attach the result to the returned
#'   `OptimizationResult`. Retrieve via [nsga_history()]. Default `FALSE`.
#' @return An `OptimizationResult` containing Pareto-optimal solutions and
#'   held-out data for downstream calibration. Use [summarize_solutions()] to
#'   inspect solutions, [fit_panel()] to fit a model, [calibrate_panel()] for
#'   NP threshold calibration, and [evaluate_panel()] for validation.
#' @export
#' @seealso [optimize_panel()], [fit_panel()], [calibrate_panel()],
#'   [evaluate_panel()], [summarize_solutions()]
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
  record_history = FALSE
) {
  .Deprecated(
    "optimize_panel",
    package = "biomarkerPanels",
    msg = paste(
      "optimize_panel_transferable() is deprecated; call optimize_panel()",
      "with fitness_mode, train_ratio, and val_ratio instead (a held-out",
      "share exists whenever train_ratio + val_ratio < 1). See ?optimize_panel."
    )
  )

  algorithm <- match.arg(algorithm)
  fitness_mode <- match.arg(fitness_mode)

  if (!.is_cohort_list(x) || length(x) < 1L) {
    stop("`x` must be a list of cohort matrices (one element per cohort).", call. = FALSE)
  }
  if (!is.list(y) || length(y) != length(x)) {
    stop("`y` must be a list of response vectors matching `x`.", call. = FALSE)
  }

  if (is.null(feature_pool)) {
    # Draw and reuse a seed so the partition used for feature selection here
    # matches the partition optimize_panel() re-derives below (it seeds again,
    # the same way, immediately before partitioning) -- otherwise the two
    # splits would differ and held-out samples could leak into the pool.
    if (is.null(seed)) {
      seed <- sample.int(.Machine$integer.max, 1L)
    }

    .validate_partition_ratios(train_ratio, val_ratio)
    set.seed(as.integer(seed))
    x_list <- lapply(x, .extract_feature_matrix, assay = assay)
    partitions <- .stratified_partition_cohorts(x_list, y, train_ratio, val_ratio)

    message("Selecting transferable features from training data...")
    fs_result <- select_transferable_features(
      partitions$train_x,
      partitions$train_y,
      n_features = n_top_features
    )
    feature_pool <- fs_result$features
    message("Selected ", length(feature_pool), " features for optimization.")
  }

  optimize_panel(
    x = x, y = y,
    objectives = objectives,
    max_features = max_features,
    feature_pool = feature_pool,
    feature_transform = feature_transform,
    feature_alignment = feature_alignment,
    constraints = constraints,
    algorithm = algorithm,
    nsga_control = nsga_control,
    assay = assay,
    seed = seed,
    fitness_mode = fitness_mode,
    n_val_splits = n_val_splits,
    train_ratio = train_ratio,
    val_ratio = val_ratio,
    np_alpha = np_alpha,
    np_delta = np_delta,
    regularized = regularized,
    regularized_alpha = regularized_alpha,
    selection_threshold = selection_threshold,
    record_history = record_history
  )
}
