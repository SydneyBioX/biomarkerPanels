#' Argument Validation for optimize_panel()
#'
#' Front-half checks for [optimize_panel()], kept in their own file so the
#' orchestration pipeline stays readable. Everything here is pure: it validates,
#' normalises, and returns the resolved values the pipeline then uses.
#'
#' @name optimize_panel_validate
#' @noRd
NULL

#' Validate and Normalise optimize_panel() Arguments
#'
#' Resolves `fitness_mode` (including the deprecated `fitness_cv` alias),
#' checks the numeric scalars, and works out the partitioning shares. Mode /
#' partition compatibility is checked here too, before any expensive work: a
#' mode that scores on a validation partition needs `val_ratio > 0`, and the
#' single-partition modes must not be handed one.
#'
#' @param fitness_mode Requested fitness strategy (may be the full default
#'   vector, as passed straight from the formals).
#' @param fitness_cv Deprecated logical alias for `fitness_mode`.
#' @param fitness_cv_folds Number of CV folds.
#' @param n_val_splits Number of rotating validation splits.
#' @param train_ratio,val_ratio Partition shares.
#' @param regularized_alpha Elastic-net mixing parameter.
#' @param selection_threshold `"adaptive"` or a numeric gate in (0, 1).
#' @param feature_transform Name of a registered feature transform.
#' @param np_alpha,np_delta Neyman-Pearson calibration parameters.
#' @param seed Optional integer seed.
#' @return List with the normalised `fitness_mode`, `fitness_cv_folds`,
#'   `n_val_splits`, `train_ratio`, `val_ratio`, `heldout_ratio`, `seed`, and
#'   the logical `partitioned` telling the caller whether to split at all.
#' @noRd
.validate_optimize_panel_args <- function(fitness_mode,
                                          fitness_cv,
                                          fitness_cv_folds,
                                          n_val_splits,
                                          train_ratio,
                                          val_ratio,
                                          regularized_alpha,
                                          selection_threshold,
                                          feature_transform,
                                          np_alpha,
                                          np_delta,
                                          seed) {
  fitness_modes <- c("cv", "in_sample", "within_cohort_val",
                     "within_cohort_rotating", "loco")

  if (!is.null(fitness_cv)) {
    if (!is.logical(fitness_cv) || length(fitness_cv) != 1L || is.na(fitness_cv)) {
      stop("`fitness_cv` must be a single logical value.", call. = FALSE)
    }
    warning(
      "`fitness_cv` is deprecated; use `fitness_mode` instead. Mapping ",
      "fitness_cv = ", fitness_cv, " to fitness_mode = ",
      if (fitness_cv) "\"cv\"" else "\"in_sample\"", ".",
      call. = FALSE
    )
    fitness_mode <- if (fitness_cv) "cv" else "in_sample"
  }
  fitness_mode <- match.arg(fitness_mode, fitness_modes)

  if (!is.character(feature_transform) || length(feature_transform) != 1L ||
      !nzchar(feature_transform)) {
    stop("`feature_transform` must be a single non-empty character string.", call. = FALSE)
  }
  if (!exists(feature_transform, envir = .transform_registry, inherits = FALSE)) {
    stop(
      "Unknown feature transform '", feature_transform, "'. ",
      "Available: ", paste(ls(.transform_registry), collapse = ", "), ". ",
      "Use register_feature_transform() to add custom transforms.",
      call. = FALSE
    )
  }

  fitness_cv_folds <- .validate_positive_integer(fitness_cv_folds, "fitness_cv_folds", min = 2L)
  .validate_probability(regularized_alpha, "regularized_alpha", bounds = "closed")
  .validate_selection_threshold(selection_threshold)
  .validate_probability(np_alpha, "np_alpha", bounds = "open")
  .validate_probability(np_delta, "np_delta", bounds = "open")

  if (fitness_mode == "within_cohort_rotating") {
    n_val_splits <- .validate_positive_integer(n_val_splits, "n_val_splits", min = 2L)
  }

  # `.run_nsga()` no longer seeds, so the check has to live here.
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("`seed` must be a single integer value.", call. = FALSE)
  }

  .validate_partition_ratios(train_ratio, val_ratio)
  heldout_ratio <- .heldout_ratio(train_ratio, val_ratio)
  partitioned <- train_ratio < 1

  needs_val <- fitness_mode %in% c("within_cohort_val", "within_cohort_rotating")
  if (needs_val && val_ratio <= 0) {
    stop("`fitness_mode = \"", fitness_mode, "\"` scores candidates on a ",
         "validation partition and therefore requires `val_ratio > 0`.",
         call. = FALSE)
  }
  if (fitness_mode %in% c("cv", "in_sample") && val_ratio > 0) {
    stop("`fitness_mode = \"", fitness_mode, "\"` trains and scores on a single ",
         "partition; set `val_ratio = 0` (validation samples would be unused).",
         call. = FALSE)
  }

  list(
    fitness_mode = fitness_mode,
    fitness_cv_folds = fitness_cv_folds,
    n_val_splits = n_val_splits,
    train_ratio = train_ratio,
    val_ratio = val_ratio,
    heldout_ratio = heldout_ratio,
    partitioned = partitioned,
    seed = seed
  )
}
