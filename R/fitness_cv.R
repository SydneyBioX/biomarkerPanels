#' CV / In-Sample NSGA Fitness
#'
#' Factory for the cross-validated (or in-sample) fitness used by
#' [optimize_panel()] by default. Sits alongside the validation,
#' rotating-validation, and LOCO factories (`R/fitness_validation.R`,
#' `R/fitness_rotating.R`, `R/fitness_loco.R`) as the fourth fitness strategy:
#' this one trains and scores on the same (optionally k-folded) data rather
#' than a held-out partition. Split out of `optimize_panel()` to keep the main
#' orchestration function readable.
#'
#' @name fitness_cv
NULL

#' Create CV / In-Sample Fitness Function (Base-Feature-First)
#'
#' Factory that returns a fitness function for rmoo that selects base
#' features, applies feature transforms on-the-fly, and scores candidates via
#' k-fold cross-validation (out-of-fold predictions) or in-sample prediction.
#'
#' @param x Feature matrix (base features).
#' @param y Response factor.
#' @param cohort Cohort factor.
#' @param feature_pool Character vector of base feature names.
#' @param max_features Maximum base features to include.
#' @param objectives Objective list from define_objectives().
#' @param constraints Constraint list.
#' @param regularized Whether to use regularized regression.
#' @param alpha Elastic net mixing parameter.
#' @param feature_transform Name of the feature transform to apply.
#' @param min_features_required Minimum number of base features.
#' @param selection_threshold Either `"adaptive"` or a fixed numeric threshold.
#' @param fitness_cv Logical; if `TRUE`, score candidates via k-fold
#'   cross-validation. Falls back to in-sample scoring (with a warning) when
#'   there are too few samples for `fitness_cv_folds` folds.
#' @param fitness_cv_folds Number of cross-validation folds when
#'   `fitness_cv = TRUE`.
#' @return List with `wrapper` (fitness for rmoo) and `evaluate` (full eval).
#' @noRd
.make_cv_fitness <- function(
  x, y, cohort,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL,
  selection_threshold = "adaptive",
  fitness_cv = TRUE,
  fitness_cv_folds = 5L
) {
  if (is.null(min_features_required)) {
    min_features_required <- .min_features_required(regularized, feature_transform)
  }

  scaffold <- .make_fitness_scaffold(
    feature_pool = feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold,
    matrices = list(x = x),
    feature_transform = feature_transform,
    objectives = objectives,
    constraints = constraints
  )
  panel_selector <- scaffold$selector
  transform_panel <- scaffold$transform
  constraint_specs <- scaffold$constraint_specs

  # Create CV folds for fitness evaluation (if enabled)
  cv_fold_ids <- NULL
  if (fitness_cv) {
    n_samples <- nrow(x)
    if (n_samples < fitness_cv_folds * 2L) {
      warning("Too few samples for ", fitness_cv_folds, "-fold CV in fitness evaluation. ",
              "Falling back to in-sample scoring.", call. = FALSE)
      fitness_cv <- FALSE
    } else {
      cv_fold_ids <- .create_stratified_folds(y, fitness_cv_folds)
    }
  }

  cv_glm_design_terms <- NULL
  if (fitness_cv && !is.null(cv_fold_ids) && !regularized) {
    cv_glm_design_terms <- .prepare_cv_glm_design_terms(cohort, cv_fold_ids)
  }
  in_sample_glm_design_terms <- NULL
  if (!fitness_cv && !regularized) {
    in_sample_glm_design_terms <- .prepare_glm_design_terms(
      cohort_train = cohort,
      cohort_new = cohort,
      predict_cohort = "observed"
    )
  }

  evaluate_candidate <- function(decision_vec = NULL, selection = NULL) {
    if (is.null(selection)) {
      selection <- panel_selector(decision_vec)
    }
    selected_base_features <- selection$base_features

    transformed <- transform_panel(selected_base_features)
    x_selected <- transformed$matrices$x
    selected_features <- transformed$features

    # Compute scores using CV (out-of-fold predictions) or in-sample
    if (fitness_cv && !is.null(cv_fold_ids)) {
      # Cross-validation based scoring: prevents overfitting
      scores <- .compute_cv_scores(
        x_selected = x_selected,
        truth = y,
        fold_ids = cv_fold_ids,
        cohort = cohort,
        regularized = regularized,
        alpha = alpha,
        glm_design_terms = cv_glm_design_terms
      )
    } else {
      # In-sample scoring (for backward compatibility or small datasets)
      scores <- if (!regularized) {
        .fit_predict_binomial_glm(
          x_train = x_selected,
          truth = y,
          x_new = x_selected,
          cohort_train = cohort,
          cohort_new = cohort,
          predict_cohort = "observed",
          design_terms = in_sample_glm_design_terms
        )
      } else {
        .default_scoring_fn(
          x_selected = x_selected,
          selected_features = selected_features,
          truth = y,
          cohort = cohort,
          regularized = regularized,
          alpha = alpha
        )
      }
    }

    if (!is.numeric(scores) || length(scores) != nrow(x)) {
      stop("Scoring must return a numeric vector matching the number of samples.",
           call. = FALSE)
    }

    # Evaluate constraints
    constraints_eval <- .evaluate_constraints(
      constraint_specs, y, scores,
      selected = selected_base_features, cohort = cohort, x = x_selected
    )
    constraint_results <- constraints_eval$results
    feasible <- constraints_eval$feasible

    # Evaluate objectives on the transformed scoring matrix.
    metrics <- .evaluate_objectives(
      objectives, y, scores,
      selected = selected_base_features, cohort = cohort, x = x_selected
    )

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

  scaffold$finalize(evaluate_candidate)
}
