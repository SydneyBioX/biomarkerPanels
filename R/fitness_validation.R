#' Validation-Based NSGA Fitness
#'
#' Factory for the train-on-training / evaluate-on-validation fitness used by
#' [optimize_panel_transferable()]. Split out of that file to keep the main
#' orchestration function readable.
#'
#' @name fitness_validation
NULL

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
  if (is.null(min_features_required)) {
    min_features_required <- .min_features_required(regularized, feature_transform)
  }

  scaffold <- .make_fitness_scaffold(
    feature_pool = feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold,
    matrices = list(train = train_x, val = val_x),
    feature_transform = feature_transform,
    objectives = objectives,
    constraints = constraints,
    cache_max_entries = cache_max_entries
  )
  panel_selector <- scaffold$selector
  transform_panel <- scaffold$transform
  constraint_specs <- scaffold$constraint_specs
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
          val_scores <- .predict_panel_model(fit, x_val_sel, cohort = val_cohort)
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
        constraints_eval <- .evaluate_constraints(
          constraint_specs, val_y, val_scores,
          selected = selected_base_features, cohort = val_cohort, x = x_val_sel
        )
        constraint_results <- constraints_eval$results
        feasible <- constraints_eval$feasible

        metrics <- .evaluate_objectives(
          objectives, val_y, val_scores,
          selected = selected_base_features, cohort = val_cohort, x = x_val_sel
        )

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

  scaffold$finalize(evaluate_candidate, cache_fitness = cache_fitness)
}
