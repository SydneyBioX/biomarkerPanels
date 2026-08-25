#' Leave-One-Cohort-Out NSGA Fitness
#'
#' Factory for the leave-one-cohort-out (LOCO) fitness used by
#' [optimize_panel_transferable()]. Split out of that file to keep the main
#' orchestration function readable.
#'
#' @name fitness_loco
NULL

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
#' @return List with `wrapper` (vectorized NSGA fitness) and `evaluate`
#'   (single-candidate evaluator returning full diagnostics).
#' @noRd
.make_loco_fitness <- function(
  x, y, cohort,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL,
  selection_threshold = "adaptive"
) {
  if (is.null(min_features_required)) {
    min_features_required <- .min_features_required(regularized, feature_transform)
  }

  cohort <- droplevels(factor(cohort))
  cohort_levels <- levels(cohort)
  if (length(cohort_levels) < 2L) {
    stop("LOCO fitness requires at least 2 cohorts.", call. = FALSE)
  }

  # Precompute row indices per cohort for speed in the fitness loop
  fold_rows <- lapply(cohort_levels, function(lv) which(cohort == lv))
  names(fold_rows) <- cohort_levels

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
          # Skip fold if training data collapses to a single class; the
          # affected cohort is then silently absent from the metrics, so warn
          # (mirrors the rotating-validation guard).
          if (length(unique(train_y)) < 2L) {
            warning(
              "LOCO fold for cohort '", lv, "' skipped: training data ",
              "contains a single class. Metrics for this candidate exclude ",
              "that cohort.",
              call. = FALSE
            )
            next
          }

          x_tr <- x_transformed[train_idx, , drop = FALSE]
          x_te <- x_transformed[test_idx, , drop = FALSE]
          coh_tr <- droplevels(cohort[train_idx])

          if (regularized) {
            fit <- .fit_final_model_regularized(x_tr, train_y, coh_tr, alpha = alpha)
            preds <- .predict_panel_model(fit, x_te, cohort = cohort[test_idx])
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

        constraints_eval <- .evaluate_constraints(
          constraint_specs, eval_truth, eval_scores,
          selected = selected_base_features, cohort = eval_cohort, x = NULL
        )
        constraint_results <- constraints_eval$results
        feasible <- constraints_eval$feasible

        metrics <- .evaluate_objectives(
          objectives, eval_truth, eval_scores,
          selected = selected_base_features, cohort = eval_cohort, x = NULL
        )

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

  scaffold$finalize(evaluate_candidate)
}
