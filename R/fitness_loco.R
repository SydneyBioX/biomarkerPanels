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
