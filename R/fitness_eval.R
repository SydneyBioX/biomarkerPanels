#' Objective and Constraint Evaluation
#'
#' Internal helpers that evaluate the registered constraint and objective
#' functions for a single scored candidate. Shared by [optimize_panel()] and the
#' validation / rotating-validation / LOCO fitness modules so the objective and
#' constraint calling conventions live in one place.
#'
#' A load-bearing design point encoded here: `metric_num_features` (and any other
#' count metric) receives the BASE features via `selected` so it counts genes,
#' not transformed columns such as pairwise ratios; performance metrics receive
#' the transformed scoring matrix via `x`. Other metrics ignore whichever
#' argument they do not need.
#'
#' @name fitness_eval
#' @noRd
NULL

#' Evaluate Hard Constraints for a Candidate
#'
#' Runs each constraint's predicate and reports both the per-constraint results
#' and whether the candidate is feasible (all constraints satisfied). Matches the
#' previous inline behavior exactly: an empty constraint set yields
#' `results = logical(0)` and `feasible = TRUE`.
#'
#' @param constraint_specs List of normalized constraint specs, each with a
#'   `$fun` predicate and a `$label`.
#' @param truth Binary response factor.
#' @param scores Numeric predicted scores.
#' @param selected Base features selected for the candidate.
#' @param cohort Optional cohort factor.
#' @param x Optional (transformed) scoring matrix.
#' @return A list with `results` (named logical vector) and `feasible` (scalar
#'   logical).
#' @noRd
.evaluate_constraints <- function(constraint_specs, truth, scores, selected,
                                  cohort = NULL, x = NULL) {
  results <- if (length(constraint_specs)) {
    setNames(
      vapply(seq_along(constraint_specs), function(j) {
        isTRUE(constraint_specs[[j]]$fun(
          truth = truth,
          scores = scores,
          selected = selected,
          cohort = cohort,
          x = x
        ))
      }, logical(1)),
      vapply(constraint_specs, `[[`, character(1), "label")
    )
  } else {
    logical(0)
  }

  list(
    results = results,
    feasible = if (length(results)) all(results) else TRUE
  )
}

#' Evaluate Objectives for a Candidate
#'
#' Computes each objective metric on the scored candidate, returning a named
#' numeric vector aligned with `objectives`.
#'
#' @inheritParams .evaluate_constraints
#' @param objectives Named list of objective specs, each with a `$fun` metric.
#' @return A named numeric vector of metric values, one per objective.
#' @noRd
.evaluate_objectives <- function(objectives, truth, scores, selected,
                                 cohort = NULL, x = NULL) {
  vapply(objectives, function(obj) {
    obj$fun(
      truth,
      scores,
      selected = selected,
      cohort = cohort,
      x = x
    )
  }, numeric(1))
}

# Large finite penalty instead of Inf — NSGA-III normalization produces NaN
# from Inf values, causing "missing value where TRUE/FALSE needed" errors.
.FITNESS_PENALTY <- 1e6

#' Build the Shared Fitness Scaffold
#'
#' Every fitness implementation (the inline [optimize_panel()] fitness and the
#' validation / rotating-validation / LOCO factories) shares the same skeleton:
#' a panel selector, an on-the-fly panel transformer, an objective cache, and
#' the candidate-to-objective conversion handed to rmoo. This factory owns that
#' skeleton; each fitness strategy supplies only its `evaluate_candidate`
#' scoring body.
#'
#' @param feature_pool Character vector of base feature names.
#' @param max_features Maximum base features to include.
#' @param min_features_required Minimum number of base features.
#' @param selection_threshold Either `"adaptive"` or a fixed numeric threshold.
#' @param matrices Named list of raw base-feature matrices handed to
#'   `.make_panel_transformer()`.
#' @param feature_transform Name of the feature transform to apply.
#' @param objectives Objective list from [define_objectives()].
#' @param constraints Constraint list.
#' @return List with `selector`, `transform`, `cache`, `directions`,
#'   `constraint_specs`, and `finalize(evaluate_candidate)`,
#'   which wraps a scoring body into the `list(wrapper, evaluate)` contract the
#'   optimizers consume.
#' @noRd
.make_fitness_scaffold <- function(feature_pool, max_features,
                                   min_features_required,
                                   selection_threshold,
                                   matrices, feature_transform,
                                   objectives, constraints) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)

  panel_selector <- .make_panel_selector(
    feature_pool = feature_pool,
    max_features = max_features,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )
  transform_panel <- .make_panel_transformer(
    matrices = matrices,
    feature_transform = feature_transform
  )
  objective_cache <- .new_fitness_cache()

  finalize <- function(evaluate_candidate) {
    evaluate_single <- function(decision_vec = NULL, selection = NULL) {
      evaluated <- evaluate_candidate(decision_vec = decision_vec,
                                      selection = selection)
      if (length(constraint_specs) && !evaluated$feasible) {
        return(rep(.FITNESS_PENALTY, length(objectives)))
      }
      .convert_metrics_to_objectives(evaluated$metrics, objective_directions,
                                     penalty = .FITNESS_PENALTY)
    }

    # rmoo fitness: receives a matrix (rows = individuals) or a single vector;
    # `...` absorbs extra arguments rmoo passes (e.g. reference_dirs).
    wrapper <- function(x, ...) {
      .evaluate_fitness_population(
        x = x,
        selector = panel_selector,
        evaluate_selection = function(selection) {
          evaluate_single(selection = selection)
        },
        n_objectives = length(objectives),
        cache = objective_cache
      )
    }

    list(wrapper = wrapper, evaluate = evaluate_candidate)
  }

  list(
    selector = panel_selector,
    transform = transform_panel,
    cache = objective_cache,
    directions = objective_directions,
    constraint_specs = constraint_specs,
    finalize = finalize
  )
}
