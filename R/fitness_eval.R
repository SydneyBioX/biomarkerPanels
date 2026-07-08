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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
