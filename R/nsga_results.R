#' NSGA Result Post-Processing
#'
#' Internal helper for turning a raw `rmoo` NSGA result into the wide-format
#' solutions data frame stored on an `OptimizationResult`.
#'
#' @name nsga_results
#' @keywords internal
NULL

#' Build the Pareto Solutions Data Frame from an NSGA Result
#'
#' Shared post-processing for [optimize_panel()] and
#' [optimize_panel_transferable()]. Extracts the rank-1 (Pareto-optimal)
#' population from an `rmoo` result, re-evaluates each candidate to recover its
#' features and objective metrics, drops constraint-infeasible solutions, filters
#' any solutions that became dominated on re-evaluation, and assembles the
#' wide-format solutions data frame (one row per solution, list columns for
#' `base_features`/`features`, one numeric column per objective).
#'
#' Re-evaluation can shift metrics relative to the NSGA run (e.g. non-deterministic
#' inner CV in glmnet), which can turn a rank-1 solution into a dominated one; the
#' post-filter removes those via [.filter_dominated()].
#'
#' @param nsga_result The object returned by `rmoo::nsga2()` / `rmoo::nsga3()`.
#' @param evaluate Function mapping a decision vector to a solution list with
#'   `feasible`, `metrics`, `base_features`, and `features` elements.
#' @param objectives Named list of objective specs (its names order the metric
#'   columns).
#' @param objective_directions Named character vector of objective directions,
#'   passed to [.filter_dominated()].
#' @param infeasible_msg Error message used when no solution is feasible.
#' @return A `data.frame` with `solution_id`, `base_features`/`features` list
#'   columns, and one column per objective.
#' @keywords internal
.build_pareto_solutions_df <- function(nsga_result, evaluate, objectives,
                                       objective_directions,
                                       infeasible_msg = "No solutions satisfied the supplied constraints.") {
  # Filter to Pareto-optimal solutions (front rank == 1). rmoo returns the full
  # population; a single solution comes back as a vector, so re-matrixify.
  optimal_idx <- which(nsga_result@front == 1)
  pareto_pop <- nsga_result@population[optimal_idx, , drop = FALSE]
  if (is.null(dim(pareto_pop))) {
    pareto_pop <- matrix(pareto_pop, nrow = 1)
  }

  solutions <- lapply(seq_len(nrow(pareto_pop)), function(i) {
    evaluate(pareto_pop[i, ])
  })

  feasible_vec <- vapply(solutions, `[[`, logical(1), "feasible")
  if (!any(feasible_vec)) {
    stop(infeasible_msg, call. = FALSE)
  }
  solutions <- solutions[feasible_vec]

  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  # Post-filter dominated solutions (see .filter_dominated docs).
  nondom_idx <- .filter_dominated(metric_matrix, objective_directions)
  solutions <- solutions[nondom_idx]
  metric_matrix <- metric_matrix[nondom_idx, , drop = FALSE]

  # Build solutions data frame in wide format (one row per solution).
  solutions_df <- data.frame(
    solution_id = seq_along(solutions),
    stringsAsFactors = FALSE
  )
  solutions_df$base_features <- I(lapply(solutions, `[[`, "base_features"))
  solutions_df$features <- I(lapply(solutions, `[[`, "features"))
  for (obj_name in names(objectives)) {
    solutions_df[[obj_name]] <- metric_matrix[, obj_name]
  }

  solutions_df
}
