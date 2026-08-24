#' Evaluate Pareto Solutions on Held-Out Data
#'
#' Fits and evaluates every Pareto-optimal solution from [optimize_panel()] on
#' held-out validation data. This provides the metrics data frame used by
#' [plot_pareto_front()] without requiring a plot.
#'
#' In addition to the always-present descriptive columns (`sensitivity`,
#' `specificity`, `auc` at a 0.5 cutoff), each solution is re-scored on the
#' objectives the optimisation actually targeted. Those are recovered from
#' `optimization_result@control$objectives` (stored by [optimize_panel()] and
#' [optimize_panel_transferable()]) and reported in columns prefixed `obj_`,
#' so `obj_auc` or `obj_sensitivity_at_specificity` are directly comparable to
#' the training-side values in `solutions(optimization_result)`.
#'
#' Objectives whose metric is cohort-aware (e.g. `cohort_auc_gap`,
#' `min_cohort_auc`) cannot be computed when the evaluation data contains a
#' single cohort; those columns are returned as `NA` with a warning. Supply
#' `x` as a list of cohorts, or pass `cohort`, to make them computable.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @param x Held-out validation data: a matrix, data.frame,
#'   `SummarizedExperiment`, or list of such objects. When `NULL` (default),
#'   automatically uses held-out data stored by
#'   [optimize_panel_transferable()].
#' @param y Held-out validation labels: a factor (or list of factors when `x`
#'   is a list). When `NULL` (default), automatically uses held-out labels
#'   stored by [optimize_panel_transferable()].
#' @param cohort Optional factor identifying cohort membership for each sample
#'   when `x` is a single matrix. Ignored when `x` is a list of cohorts. When
#'   `x`/`y` are taken from the stored held-out partition, the stored cohort
#'   labels are used unless overridden here.
#' @param objectives Optional named list of objective descriptors (as returned
#'   by [define_objectives()]) overriding the objectives stored on
#'   `optimization_result`. Set to `NA` to skip objective re-scoring entirely
#'   and return only the descriptive columns.
#' @param on_error One of `"warn"` or `"stop"`. Controls behaviour when a
#'   solution fails evaluation. `"warn"` (default) skips the solution with a
#'   warning; `"stop"` raises an error immediately.
#' @param regularized Logical; passed to [fit_panel()]. Default `TRUE`.
#' @param verbose Logical; print progress messages. Default `interactive()`.
#' @return A data frame with one row per successfully evaluated solution and
#'   columns `solution_id`, `n_features`, `n_base_features`, `sensitivity`,
#'   `specificity`, `auc`, plus one `obj_<name>` column per optimisation
#'   objective. Attributes `objective_labels` and `objective_directions` map
#'   the `obj_` columns back to their metric labels and optimisation
#'   directions.
#' @export
#' @seealso [optimize_panel()], [optimize_panel_transferable()], [fit_panel()],
#'   [evaluate_panel()], [plot_pareto_front()]
#' @examples
#' \dontrun{
#' opt <- optimize_panel(x_train, y_train, objectives = define_objectives())
#' evaluate_pareto_solutions(opt, x_test, y_test)
#' }
evaluate_pareto_solutions <- function(optimization_result,
                                      x = NULL,
                                      y = NULL,
                                      cohort = NULL,
                                      objectives = NULL,
                                      on_error = c("warn", "stop"),
                                      regularized = TRUE,
                                      verbose = interactive()) {
  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult.", call. = FALSE)
  }
  on_error <- match.arg(on_error)
  ctrl <- optimization_result@control

  # Resolve held-out data
  if (is.null(x) || is.null(y)) {
    if (!is.null(ctrl$heldout_x) && !is.null(ctrl$heldout_y)) {
      x <- ctrl$heldout_x
      y <- ctrl$heldout_y
      if (is.null(cohort)) cohort <- ctrl$heldout_cohort
      if (verbose) message("Using held-out data stored in OptimizationResult.")
    } else {
      stop("Provide held-out `x` and `y`, or use a result from ",
           "optimize_panel_transferable().", call. = FALSE)
    }
  }
  if (is.list(x) && !is.null(cohort)) {
    warning("`cohort` is ignored when `x` is supplied as a list of cohorts.",
            call. = FALSE)
    cohort <- NULL
  }

  n_sol <- n_solutions(optimization_result)
  if (n_sol == 0L) {
    stop("OptimizationResult contains no solutions.", call. = FALSE)
  }

  # Recover the objectives the search optimised, then make each one safe to
  # evaluate on data it may not support (single cohort, degenerate labels).
  eval_objectives <- .resolve_stored_objectives(optimization_result, objectives)
  eval_objectives <- .prepare_pareto_objectives(
    eval_objectives, .count_eval_cohorts(x, cohort)
  )
  obj_names <- names(eval_objectives)
  obj_cols <- if (length(obj_names)) {
    make.names(paste0("obj_", obj_names), unique = TRUE)
  } else {
    character()
  }

  # evaluate_panel() requires a non-empty objective list; the descriptive
  # sensitivity/specificity/auc columns are read off the ROC summary instead,
  # so they stay comparable across calls regardless of the objectives used.
  panel_objectives <- if (length(eval_objectives)) {
    eval_objectives
  } else {
    define_objectives(metrics = "auc")
  }

  # Evaluate each solution: fit -> evaluate -> extract metrics
  rows <- vector("list", n_sol)
  for (i in seq_len(n_sol)) {
    sid <- solutions(optimization_result)$solution_id[i]
    if (verbose) message(sprintf("[%d/%d] solution %d", i, n_sol, sid))

    rows[[i]] <- tryCatch({
      panel <- fit_panel(optimization_result, solution_id = sid,
                         regularized = regularized)
      ev <- evaluate_panel(panel, x, y, objectives = panel_objectives,
                           cohort = cohort)

      row <- data.frame(
        solution_id = sid,
        n_features = length(panel@features),
        n_base_features = length(panel@base_features),
        sensitivity = ev$roc$highlight$sensitivity,
        specificity = ev$roc$highlight$specificity,
        auc = ev$roc$auc,
        stringsAsFactors = FALSE
      )
      for (j in seq_along(obj_names)) {
        value <- ev$metrics[[obj_names[j]]]
        row[[obj_cols[j]]] <- if (is.null(value)) NA_real_ else as.numeric(value)
      }
      row
    }, error = function(e) {
      if (on_error == "stop") {
        stop("Failed on solution ", sid, ": ", conditionMessage(e),
             call. = FALSE)
      }
      warning("Solution ", sid, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    stop("All solutions failed evaluation.", call. = FALSE)
  }
  out <- do.call(rbind, rows)

  if (length(obj_names)) {
    attr(out, "objective_labels") <- stats::setNames(
      vapply(eval_objectives, function(o) as.character(o$label)[1L],
             character(1)),
      obj_cols
    )
    attr(out, "objective_directions") <- stats::setNames(
      vapply(eval_objectives, function(o) as.character(o$direction)[1L],
             character(1)),
      obj_cols
    )
  }
  out
}

#' Recover the Objectives an Optimization Run Targeted
#'
#' Returns the objective descriptors stored on an `OptimizationResult`, falling
#' back to rebuilding them from `control$objective_directions` for results
#' created before objectives were stored. Rebuilt objectives use registry
#' defaults, so any metric parameters supplied at optimisation time (e.g. the
#' target specificity of `sensitivity_at_specificity`) are lost; this is
#' warned about.
#'
#' @param optimization_result An `OptimizationResult`.
#' @param objectives User override: `NULL` to use the stored objectives, `NA`
#'   to skip objective re-scoring, or a named list of objective descriptors.
#' @return A (possibly empty) named list of objective descriptors.
#' @noRd
.resolve_stored_objectives <- function(optimization_result, objectives = NULL) {
  if (length(objectives) == 1L && !is.list(objectives) && is.na(objectives)) {
    return(list())
  }
  if (!is.null(objectives)) {
    if (!is.list(objectives) || is.null(names(objectives)) ||
        !all(nzchar(names(objectives)))) {
      stop("`objectives` must be a named list of objective descriptors.",
           call. = FALSE)
    }
    return(objectives)
  }

  ctrl <- optimization_result@control
  stored <- ctrl$objectives
  if (is.list(stored) && length(stored) && !is.null(names(stored))) {
    return(stored)
  }

  # Backward compatibility: older results kept only the direction vector.
  directions <- ctrl$objective_directions
  if (is.null(directions) || !length(directions)) {
    return(list())
  }
  registry <- metric_registry()
  known <- intersect(names(directions), names(registry))
  dropped <- setdiff(names(directions), known)
  if (length(dropped)) {
    warning("Cannot rebuild objective(s) not in the metric registry: ",
            paste(dropped, collapse = ", "), ".", call. = FALSE)
  }
  if (!length(known)) {
    return(list())
  }
  warning("OptimizationResult stores no objective definitions; rebuilding ",
          "from objective_directions using registry defaults. Any metric ",
          "parameters used during optimisation (e.g. target specificity) are ",
          "not recoverable.", call. = FALSE)
  build_objectives(known, directions = directions[known])
}

#' Count Cohorts Available in Evaluation Data
#'
#' @param x Evaluation data (matrix-like or list of cohorts).
#' @param cohort Optional cohort vector for single-matrix input.
#' @return Integer count of distinct cohorts.
#' @noRd
.count_eval_cohorts <- function(x, cohort = NULL) {
  if (is.list(x) && !is.data.frame(x)) {
    return(length(x))
  }
  if (is.null(cohort)) {
    return(1L)
  }
  nlevels(factor(cohort))
}

#' Make Objectives Safe to Evaluate on Held-Out Data
#'
#' Wraps each objective so that (a) cohort-aware metrics return `NA` when the
#' evaluation data holds a single cohort — where they would otherwise report a
#' degenerate value such as a zero cohort gap — and (b) an error inside one
#' metric yields `NA` for that objective rather than discarding the whole
#' solution.
#'
#' @param objectives Named list of objective descriptors.
#' @param n_cohorts Number of cohorts present in the evaluation data.
#' @return The objective list with wrapped `fun` elements.
#' @noRd
.prepare_pareto_objectives <- function(objectives, n_cohorts) {
  if (!length(objectives)) {
    return(objectives)
  }
  registry <- metric_registry()
  unavailable <- character()

  prepared <- lapply(names(objectives), function(name) {
    obj <- objectives[[name]]
    entry <- registry[[name]]
    needs_cohort <- !is.null(entry) &&
      "cohort" %in% names(formals(entry$fun))
    if (needs_cohort && n_cohorts < 2L) {
      unavailable <<- c(unavailable, name)
      obj$fun <- function(truth, estimate, selected = NULL, ...) NA_real_
      return(obj)
    }
    inner <- obj$fun
    obj$fun <- function(truth, estimate, selected = NULL, ...) {
      value <- tryCatch(
        inner(truth, estimate, selected = selected, ...),
        error = function(e) NA_real_
      )
      if (!is.numeric(value) || length(value) != 1L) NA_real_ else as.numeric(value)
    }
    obj
  })
  names(prepared) <- names(objectives)

  if (length(unavailable)) {
    warning("Cohort-aware objective(s) cannot be computed on single-cohort ",
            "evaluation data and are returned as NA: ",
            paste(unavailable, collapse = ", "), ".", call. = FALSE)
  }
  prepared
}
