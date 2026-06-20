#' Evaluate Pareto Solutions on Held-Out Data
#'
#' Fits and evaluates every Pareto-optimal solution from [optimize_panel()] on
#' held-out validation data. This provides the metrics data frame used by
#' [plot_pareto_front()] without requiring a plot.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @param x Held-out validation data: a matrix, data.frame,
#'   `SummarizedExperiment`, or list of such objects. When `NULL` (default),
#'   automatically uses held-out data stored by
#'   [optimize_panel_transferable()].
#' @param y Held-out validation labels: a factor (or list of factors when `x`
#'   is a list). When `NULL` (default), automatically uses held-out labels
#'   stored by [optimize_panel_transferable()].
#' @param on_error One of `"warn"` or `"stop"`. Controls behaviour when a
#'   solution fails evaluation. `"warn"` (default) skips the solution with a
#'   warning; `"stop"` raises an error immediately.
#' @param regularized Logical; passed to [fit_panel()]. Default `TRUE`.
#' @param verbose Logical; print progress messages. Default `interactive()`.
#' @return A data frame with one row per successfully evaluated solution and
#'   columns `solution_id`, `n_features`, `n_base_features`, `sensitivity`,
#'   `specificity`, and `auc`.
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
                                      on_error = c("warn", "stop"),
                                      regularized = TRUE,
                                      verbose = interactive()) {
  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult.", call. = FALSE)
  }
  on_error <- match.arg(on_error)

  # Resolve held-out data
  if (is.null(x) || is.null(y)) {
    ctrl <- optimization_result@control
    if (!is.null(ctrl$heldout_x) && !is.null(ctrl$heldout_y)) {
      x <- ctrl$heldout_x
      y <- ctrl$heldout_y
      if (verbose) message("Using held-out data stored in OptimizationResult.")
    } else {
      stop("Provide held-out `x` and `y`, or use a result from ",
           "optimize_panel_transferable().", call. = FALSE)
    }
  }

  n_sol <- n_solutions(optimization_result)
  if (n_sol == 0L) {
    stop("OptimizationResult contains no solutions.", call. = FALSE)
  }

  # Evaluate each solution: fit -> evaluate -> extract metrics
  rows <- vector("list", n_sol)
  for (i in seq_len(n_sol)) {
    sid <- solutions(optimization_result)$solution_id[i]
    if (verbose) message(sprintf("[%d/%d] solution %d", i, n_sol, sid))

    rows[[i]] <- tryCatch({
      panel <- fit_panel(optimization_result, solution_id = sid,
                         regularized = regularized)
      ev <- evaluate_panel(panel, x, y)

      data.frame(
        solution_id = sid,
        n_features = length(panel@features),
        n_base_features = length(panel@base_features),
        sensitivity = ev$metrics[["sensitivity"]],
        specificity = ev$metrics[["specificity"]],
        auc = ev$roc$auc,
        stringsAsFactors = FALSE
      )
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
  do.call(rbind, rows)
}
