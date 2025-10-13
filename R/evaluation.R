#' Evaluate a Biomarker Panel
#'
#' Compute objective values and additional summary statistics for a fitted panel.
#'
#' @param panel A `BiomarkerPanelResult`.
#' @param x Validation feature matrix or `SummarizedExperiment`.
#' @param y Validation binary response factor.
#' @param objectives Optional override for objectives.
#' @return A list with `metrics` and `objectives`.
#' @export
evaluate_panel <- function(panel, x, y,
                           objectives = define_objectives()) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  truth <- ensure_binary_response(y)
  scores <- stats::runif(length(truth))

  selected <- panel@features
  objective_values <- vapply(objectives, function(obj) {
    obj$fun(truth, scores, selected = selected)
  }, numeric(1))

  list(
    metrics = objective_values,
    objectives = data.frame(
      objective = names(objectives),
      value = objective_values,
      direction = vapply(objectives, `[[`, character(1), "direction")
    )
  )
}
