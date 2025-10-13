#' Plot Pareto Fronts
#'
#' Visualize trade-offs between two objectives for candidate biomarker panels.
#'
#' @param results Data frame with columns `objective_x`, `objective_y`, `label`.
#' @param xlab,ylab Axis labels.
#' @return A `ggplot` object.
#' @export
plot_pareto_front <- function(results,
                              xlab = "Objective X",
                              ylab = "Objective Y") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_pareto_front()", call. = FALSE)
  }

  ggplot2::ggplot(results, ggplot2::aes(x = objective_x, y = objective_y, label = label)) +
    ggplot2::geom_point(color = "#005B96") +
    ggplot2::geom_text(vjust = -0.4, size = 3) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12)
}
