#' Define Optimization Objectives
#'
#' Create a standardized list of objectives used by the optimization routine.
#' Each objective includes a label, the direction (`"maximize"`/`"minimize"`),
#' and an evaluation function that consumes predictions and true labels.
#'
#' @param objectives Named list or character vector of built-in objectives.
#' @param custom Named list of user-supplied functions returning scalar scores.
#' @return A list of objective descriptors consumed by `optimize_panel()`.
#' @export
define_objectives <- function(objectives = c("sensitivity", "specificity"),
                              custom = NULL) {
  base <- list()

  if ("sensitivity" %in% objectives) {
    base$sensitivity <- list(
      label = "Sensitivity",
      direction = "maximize",
      fun = function(truth, estimate) {
        tp <- sum(truth & estimate)
        fn <- sum(truth & !estimate)
        if ((tp + fn) == 0) return(NA_real_)
        tp / (tp + fn)
      }
    )
  }

  if ("specificity" %in% objectives) {
    base$specificity <- list(
      label = "Specificity",
      direction = "maximize",
      fun = function(truth, estimate) {
        tn <- sum(!truth & !estimate)
        fp <- sum(!truth & estimate)
        if ((tn + fp) == 0) return(NA_real_)
        tn / (tn + fp)
      }
    )
  }

  if (!is.null(custom)) {
    stopifnot(is.list(custom))
    base <- c(base, custom)
  }

  base
}
