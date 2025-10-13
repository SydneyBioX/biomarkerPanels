#' Define optimization objectives.
#'
#' Assemble a list of objective descriptors ready for consumption by the
#' optimization engine. Built-in losses are pulled from the registry created in
#' [`build_objectives()`]; custom entries may be appended for experimental goals.
#'
#' @param losses Character vector of registered loss names (default:
#'   `c("sensitivity", "specificity")`).
#' @param custom Optional named list of objective descriptors with elements
#'   `label`, `direction`, and `fun`.
#' @param params Named list of additional argument lists for specific losses.
#' @param directions Optional named vector overriding loss directions.
#' @return Named list of objective descriptors.
#' @export
define_objectives <- function(losses = c("sensitivity", "specificity"),
                              custom = NULL,
                              params = list(),
                              directions = NULL) {
  objectives <- build_objectives(losses, params = params, directions = directions)

  if (!is.null(custom)) {
    stopifnot(is.list(custom))
    valid_custom <- vapply(custom, function(entry) {
      is.list(entry) &&
        !is.null(entry$label) &&
        !is.null(entry$direction) &&
        is.function(entry$fun)
    }, logical(1))
    if (!all(valid_custom)) {
      stop("All custom objectives must supply label, direction, and fun.", call. = FALSE)
    }
    objectives <- c(objectives, custom)
  }

  objectives
}
