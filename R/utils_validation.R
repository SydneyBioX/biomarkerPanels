#' Validate Binary Response
#'
#' Ensures that the response vector conforms to the required binary labels.
#'
#' @param y A factor or character vector.
#' @return Factor with levels `c("No", "Yes")`.
#' @export
ensure_binary_response <- function(y) {
  if (is.null(y)) {
    stop("`y` cannot be NULL.", call. = FALSE)
  }

  if (is.character(y)) {
    y <- factor(y, levels = c("No", "Yes"))
  }

  if (!is.factor(y) || !identical(levels(y), c("No", "Yes"))) {
    stop("`y` must be a factor with levels 'No' and 'Yes'.", call. = FALSE)
  }

  y
}
