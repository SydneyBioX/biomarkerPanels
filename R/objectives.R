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
#' @param params Named list of additional argument lists for specific losses
#'   (e.g., custom cutoff probabilities for sensitivity).
#' @param directions Optional named vector overriding loss directions.
#' @param cutoff_strategy Strategy for computing classification cutoff. One of:
#'   \describe{
#'     \item{"fixed"}{Use `cutoff_prob` parameter (default 0.5)}
#'     \item{"prevalence"}{Cutoff equals training class prevalence}
#'     \item{"youden"}{Optimal Youden's J point from ROC curve}
#'   }
#' @param cutoff_prob Fixed cutoff probability (default 0.5). Only used when
#'   `cutoff_strategy = "fixed"`.
#' @return Named list of objective descriptors.
#' @export
define_objectives <- function(losses = c("sensitivity", "specificity"),
                              custom = NULL,
                              params = list(),
                              directions = NULL,
                              cutoff_strategy = c("fixed", "prevalence", "youden"),
                              cutoff_prob = 0.5) {
  cutoff_strategy <- match.arg(cutoff_strategy)

  # Add cutoff_strategy to params for losses that use cutoffs
  cutoff_losses <- c("sensitivity", "specificity", "balanced_accuracy",
                     "min_cohort_sensitivity", "min_cohort_specificity",
                     "cohort_sensitivity_gap")

  for (loss_name in intersect(losses, cutoff_losses)) {
    if (is.null(params[[loss_name]])) {
      params[[loss_name]] <- list()
    }
    # Only set if not already specified
    if (is.null(params[[loss_name]]$cutoff_strategy)) {
      params[[loss_name]]$cutoff_strategy <- cutoff_strategy
    }
    if (is.null(params[[loss_name]]$cutoff_prob)) {
      params[[loss_name]]$cutoff_prob <- cutoff_prob
    }
  }

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

#' Minimum metric constraint constructor.
#'
#' Create a boolean constraint that requires a registered loss to meet or exceed
#' (or undercut, for minimised losses) a specified threshold during optimisation.
#' These constraints can be supplied to [optimize_panel()] via its `constraints`
#' argument.
#'
#' @param loss Character scalar naming a registered loss function.
#' @param threshold Numeric scalar describing the required metric level.
#' @param params Optional named list of additional parameters forwarded to the
#'   loss function (e.g., cutoff probabilities for sensitivity).
#' @param label Optional human-readable label for the constraint; defaults to
#'   `paste0("min_", loss, "_", threshold)`.
#' @return A constraint descriptor (list) with elements `label`, `fun`,
#'   `threshold`, `loss`, and `direction`.
#' @export
min_metric_constraint <- function(loss,
                                  threshold,
                                  params = list(),
                                  label = NULL) {
  stopifnot(is.character(loss), length(loss) == 1L, nzchar(loss))
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("`threshold` must be a finite numeric scalar.", call. = FALSE)
  }

  # Check if loss is cutoff-dependent and warn

  registry <- loss_registry()
  if (!loss %in% names(registry)) {
    stop(sprintf("Loss '%s' is not registered.", loss), call. = FALSE)
  }
  if (isTRUE(registry[[loss]]$cutoff_dependent)) {
    warning(
      sprintf("Loss '%s' depends on a probability cutoff. ", loss),
      "Constraint behavior depends on the cutoff_prob parameter (default 0.5), ",
      "which may not reflect the ROC trade-off you want. Consider using ",
      "cutoff-free metrics like 'auc', 'pauc', or 'specificity_at_sensitivity' ",
      "for constraints.",
      call. = FALSE
    )
  }

  if (!is.list(params)) {
    params <- as.list(params)
  }
  param_list <- list()
  if (length(params)) {
    param_list[[loss]] <- params
  }
  entry <- build_objectives(losses = loss, params = param_list)[[loss]]
  direction <- entry$direction

  if (is.null(label) || !nzchar(label)) {
    label <- sprintf(
      "min_%s_%s",
      loss,
      format(threshold, trim = TRUE, scientific = FALSE)
    )
  }

  fun <- function(truth, scores, selected = NULL, cohort = NULL, x = NULL, ...) {
    value <- entry$fun(
      truth = truth,
      estimate = scores,
      selected = selected,
      cohort = cohort,
      x = x,
      ...
    )
    if (is.na(value)) {
      return(FALSE)
    }
    if (direction == "maximize") {
      value >= threshold
    } else {
      value <= threshold
    }
  }

  list(
    label = label,
    fun = fun,
    threshold = threshold,
    loss = loss,
    direction = direction
  )
}