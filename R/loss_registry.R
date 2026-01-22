#' Loss Function Registry
#'
#' Infrastructure for registering, retrieving, and composing loss functions
#' into multi-objective optimization targets.
#'
#' @name loss_registry
NULL

.loss_registry <- new.env(parent = emptyenv())

#' Register a Default Loss Function
#'
#' Internal helper to register built-in loss functions in the loss registry.
#' Used during package loading to populate `.loss_registry`.
#'
#' @param name Character name for the loss function.
#' @param fun The loss function itself.
#' @param direction Either "maximize" or "minimize".
#' @param label Human-readable label for the metric.
#' @param cutoff_dependent Logical; whether the metric requires a probability
#'   cutoff for classification.
#' @return Invisible NULL; called for side effect of registering the loss.
#' @keywords internal
.register_default_loss <- function(name, fun, direction, label,
                                    cutoff_dependent = FALSE) {
  assign(name, list(
    fun = fun,
    direction = direction,
    label = label,
    cutoff_dependent = cutoff_dependent
  ), envir = .loss_registry)
}

# -----------------------------------------------------------------------------
# Register Default Loss Functions
# -----------------------------------------------------------------------------

# Cutoff-dependent metrics (require probability threshold for classification)
.register_default_loss(
  "sensitivity", loss_sensitivity, "maximize", "Sensitivity",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "specificity", loss_specificity, "maximize", "Specificity",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "f1", loss_f1, "maximize", "F1 Score",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "precision", loss_precision, "maximize", "Precision",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "npv", loss_npv, "maximize", "Negative Predictive Value",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "balanced_accuracy", loss_balanced_accuracy, "maximize", "Balanced Accuracy",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "min_cohort_sensitivity", loss_min_cohort_sensitivity, "maximize",
  "Minimum Cohort Sensitivity",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "min_cohort_specificity", loss_min_cohort_specificity, "maximize",
  "Minimum Cohort Specificity",
  cutoff_dependent = TRUE
)
.register_default_loss(
  "cohort_sensitivity_gap", loss_cohort_sensitivity_gap, "minimize",
  "Cohort Sensitivity Gap",
  cutoff_dependent = TRUE
)

# Cutoff-free metrics (do not depend on probability threshold)
.register_default_loss(
  "auc", loss_auc, "maximize", "Area Under ROC Curve"
)
.register_default_loss(
  "num_features", loss_num_features, "minimize", "Number of Features"
)
.register_default_loss(
  "max_cohort_brier", loss_max_cohort_brier, "minimize",
  "Maximum Cohort Brier Score"
)
.register_default_loss(
  "max_cohort_mean_shift", loss_max_cohort_mean_shift, "minimize",
  "Maximum Cohort Mean Shift"
)
.register_default_loss(
  "pauc", loss_pauc, "maximize", "Partial AUC (High Sensitivity)"
)
.register_default_loss(
  "specificity_at_sensitivity", loss_specificity_at_sensitivity, "maximize",
  "Specificity at Target Sensitivity"
)

#' Register a loss function.
#'
#' @param name Unique identifier for the loss.
#' @param fun Function implementing the loss. Must accept at least the arguments
#'   `truth`, `scores`, and `selected`.
#' @param direction Either `"maximize"` or `"minimize"`.
#' @param label Human-readable name.
#' @param cutoff_dependent Logical; if `TRUE`, the loss depends on a probability
#'   cutoff for classification (e.g., sensitivity, specificity). Cutoff-dependent
#'   losses may not behave as expected in [min_metric_constraint()] because
#'   their values depend on the arbitrary `cutoff_prob` parameter.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#' @export
register_loss_function <- function(name, fun,
                                   direction = c("maximize", "minimize"),
                                   label = name,
                                   cutoff_dependent = FALSE,
                                   overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  direction <- match.arg(direction)
  if (!overwrite && exists(name, envir = .loss_registry, inherits = FALSE)) {
    stop(sprintf("Loss function '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(
    fun = fun,
    direction = direction,
    label = label,
    cutoff_dependent = cutoff_dependent
  ), envir = .loss_registry)
  invisible(name)
}

#' List registered loss functions.
#'
#' @return Named list of loss registrations (`fun`, `direction`, `label`).
#' @export
loss_registry <- function() {
  if (length(ls(.loss_registry)) == 0L) {
    return(list())
  }
  mget(ls(.loss_registry), envir = .loss_registry)
}

#' Build objective descriptors from registered losses.
#'
#' Creates a list compatible with [define_objectives()] where each entry contains
#' the loss label, optimization direction, and an evaluation function with the
#' signature `(truth, estimate, selected = NULL)`. Additional parameters may be
#' supplied per loss via the `params` list (e.g., custom cutoffs).
#'
#' @param losses Character vector of registered loss names.
#' @param params Named list mapping loss names to argument lists.
#' @param directions Optional named character vector overriding the optimization
#'   direction for subsets of `losses`.
#' @return Named list of objective descriptors.
#' @export
build_objectives <- function(losses,
                             params = list(),
                             directions = NULL) {
  stopifnot(is.character(losses), length(losses) >= 1L)
  registry <- loss_registry()
  missing <- setdiff(losses, names(registry))
  if (length(missing)) {
    stop(sprintf("Loss function(s) not registered: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  objs <- lapply(losses, function(name) {
    entry <- registry[[name]]
    extras <- params[[name]]
    if (is.null(extras)) {
      extras <- list()
    } else if (!is.list(extras)) {
      extras <- as.list(extras)
    }
    direction <- entry$direction
    if (!is.null(directions) && !is.null(directions[[name]])) {
      direction <- match.arg(directions[[name]], c("maximize", "minimize"))
    }
    fun <- entry$fun
    fun_formals <- names(formals(fun))
    has_dots <- any(fun_formals == "...")
    wrapper <- function(truth, estimate, selected = NULL, ...) {
      dots <- list(...)
      extras_filtered <- extras
      if (!has_dots) {
        if (length(extras_filtered)) {
          extras_filtered <- extras_filtered[names(extras_filtered) %in% fun_formals]
        }
        if (length(dots)) {
          dots <- dots[names(dots) %in% fun_formals]
        }
      }
      args <- c(
        list(truth = truth, scores = estimate, selected = selected),
        extras_filtered,
        dots
      )
      do.call(fun, args)
    }
    list(
      label = entry$label,
      direction = direction,
      fun = wrapper
    )
  })
  names(objs) <- losses
  objs
}
