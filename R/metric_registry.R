#' Metric Function Registry
#'
#' Infrastructure for registering, retrieving, and composing metric functions
#' into multi-objective optimization targets.
#'
#' @name metric_registry
NULL

.metric_registry <- new.env(parent = emptyenv())

#' Register a Default Metric Function
#'
#' Internal helper to register built-in metric functions in the metric registry.
#' Used during package loading to populate `.metric_registry`.
#'
#' @param name Character name for the metric function.
#' @param fun The metric function itself.
#' @param direction Either "maximize" or "minimize".
#' @param label Human-readable label for the metric.
#' @param cutoff_dependent Logical; whether the metric requires a probability
#'   cutoff for classification.
#' @return Invisible NULL; called for side effect of registering the metric.
#' @keywords internal
.register_default_metric <- function(name, fun, direction, label,
                                      cutoff_dependent = FALSE) {
  assign(name, list(
    fun = fun,
    direction = direction,
    label = label,
    cutoff_dependent = cutoff_dependent
  ), envir = .metric_registry)
}

# -----------------------------------------------------------------------------
# Register Default Loss Functions
# -----------------------------------------------------------------------------

# Cutoff-dependent metrics (require probability threshold for classification)
.register_default_metric(
  "sensitivity", metric_sensitivity, "maximize", "Sensitivity",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "specificity", metric_specificity, "maximize", "Specificity",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "f1", metric_f1, "maximize", "F1 Score",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "precision", metric_precision, "maximize", "Precision",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "npv", metric_npv, "maximize", "Negative Predictive Value",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "balanced_accuracy", metric_balanced_accuracy, "maximize", "Balanced Accuracy",
  cutoff_dependent = TRUE
)
.register_default_metric(
  "min_cohort_auc", metric_min_cohort_auc, "maximize",
  "Minimum Cohort AUC"
)
.register_default_metric(
  "cohort_auc_gap", metric_cohort_auc_gap, "minimize",
  "Cohort AUC Gap"
)
.register_default_metric(
  "cohort_auc_var", metric_cohort_auc_var, "minimize",
  "Cohort AUC Variance"
)

# Cutoff-free metrics (do not depend on probability threshold)
.register_default_metric(
  "auc", metric_auc, "maximize", "Area Under ROC Curve"
)
.register_default_metric(
  "num_features", metric_num_features, "minimize", "Number of Features"
)
.register_default_metric(
  "max_cohort_brier", metric_max_cohort_brier, "minimize",
  "Maximum Cohort Brier Score"
)
.register_default_metric(
  "pauc", metric_pauc, "maximize", "Partial AUC (High Sensitivity)"
)
.register_default_metric(
  "specificity_at_sensitivity", metric_specificity_at_sensitivity, "maximize",
  "Specificity at Target Sensitivity"
)
.register_default_metric(
  "cohort_leakage", metric_cohort_leakage, "minimize",
  "Cohort Leakage (Adjusted R-squared)"
)
.register_default_metric(
  "conditional_score_shift", metric_conditional_score_shift, "minimize",
  "Class-Conditional Score Shift (W1)"
)
.register_default_metric(
  "easy_hard_accuracy", metric_easy_hard_accuracy, "maximize",
  "Easy-Hard Balanced Accuracy",
  cutoff_dependent = TRUE
)

#' Register a metric function.
#'
#' @param name Unique identifier for the metric.
#' @param fun Function implementing the metric. Must accept at least the arguments
#'   `truth`, `scores`, and `selected`.
#' @param direction Either `"maximize"` or `"minimize"`.
#' @param label Human-readable name.
#' @param cutoff_dependent Logical; if `TRUE`, the metric depends on a probability
#'   cutoff for classification (e.g., sensitivity, specificity). Cutoff-dependent
#'   metrics may not behave as expected in [min_metric_constraint()] because
#'   their values depend on the arbitrary `cutoff_prob` parameter.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#' @export
register_metric <- function(name, fun,
                             direction = c("maximize", "minimize"),
                             label = name,
                             cutoff_dependent = FALSE,
                             overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  direction <- match.arg(direction)
  if (!overwrite && exists(name, envir = .metric_registry, inherits = FALSE)) {
    stop(sprintf("Metric '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(
    fun = fun,
    direction = direction,
    label = label,
    cutoff_dependent = cutoff_dependent
  ), envir = .metric_registry)
  invisible(name)
}

#' List registered metric functions.
#'
#' @return Named list of metric registrations (`fun`, `direction`, `label`).
#' @export
metric_registry <- function() {
  if (length(ls(.metric_registry)) == 0L) {
    return(list())
  }
  mget(ls(.metric_registry), envir = .metric_registry)
}

#' Build objective descriptors from registered metrics.
#'
#' Creates a list compatible with [define_objectives()] where each entry contains
#' the metric label, optimization direction, and an evaluation function with the
#' signature `(truth, estimate, selected = NULL)`. Additional parameters may be
#' supplied per metric via the `params` list (e.g., custom cutoffs).
#'
#' @param metrics Character vector of registered metric names.
#' @param params Named list mapping metric names to argument lists.
#' @param directions Optional named character vector overriding the optimization
#'   direction for subsets of `metrics`.
#' @return Named list of objective descriptors.
#' @export
build_objectives <- function(metrics,
                             params = list(),
                             directions = NULL) {
  stopifnot(is.character(metrics), length(metrics) >= 1L)
  registry <- metric_registry()
  missing <- setdiff(metrics, names(registry))
  if (length(missing)) {
    stop(sprintf("Metric(s) not registered: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  objs <- lapply(metrics, function(name) {
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
  names(objs) <- metrics
  objs
}
