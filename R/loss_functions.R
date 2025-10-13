#' Default loss functions for biomarker panel optimization.
#'
#' These functions compute scalar metrics such as sensitivity, specificity,
#' concordance index, and panel size. They share a common signature so that the
#' package can compose them into multi-objective optimization routines.
#'
#' Each loss accepts a binary response (`truth`), numeric prediction scores
#' (`scores`), and optionally the currently selected feature set (`selected`).
#' Additional parameters (for example, classification thresholds) can be supplied
#' via [`build_objectives()`].
#'
#' @name loss_functions
NULL

#' Sensitivity (True Positive Rate)
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param threshold Classification threshold applied to `scores`.
#' @param positive Label treated as the positive ("event") class.
#' @return Sensitivity between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_sensitivity <- function(truth, scores = NULL, selected = NULL,
                             threshold = 0.5, positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute sensitivity.", call. = FALSE)
  }
  positives <- truth == positive
  predicted <- scores >= threshold
  tp <- sum(positives & predicted)
  fn <- sum(positives & !predicted)
  if ((tp + fn) == 0) {
    return(NA_real_)
  }
  tp / (tp + fn)
}

#' Specificity (True Negative Rate)
#'
#' @inheritParams loss_sensitivity
#' @return Specificity between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_specificity <- function(truth, scores = NULL, selected = NULL,
                             threshold = 0.5, positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute specificity.", call. = FALSE)
  }
  negatives <- truth != positive
  predicted <- scores >= threshold
  tn <- sum(negatives & !predicted)
  fp <- sum(negatives & predicted)
  if ((tn + fp) == 0) {
    return(NA_real_)
  }
  tn / (tn + fp)
}

#' Concordance Index (c-index)
#'
#' Computes the probability that a random positive instance receives a higher
#' score than a random negative instance. Equivalent to the AUC for binary
#' outcomes.
#'
#' @inheritParams loss_sensitivity
#' @return Concordance index between 0 and 1, or `NA_real_` if undefined.
#' @export
loss_c_index <- function(truth, scores = NULL, selected = NULL,
                         positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute c-index.", call. = FALSE)
  }
  pos_scores <- scores[truth == positive]
  neg_scores <- scores[truth != positive]
  if (length(pos_scores) == 0L || length(neg_scores) == 0L) {
    return(NA_real_)
  }
  pair_matrix <- outer(pos_scores, neg_scores, FUN = "-")
  concordant <- sum(pair_matrix > 0)
  ties <- sum(pair_matrix == 0)
  total <- length(pair_matrix)
  (concordant + 0.5 * ties) / total
}

#' Panel Size Penalty
#'
#' @inheritParams loss_sensitivity
#' @param selected Vector of selected features (character, numeric, or logical).
#' @return Count of selected biomarkers.
#' @export
loss_num_features <- function(truth = NULL, scores = NULL, selected = NULL, ...) {
  if (is.null(selected)) {
    return(0)
  }
  if (is.logical(selected)) {
    return(sum(selected))
  }
  length(selected)
}

#' Balanced Accuracy
#'
#' @inheritParams loss_sensitivity
#' @return Mean of sensitivity and specificity.
#' @export
loss_balanced_accuracy <- function(truth, scores = NULL, selected = NULL,
                                   threshold = 0.5, positive = "Yes") {
  sens <- loss_sensitivity(truth, scores, selected, threshold = threshold, positive = positive)
  spec <- loss_specificity(truth, scores, selected, threshold = threshold, positive = positive)
  if (is.na(sens) || is.na(spec)) {
    return(NA_real_)
  }
  (sens + spec) / 2
}

#' Minimum Cohort Sensitivity
#'
#' Computes sensitivity within each cohort and returns the minimum value to
#' capture worst-case performance.
#'
#' @inheritParams loss_sensitivity
#' @param cohort Factor indicating cohort membership.
#' @return Sensitivity of the weakest cohort.
#' @export
loss_min_cohort_sensitivity <- function(truth, scores = NULL, selected = NULL,
                                        threshold = 0.5, positive = "Yes",
                                        cohort = NULL) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute sensitivity.", call. = FALSE)
  }
  if (is.null(cohort)) {
    return(loss_sensitivity(truth, scores, selected,
                            threshold = threshold, positive = positive))
  }
  cohort <- factor(cohort)
  if (length(cohort) != length(truth)) {
    stop("Length of `cohort` must match `truth`.", call. = FALSE)
  }
  values <- vapply(levels(cohort), function(level) {
    idx <- !is.na(cohort) & cohort == level
    if (!any(idx)) {
      return(NA_real_)
    }
    loss_sensitivity(
      truth[idx],
      scores[idx],
      selected = selected,
      threshold = threshold,
      positive = positive
    )
  }, numeric(1))
  if (all(is.na(values))) {
    return(NA_real_)
  }
  min(values, na.rm = TRUE)
}

#' Minimum Cohort Specificity
#'
#' Computes specificity within each cohort and returns the minimum value.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Specificity of the weakest cohort.
#' @export
loss_min_cohort_specificity <- function(truth, scores = NULL, selected = NULL,
                                        threshold = 0.5, positive = "Yes",
                                        cohort = NULL) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute specificity.", call. = FALSE)
  }
  if (is.null(cohort)) {
    return(loss_specificity(truth, scores, selected,
                            threshold = threshold, positive = positive))
  }
  cohort <- factor(cohort)
  if (length(cohort) != length(truth)) {
    stop("Length of `cohort` must match `truth`.", call. = FALSE)
  }
  values <- vapply(levels(cohort), function(level) {
    idx <- !is.na(cohort) & cohort == level
    if (!any(idx)) {
      return(NA_real_)
    }
    loss_specificity(
      truth[idx],
      scores[idx],
      selected = selected,
      threshold = threshold,
      positive = positive
    )
  }, numeric(1))
  if (all(is.na(values))) {
    return(NA_real_)
  }
  min(values, na.rm = TRUE)
}

#' Cohort Sensitivity Range
#'
#' Difference between maximum and minimum cohort sensitivities. Smaller values
#' indicate more uniform transfer across cohorts.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Sensitivity range across cohorts.
#' @export
loss_cohort_sensitivity_gap <- function(truth, scores = NULL, selected = NULL,
                                        threshold = 0.5, positive = "Yes",
                                        cohort = NULL) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute sensitivity.", call. = FALSE)
  }
  if (is.null(cohort)) {
    sens <- loss_sensitivity(truth, scores, selected,
                             threshold = threshold, positive = positive)
    return(0 * sens)
  }
  cohort <- factor(cohort)
  if (length(cohort) != length(truth)) {
    stop("Length of `cohort` must match `truth`.", call. = FALSE)
  }
  values <- vapply(levels(cohort), function(level) {
    idx <- !is.na(cohort) & cohort == level
    if (!any(idx)) {
      return(NA_real_)
    }
    loss_sensitivity(
      truth[idx],
      scores[idx],
      selected = selected,
      threshold = threshold,
      positive = positive
    )
  }, numeric(1))
  if (all(is.na(values))) {
    return(NA_real_)
  }
  max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
}

#' Maximum Cohort Brier Score
#'
#' Computes the Brier score (mean squared error on probabilities) within each
#' cohort and returns the maximum, highlighting the worst calibrated cohort.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @return Maximum Brier score across cohorts.
#' @export
loss_max_cohort_brier <- function(truth, scores = NULL, selected = NULL,
                                  positive = "Yes", cohort = NULL) {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute the Brier score.", call. = FALSE)
  }
  target <- as.numeric(truth == positive)
  cohort_vec <- factor(if (is.null(cohort)) rep("cohort_01", length(truth)) else cohort)
  if (length(cohort_vec) != length(truth)) {
    stop("Length of `cohort` must match `truth`.", call. = FALSE)
  }
  values <- vapply(levels(cohort_vec), function(level) {
    idx <- !is.na(cohort_vec) & cohort_vec == level
    if (!any(idx)) {
      return(NA_real_)
    }
    mean((scores[idx] - target[idx])^2)
  }, numeric(1))
  max(values, na.rm = TRUE)
}

#' Maximum Mean Shift Across Cohorts
#'
#' Computes pairwise distances between cohort-specific mean expression vectors
#' for the selected features and returns the maximum distance. Encourages panels
#' whose selected biomarkers exhibit similar distributions across cohorts.
#'
#' @inheritParams loss_min_cohort_sensitivity
#' @param x Matrix of selected feature values (samples x features).
#' @return Maximum pairwise distance between cohort means.
#' @export
loss_max_cohort_mean_shift <- function(truth, scores = NULL, selected = NULL,
                                       cohort = NULL, x = NULL) {
  if (is.null(cohort) || is.null(x) || ncol(x) == 0L) {
    return(0)
  }
  cohort <- factor(cohort)
  if (length(cohort) != nrow(x)) {
    stop("Length of `cohort` must match the number of rows in `x`.", call. = FALSE)
  }
  if (nlevels(cohort) <= 1L) {
    return(0)
  }
  means <- lapply(levels(cohort), function(level) {
    idx <- !is.na(cohort) & cohort == level
    if (!any(idx)) {
      return(NULL)
    }
    colMeans(x[idx, , drop = FALSE])
  })
  means <- Filter(Negate(is.null), means)
  if (length(means) <= 1L) {
    return(0)
  }
  mean_matrix <- do.call(rbind, means)
  if (nrow(mean_matrix) <= 1L) {
    return(0)
  }
  max(stats::dist(mean_matrix))
}

.loss_registry <- new.env(parent = emptyenv())

.register_default_loss <- function(name, fun, direction, label) {
  assign(name, list(fun = fun, direction = direction, label = label),
         envir = .loss_registry)
}

.register_default_loss(
  "sensitivity", loss_sensitivity, "maximize", "Sensitivity"
)
.register_default_loss(
  "specificity", loss_specificity, "maximize", "Specificity"
)
.register_default_loss(
  "c_index", loss_c_index, "maximize", "Concordance Index"
)
.register_default_loss(
  "num_features", loss_num_features, "minimize", "Number of Features"
)
.register_default_loss(
  "balanced_accuracy", loss_balanced_accuracy, "maximize", "Balanced Accuracy"
)
.register_default_loss(
  "min_cohort_sensitivity", loss_min_cohort_sensitivity, "maximize",
  "Minimum Cohort Sensitivity"
)
.register_default_loss(
  "min_cohort_specificity", loss_min_cohort_specificity, "maximize",
  "Minimum Cohort Specificity"
)
.register_default_loss(
  "cohort_sensitivity_gap", loss_cohort_sensitivity_gap, "minimize",
  "Cohort Sensitivity Gap"
)
.register_default_loss(
  "max_cohort_brier", loss_max_cohort_brier, "minimize",
  "Maximum Cohort Brier Score"
)
.register_default_loss(
  "max_cohort_mean_shift", loss_max_cohort_mean_shift, "minimize",
  "Maximum Cohort Mean Shift"
)

#' Register a loss function.
#'
#' @param name Unique identifier for the loss.
#' @param fun Function implementing the loss. Must accept at least the arguments
#'   `truth`, `scores`, and `selected`.
#' @param direction Either `"maximize"` or `"minimize"`.
#' @param label Human-readable name.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#' @export
register_loss_function <- function(name, fun,
                                   direction = c("maximize", "minimize"),
                                   label = name, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  direction <- match.arg(direction)
  if (!overwrite && exists(name, envir = .loss_registry, inherits = FALSE)) {
    stop(sprintf("Loss function '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(fun = fun, direction = direction, label = label),
         envir = .loss_registry)
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
#' supplied per loss via the `params` list (e.g., custom thresholds).
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
