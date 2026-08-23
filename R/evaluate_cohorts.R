#' Per-Cohort Panel Evaluation
#'
#' Functions for evaluating biomarker panel performance on a per-cohort basis,
#' enabling analysis of transferability across independent datasets.
#'
#' @name evaluate_cohorts
NULL

#' Evaluate a Biomarker Panel by Cohort
#'
#' Compute per-cohort performance metrics for a fitted panel on validation data.
#' Unlike [evaluate_panel()] which returns aggregate metrics, this function
#' returns a data frame with one row per cohort, enabling comparison of panel
#' transferability across independent datasets.
#'
#' @param panel A `BiomarkerPanelResult`.
#' @param x Validation feature matrix, or a list of matrices for multi-cohort
#'   evaluation. Names of the list become cohort identifiers.
#' @param y Validation binary response factor, or a list of factors when `x` is
#'   a list.
#' @param metrics Character vector of registered metric names to compute (any
#'   name in [metric_registry()], e.g. `"pauc"` or
#'   `"specificity_at_sensitivity"`). Defaults to
#'   `c("sensitivity", "specificity", "auc")`.
#' @param cutoff_prob Classification probability threshold. Defaults to `0.5`.
#' @param positive Label treated as the positive class. Defaults to the value
#'   stored in the panel, or `"Yes"`.
#' @param feature_transform Optional feature transform name. Defaults to the
#'   transform stored in the fitted panel.
#' @param assay For `SummarizedExperiment` inputs, assay name or index.
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{cohort}{Character; cohort identifier}
#'     \item{n_samples}{Integer; number of samples in cohort}
#'     \item{n_positive}{Integer; number of positive samples}
#'     \item{n_negative}{Integer; number of negative samples}
#'     \item{...}{One column per metric in `metrics` argument}
#'   }
#'   Attributes: `cutoff_prob`, `positive`, `features`.
#' @export
#' @seealso [evaluate_panel()] for aggregate evaluation,
#'   [plot_cohort_comparison()] for visualization.
#' @examples
#' # With a fitted panel and multi-cohort validation data:
#' # cohort_metrics <- evaluate_panel_by_cohort(panel, x_list, y_list)
#' # plot_cohort_comparison(cohort_metrics)
evaluate_panel_by_cohort <- function(panel,
                                     x,
                                     y,
                                     metrics = c("sensitivity", "specificity", "auc"),
                                     cutoff_prob = 0.5,
                                     positive = NULL,
                                     feature_transform = NULL,
                                     assay = NULL) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  # Get base and transformed features from panel
  base_features <- panel@base_features
  features <- panel@features
  if (is.null(base_features) || length(base_features) == 0L) {
    # Backward compatibility
    base_features <- features
  }
  if (length(features) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  stored_model <- panel@model

  if (is.null(positive)) {
    stored_positive <- panel@control$positive_class
    positive <- if (!is.null(stored_positive)) stored_positive else "Yes"
  }

  stored_transform <- panel@control$feature_transform
  if (is.null(feature_transform)) {
    feature_transform <- if (is.null(stored_transform)) "none" else stored_transform
  }

  metric_objs <- build_objectives(metrics)

  # Normalize single-cohort input to list format

  if (!is.list(x) || inherits(x, "SummarizedExperiment")) {
    x <- list(x)
    y <- list(y)
  }

  if (length(x) != length(y)) {
    stop("`x` and `y` lists must have the same length.", call. = FALSE)
  }

  cohort_names <- .default_cohort_names(x)

  # Process each cohort independently
  cohort_results <- lapply(seq_along(x), function(i) {
    xi <- .extract_feature_matrix(x[[i]], assay = assay)
    yi <- ensure_binary_response(y[[i]], positive = positive)

    if (nrow(xi) != length(yi)) {
      stop(sprintf("Cohort '%s': `x` and `y` must have matching sample sizes.",
                   cohort_names[i]), call. = FALSE)
    }

    xi <- .ensure_feature_colnames(xi)

    # Validate base features and apply the panel's transform on the fly.
    x_selected <- .prepare_scoring_matrix(
      xi, base_features, feature_transform,
      context = sprintf("cohort '%s'", cohort_names[i])
    )$x_selected

    # CPOP-style panels store a subset of all possible pair features; subset
    # the transformed matrix to match the model's expected columns.
    if (feature_transform != "none" && length(base_features) >= 2L &&
        length(features) && all(features %in% colnames(x_selected))) {
      x_selected <- x_selected[, features, drop = FALSE]
    }

    # Generate scores using stored model
    cohort_vec <- factor(rep(cohort_names[i], nrow(x_selected)))
    if (!is.null(stored_model)) {
      scores <- tryCatch({
        .predict_panel_model(
          stored_model, x_selected,
          cohort = cohort_vec, expected_features = features
        )
      }, error = function(e) {
        stop(sprintf("Cohort '%s': Failed to generate predictions: %s",
                     cohort_names[i], conditionMessage(e)), call. = FALSE)
      })
    } else {
      stop(
        "Cohort '", cohort_names[i], "': panel has no fitted model. ",
        "Use fit_panel() to fit a model before evaluate_panel_by_cohort() ",
        "(refitting on the evaluation data would report leaky in-sample metrics).",
        call. = FALSE
      )
    }

    n_pos <- sum(yi == positive)
    n_neg <- sum(yi != positive)

    # Any registered metric works here; the wrapper drops arguments a metric
    # does not accept (e.g. cutoff_prob for the threshold-free metrics).
    metric_values <- vapply(metrics, function(m) {
      metric_objs[[m]]$fun(yi, scores,
                           cutoff_prob = cutoff_prob, positive = positive)
    }, numeric(1))

    c(
      list(
        cohort = cohort_names[i],
        n_samples = length(yi),
        n_positive = n_pos,
        n_negative = n_neg
      ),
      as.list(metric_values)
    )
  })

  result <- do.call(rbind.data.frame, lapply(cohort_results, as.data.frame,
                                              stringsAsFactors = FALSE))
  rownames(result) <- NULL


  attr(result, "cutoff_prob") <- cutoff_prob
  attr(result, "positive") <- positive
  attr(result, "features") <- features

  result
}
