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
#' @param metrics Character vector of metric names to compute. Defaults to
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

  # Normalize single-cohort input to list format

  if (!is.list(x) || inherits(x, "SummarizedExperiment")) {
    x <- list(x)
    y <- list(y)
  }

  if (length(x) != length(y)) {
    stop("`x` and `y` lists must have the same length.", call. = FALSE)
  }

  cohort_names <- names(x)
  if (is.null(cohort_names) || any(cohort_names == "")) {
    cohort_names <- sprintf("cohort_%02d", seq_along(x))
  }

  # Process each cohort independently
  cohort_results <- lapply(seq_along(x), function(i) {
    xi <- .extract_feature_matrix(x[[i]], assay = assay)
    yi <- ensure_binary_response(y[[i]], positive = positive)

    if (nrow(xi) != length(yi)) {
      stop(sprintf("Cohort '%s': `x` and `y` must have matching sample sizes.",
                   cohort_names[i]), call. = FALSE)
    }

    if (is.null(colnames(xi))) {
      colnames(xi) <- sprintf("feature_%04d", seq_len(ncol(xi)))
    }

    # Validate base features are present
    if (!all(base_features %in% colnames(xi))) {
      missing <- setdiff(base_features, colnames(xi))
      stop(sprintf("Cohort '%s': Base feature(s) not found: %s",
                   cohort_names[i], paste(missing, collapse = ", ")),
           call. = FALSE)
    }

    # Extract base features and apply transform
    x_base <- xi[, base_features, drop = FALSE]
    if (feature_transform != "none" && length(base_features) >= 2L) {
      x_selected <- .apply_feature_transform_single(x_base, feature_transform)
    } else {
      x_selected <- x_base
    }

    # Generate scores using stored model
    cohort_vec <- factor(rep(cohort_names[i], nrow(x_selected)))
    if (!is.null(stored_model)) {
      scores <- tryCatch({
        if (inherits(stored_model, "cv.glmnet")) {
          .predict_glmnet_model(stored_model, x_selected, cohort_vec)
        } else {
          newdata <- as.data.frame(x_selected, check.names = TRUE)
          expected_names <- make.names(features, unique = TRUE)
          names(newdata) <- expected_names

          model_cols <- names(stored_model$model)
          feature_cols <- setdiff(model_cols, c(".response", ".cohort"))

          if (".cohort" %in% names(stored_model$model)) {
            train_cohort_levels <- levels(stored_model$model$.cohort)
            newdata$.cohort <- factor(rep(train_cohort_levels[1], nrow(newdata)),
                                      levels = train_cohort_levels)
          }

          as.numeric(stats::predict(stored_model, newdata = newdata, type = "response"))
        }
      }, error = function(e) {
        stop(sprintf("Cohort '%s': Failed to generate predictions: %s",
                     cohort_names[i], conditionMessage(e)), call. = FALSE)
      })
    } else {
      scores <- .default_scoring_fn(x_selected, features, yi)
    }

    n_pos <- sum(yi == positive)
    n_neg <- sum(yi != positive)

    metric_values <- vapply(metrics, function(m) {
      switch(m,
        sensitivity = metric_sensitivity(truth = yi, scores = scores,
                                        cutoff_prob = cutoff_prob, positive = positive),
        specificity = metric_specificity(truth = yi, scores = scores,
                                        cutoff_prob = cutoff_prob, positive = positive),
        auc = metric_auc(truth = yi, scores = scores, positive = positive),
        stop(sprintf("Unknown metric: %s", m), call. = FALSE)
      )
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
