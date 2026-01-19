#' Evaluate a Biomarker Panel
#'
#' Compute objective values and additional summary statistics for a fitted panel
#' on validation data. This function extracts the selected features from the
#' validation set and computes predictions using the default scoring function:
#' a logistic regression fit on the selected features (with a centred/standardised
#' fallback if the model fails to converge).
#' When evaluating multiple cohorts, feature alignment is performed via the
#' simple intersection of shared column names; future releases will offer more
#' flexible strategies.
#'
#' @param panel A `BiomarkerPanelResult`.
#' @param x Validation feature matrix, `data.frame`, `SummarizedExperiment`, or a
#'   list of such objects for multi-cohort evaluation.
#' @param y Validation binary response factor (or list of factors when `x` is a
#'   list).
#' @param objectives Optional override for objectives. Defaults to sensitivity,
#'   specificity, and auc.
#' @param cohort Optional factor identifying cohort membership for each sample
#'   when `x` represents a single cohort. Ignored when `x` is provided as a list.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @param scoring_fn Optional custom scoring function. If `NULL`, uses the
#'   default row-mean aggregation. Must have signature
#'   `function(x_selected, selected_features, truth, ...)`.
#' @param cohort_aggregator Optional transformation applied to `x` before
#'   evaluation. Defaults to the aggregator stored in the fitted panel,
#'   keeping training and validation pipelines aligned. Must be a registered
#'   aggregator name (see [`aggregator_registry()`]). Built-in options include
#'   `"pairwise_ratios"`, `"pairwise_log_ratios"`, `"reference_norm"`, and `"none"`.
#' @param cutoff_prob Classification probability threshold used for confusion
#'   matrix summaries and highlight point on the ROC curve. Defaults to `0.5`.
#' @param positive Label treated as the positive class when computing confusion
#'   matrices. Defaults to the positive class stored during training (typically
#'   `"Yes"`). Specify explicitly to override.
#' @return A list with `metrics` (named numeric vector), `objectives`
#'   (data.frame with columns objective, value, direction), `confusion`
#'   (2x2 matrix of counts), and `roc` (list containing a `curve` data.frame,
#'   `highlight` point, optional `plot`, and `evalm` helper data).
#' @export
evaluate_panel <- function(panel, x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "auc")
                           ),
                           cohort = NULL,
                           assay = NULL,
                           scoring_fn = NULL,
                           cohort_aggregator = NULL,
                           cutoff_prob = 0.5,
                           positive = NULL) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  # Use stored positive class from training if not explicitly specified
  if (is.null(positive)) {
    stored_positive <- panel@control$positive_class
    positive <- if (!is.null(stored_positive)) stored_positive else "Yes"
  }

  stored_aggregator <- panel@control$cohort_aggregator
  if (is.null(cohort_aggregator)) {
    cohort_aggregator <- if (is.null(stored_aggregator)) "none" else stored_aggregator
  } else {
    if (!is.character(cohort_aggregator) || length(cohort_aggregator) != 1L) {
      stop("`cohort_aggregator` must be a single character string.", call. = FALSE)
    }
    if (!exists(cohort_aggregator, envir = .aggregator_registry, inherits = FALSE)) {
      available <- ls(.aggregator_registry)
      stop(
        "Unknown aggregator '", cohort_aggregator, "'. ",
        "Available: ", paste(available, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (is.list(x)) {
    if (!is.null(cohort)) {
      warning("`cohort` argument is ignored when `x` is supplied as a list of cohorts.",
              call. = FALSE)
    }
    prepared <- .prepare_cohort_inputs(x, y, assay = assay, aggregator = cohort_aggregator)
    x_mat <- prepared$x
    truth <- prepared$truth
    cohort_vec <- prepared$cohort
  } else {
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    if (is.null(colnames(x_mat))) {
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }
    x_mat <- .apply_cohort_aggregator(list(x_mat), cohort_aggregator)[[1]]
    if (is.null(colnames(x_mat))) {
      stop("`x` must have column names in order to align with panel features.",
           call. = FALSE)
    }
    if (is.null(cohort)) {
      cohort_vec <- factor(rep("cohort_01", nrow(x_mat)), levels = "cohort_01")
    } else {
      if (length(cohort) != nrow(x_mat)) {
        stop("Length of `cohort` must match the number of samples in `x`.", call. = FALSE)
      }
      cohort_vec <- factor(cohort)
    }
  }

  selected <- panel@features
  if (length(selected) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  if (!all(selected %in% colnames(x_mat))) {
    missing <- setdiff(selected, colnames(x_mat))
    stop("Selected feature(s) not found in validation data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  x_selected <- x_mat[, selected, drop = FALSE]

  # Check if panel has a stored model - use it for prediction if available
  stored_model <- panel@model
  if (!is.null(stored_model) && is.null(scoring_fn)) {
    # Use the stored model for prediction (true out-of-sample evaluation)
    scores <- tryCatch({
      # Check if this is a glmnet model (regularized) or glm model (unregularized)
      if (inherits(stored_model, "cv.glmnet")) {
        # Regularized model: use glmnet prediction
        scores <- .predict_glmnet_model(
          model = stored_model,
          x_selected = x_selected,
          cohort_vec = cohort_vec
        )
      } else {
        # Standard GLM model: use standard prediction
        newdata <- as.data.frame(x_selected, check.names = TRUE)

        # Get expected column names from the model (excluding .response and .cohort)
        model_cols <- names(stored_model$model)
        feature_cols <- setdiff(model_cols, c(".response", ".cohort"))

        # Robust column name matching: the model was fit with check.names=TRUE,
        # so we need to match original feature names to their transformed versions.
        # Create the expected mapping by applying check.names to panel features.
        expected_names <- make.names(selected, unique = TRUE)

        # Verify the mapping matches what the model expects
        if (!setequal(expected_names, feature_cols)) {
          # Names don't match - this could happen if features were renamed
          # differently during training. Try position-based mapping as fallback.
          if (length(feature_cols) == ncol(x_selected)) {
            warning(
              "Feature name mismatch between training and validation. ",
              "Using positional mapping (verify features are aligned).",
              call. = FALSE
            )
            names(newdata) <- feature_cols
          } else {
            stop("Cannot match validation features to model features.")
          }
        } else {
          # Names match - ensure same order as model expects
          names(newdata) <- expected_names
          # Reorder to match model's feature order
          newdata <- newdata[, feature_cols, drop = FALSE]
        }

        # Add cohort if the model was trained with cohort
        if (".cohort" %in% names(stored_model$model)) {
          if (length(unique(cohort_vec)) > 1L) {
            newdata$.cohort <- factor(cohort_vec)
          } else {
            # If validation has only one cohort, use the first training cohort level
            train_cohort_levels <- levels(stored_model$model$.cohort)
            newdata$.cohort <- factor(rep(train_cohort_levels[1], nrow(newdata)),
                                      levels = train_cohort_levels)
          }
        }
        scores <- stats::predict(stored_model, newdata = newdata, type = "response")
      }

      if (length(scores) != nrow(x_selected) || anyNA(scores)) {
        stop("Invalid predictions from stored model.")
      }
      as.numeric(scores)
    }, error = function(e) {
      stop(
        "Failed to generate predictions from stored model: ",
        conditionMessage(e),
        "\nPrediction on validation data must use the model trained on training data. ",
        "Refitting on validation data would invalidate the evaluation. ",
        "Ensure validation features match training features.",
        call. = FALSE
      )
    })
  } else if (is.null(scoring_fn)) {
    scoring_fn <- .default_scoring_fn
  }

  # If scores not yet computed (no stored model or it failed), use scoring_fn
 if (!exists("scores") || is.null(scores)) {
    if (!is.function(scoring_fn)) {
      stop("`scoring_fn` must be a function.", call. = FALSE)
    }

    score_args <- list(
      x_selected = x_selected,
      selected_features = selected,
      truth = truth,
      cohort = cohort_vec
    )
    scores <- do.call(scoring_fn, score_args)
  }

  if (!is.numeric(scores) || length(scores) != length(truth)) {
    stop("`scoring_fn` must return a numeric vector matching the number of samples.",
         call. = FALSE)
  }

  objective_values <- vapply(objectives, function(obj) {
    obj$fun(
      truth,
      scores,
      selected = selected,
      cohort = cohort_vec,
      x = x_selected
    )
  }, numeric(1))

  sensitivity_point <- loss_sensitivity(
    truth,
    scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )
  specificity_point <- loss_specificity(
    truth,
    scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )

  confusion <- .compute_confusion_matrix(
    truth = truth,
    scores = scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )

  roc_curve <- .compute_roc_curve(
    truth = truth,
    scores = scores,
    positive = positive
  )

  highlight_point <- data.frame(
    threshold = cutoff_prob,
    tpr = sensitivity_point,
    fpr = 1 - specificity_point,
    sensitivity = sensitivity_point,
    specificity = specificity_point,
    stringsAsFactors = FALSE
  )

  roc_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    roc_plot <- ggplot2::ggplot(roc_curve, ggplot2::aes(x = fpr, y = tpr)) +
      ggplot2::geom_path(color = "#377EB8", linewidth = 0.8) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(
        data = highlight_point,
        ggplot2::aes(x = fpr, y = tpr),
        color = "#E41A1C",
        size = 2.5
      ) +
      ggplot2::labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        title = "ROC Curve"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal()
  }

  roc_auc <- loss_auc(truth, scores, positive = positive)

  list(
    metrics = objective_values,
    objectives = data.frame(
      objective = names(objectives),
      value = objective_values,
      direction = vapply(objectives, `[[`, character(1), "direction"),
      stringsAsFactors = FALSE
    ),
    confusion = confusion,
    roc = list(
      curve = roc_curve,
      highlight = highlight_point,
      auc = roc_auc,
      plot = roc_plot,
      evalm = data.frame(
        obs = truth,
        prob = scores,
        stringsAsFactors = FALSE
      )
    ),
    cutoff = cutoff_prob,
    scores = scores
  )
}

.compute_confusion_matrix <- function(truth, scores, cutoff_prob, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }
  negative <- setdiff(levels_truth, positive)
  if (!length(negative)) {
    stop("Unable to determine negative class for confusion matrix.", call. = FALSE)
  }
  negative <- negative[1]
  predicted_positive <- scores >= cutoff_prob
  predicted <- ifelse(predicted_positive, positive, negative)
  predicted <- factor(predicted, levels = levels_truth)
  table_res <- table(
    truth = truth,
    predicted = predicted
  )
  structure(
    as.matrix(table_res),
    cutoff = cutoff_prob,
    positive = positive,
    negative = negative
  )
}

#' Compute ROC Curve
#'
#' Computes the ROC curve for a set of predictions. Uses C++ for performance.
#'
#' @param truth A binary response vector or factor.
#' @param scores Numeric vector of predicted scores or probabilities.
#' @param positive The label for the positive class.
#' @return A data.frame with columns `threshold`, `tpr`, `fpr`, `sensitivity`,
#'   and `specificity`.
#' @keywords internal
.compute_roc_curve <- function(truth, scores, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }

  # Convert truth to logical vector for C++
  is_positive <- truth == positive

  # Call C++ implementation
  df <- .compute_roc_curve_cpp(as.numeric(scores), is_positive)

  # Sort by fpr, tpr for consistency with pure R version
  df[order(df$fpr, df$tpr), , drop = FALSE]
}

#' Predict Using Stored glmnet Model
#'
#' Helper function to make predictions from a cv.glmnet model stored in a
#' BiomarkerPanelResult. Handles cohort dummy variables using the metadata
#' stored during training.
#'
#' @param model A cv.glmnet model object with biomarkerPanels_meta attribute.
#' @param x_selected Matrix of selected features for prediction.
#' @param cohort_vec Factor of cohort membership for samples.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_glmnet_model <- function(model, x_selected, cohort_vec) {
  x_mat <- as.matrix(x_selected)

  # Get metadata stored during training

  meta <- model$biomarkerPanels_meta

  # Add cohort dummies if the model was trained with them
  if (!is.null(meta$cohort_info)) {
    cohort_info <- meta$cohort_info
    train_levels <- cohort_info$levels
    n_dummies <- cohort_info$n_dummies

    # Create cohort factor with training levels
    cohort_factor <- factor(cohort_vec, levels = train_levels)

    # Handle unseen cohort levels - this is a data integrity error
    if (any(is.na(cohort_factor))) {
      unseen <- unique(as.character(cohort_vec)[is.na(cohort_factor)])
      stop(
        "Validation data contains cohort levels not seen during training: ",
        paste(unseen, collapse = ", "), ".\n",
        "Training cohorts: ", paste(train_levels, collapse = ", "), ".\n",
        "Cannot make valid predictions for unseen cohorts. ",
        "Either include these cohorts in training or exclude them from validation.",
        call. = FALSE
      )
    }

    cohort_dummies <- stats::model.matrix(~ cohort_factor - 1)[, -1, drop = FALSE]

    # Ensure correct number of dummy columns
    if (ncol(cohort_dummies) != n_dummies) {
      # Pad with zeros if needed (for cohorts not in validation data)
      if (ncol(cohort_dummies) < n_dummies) {
        padding <- matrix(0, nrow = nrow(cohort_dummies),
                          ncol = n_dummies - ncol(cohort_dummies))
        cohort_dummies <- cbind(cohort_dummies, padding)
      } else {
        cohort_dummies <- cohort_dummies[, seq_len(n_dummies), drop = FALSE]
      }
    }

    colnames(cohort_dummies) <- paste0(".cohort_", seq_len(n_dummies))
    x_mat <- cbind(x_mat, cohort_dummies)
  }

  # Predict using lambda.min
  preds <- stats::predict(model, newx = x_mat, s = "lambda.min", type = "response")[, 1]

  preds
}

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
#' @param cohort_aggregator Optional aggregator name. Defaults to the aggregator
#'   stored in the fitted panel.
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
                                     cohort_aggregator = NULL,
                                     assay = NULL) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  features <- panel@features
  if (length(features) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  stored_model <- panel@model

  if (is.null(positive)) {
    stored_positive <- panel@control$positive_class
    positive <- if (!is.null(stored_positive)) stored_positive else "Yes"
  }

  stored_aggregator <- panel@control$cohort_aggregator
  if (is.null(cohort_aggregator)) {
    cohort_aggregator <- if (is.null(stored_aggregator)) "none" else stored_aggregator
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

    xi <- .apply_cohort_aggregator(list(xi), cohort_aggregator)[[1]]

    if (!all(features %in% colnames(xi))) {
      missing <- setdiff(features, colnames(xi))
      stop(sprintf("Cohort '%s': Feature(s) not found: %s",
                   cohort_names[i], paste(missing, collapse = ", ")),
           call. = FALSE)
    }

    x_selected <- xi[, features, drop = FALSE]

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
        sensitivity = loss_sensitivity(truth = yi, scores = scores,
                                        cutoff_prob = cutoff_prob, positive = positive),
        specificity = loss_specificity(truth = yi, scores = scores,
                                        cutoff_prob = cutoff_prob, positive = positive),
        auc = loss_auc(truth = yi, scores = scores, positive = positive),
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

# ==============================================================================
# PURE R REFERENCE IMPLEMENTATIONS (for testing)
# TODO: Remove these once Rcpp equivalents are fully debugged and validated.
# See tests/testthat/test-rcpp-equivalence.R for equivalence tests.
# ==============================================================================

#' Pure R Reference Implementation of ROC Curve Computation
#'
#' Kept for regression testing against the C++ implementation.
#' TODO: Remove once Rcpp equivalents are fully debugged and validated.
#'
#' @inheritParams .compute_roc_curve
#' @keywords internal
.compute_roc_curve_pure_r <- function(truth, scores, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }
  negative <- setdiff(levels_truth, positive)
  if (!length(negative)) {
    stop("Unable to determine negative class for ROC.", call. = FALSE)
  }
  negative <- negative[1]

  thresholds <- sort(unique(scores), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  pos_total <- sum(truth == positive)
  neg_total <- sum(truth == negative)

  compute_point <- function(thresh) {
    predicted_positive <- scores >= thresh
    tp <- sum(predicted_positive & truth == positive)
    fp <- sum(predicted_positive & truth == negative)
    tpr <- if (pos_total == 0) NA_real_ else tp / pos_total
    fpr <- if (neg_total == 0) NA_real_ else fp / neg_total
    c(tpr = tpr, fpr = fpr)
  }

  mat <- vapply(thresholds, compute_point, numeric(2))
  df <- data.frame(
    threshold = thresholds,
    tpr = mat["tpr", ],
    fpr = mat["fpr", ],
    sensitivity = mat["tpr", ],
    specificity = 1 - mat["fpr", ],
    stringsAsFactors = FALSE
  )
  df[order(df$fpr, df$tpr), , drop = FALSE]
}
