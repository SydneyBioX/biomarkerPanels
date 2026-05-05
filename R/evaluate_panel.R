#' Core Panel Evaluation
#'
#' Functions for evaluating biomarker panel performance including confusion
#' matrices, ROC curves, and objective computation.
#'
#' @name evaluate_panel
NULL

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
#' @param feature_transform Optional transformation applied to base features
#'   before evaluation. Defaults to the transform stored in the fitted panel,
#'   keeping training and validation pipelines aligned. Must be a registered
#'   transform name (see [`feature_transform_registry()`]). Built-in options
#'   include `"pairwise_ratios"`, `"pairwise_log_ratios"`, `"reference_norm"`,
#'   and `"none"`.
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
                             metrics = c("sensitivity", "specificity", "auc")
                           ),
                           cohort = NULL,
                           assay = NULL,
                           scoring_fn = NULL,
                           feature_transform = NULL,
                           cutoff_prob = 0.5,
                           positive = NULL) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  if (!is.numeric(cutoff_prob) || length(cutoff_prob) != 1L ||
      is.na(cutoff_prob) || cutoff_prob <= 0 || cutoff_prob >= 1) {
    stop("`cutoff_prob` must be a single numeric value in (0, 1).", call. = FALSE)
  }

  # Use stored positive class from training if not explicitly specified
  if (is.null(positive)) {
    stored_positive <- panel@control$positive_class
    positive <- if (!is.null(stored_positive)) stored_positive else "Yes"
  }

  # Get feature transform from stored control or use provided value
  stored_transform <- panel@control$feature_transform
  if (is.null(feature_transform)) {
    feature_transform <- if (is.null(stored_transform)) "none" else stored_transform
  } else {
    if (!is.character(feature_transform) || length(feature_transform) != 1L) {
      stop("`feature_transform` must be a single character string.", call. = FALSE)
    }
    if (!exists(feature_transform, envir = .transform_registry, inherits = FALSE)) {
      available <- ls(.transform_registry)
      stop(
        "Unknown feature transform '", feature_transform, "'. ",
        "Available: ", paste(available, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # Get base features from panel (original features before transform)
  base_features <- panel@base_features
  if (is.null(base_features) || length(base_features) == 0L) {
    # Backward compatibility: if no base_features, use features directly
    base_features <- panel@features
  }

  # Prepare validation data - get RAW features (no transform yet)
  if (is.list(x)) {
    if (!is.null(cohort)) {
      warning("`cohort` argument is ignored when `x` is supplied as a list of cohorts.",
              call. = FALSE)
    }
    prepared <- .prepare_cohort_inputs(x, y, assay = assay, transform = "none")
    x_raw <- prepared$x
    truth <- prepared$truth
    cohort_vec <- prepared$cohort
  } else {
    x_raw <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_raw) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    if (is.null(colnames(x_raw))) {
      colnames(x_raw) <- sprintf("feature_%04d", seq_len(ncol(x_raw)))
    }
    if (is.null(cohort)) {
      cohort_vec <- factor(rep("cohort_01", nrow(x_raw)), levels = "cohort_01")
    } else {
      if (length(cohort) != nrow(x_raw)) {
        stop("Length of `cohort` must match the number of samples in `x`.", call. = FALSE)
      }
      cohort_vec <- factor(cohort)
    }
  }

  # Validate base features are present in validation data
  if (!all(base_features %in% colnames(x_raw))) {
    missing <- setdiff(base_features, colnames(x_raw))
    stop("Base feature(s) not found in validation data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Extract base features from validation data
  x_base_selected <- x_raw[, base_features, drop = FALSE]

  # Apply feature transform to get transformed features
  if (feature_transform != "none" && length(base_features) >= 2L) {
    x_selected <- .apply_feature_transform_single(x_base_selected, feature_transform)
    selected <- colnames(x_selected)
  } else {
    x_selected <- x_base_selected
    selected <- base_features
  }

  # Validate transformed features match what the model expects
  if (length(selected) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  # Verify feature names match model expectations
  expected_features <- panel@features
  if (!setequal(selected, expected_features)) {
    warning(
      "Transformed feature names differ from model training. ",
      "Expected: ", paste(expected_features, collapse = ", "), ". ",
      "Got: ", paste(selected, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # Determine how to compute scores
  stored_model <- panel@model
  scores <- NULL

  if (!is.null(scoring_fn)) {
    # User provided custom scoring function - use it
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
  } else if (!is.null(stored_model)) {
    # Use the stored model for prediction (true out-of-sample evaluation)
    scores <- tryCatch({
      # Check if this is a glmnet model (regularized) or glm model (unregularized)
      if (inherits(stored_model, "cv.glmnet")) {
        # Regularized model: use glmnet prediction
        preds <- .predict_glmnet_model(
          model = stored_model,
          x_selected = x_selected,
          cohort_vec = cohort_vec
        )
      } else if (inherits(stored_model, "npc")) {
        # NP classifier: use nproc prediction with label inversion
        preds <- .predict_np_model(stored_model, x_selected)
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
          stop(
            "Feature name mismatch between training and validation. ",
            "Expected: ", paste(head(feature_cols, 10), collapse = ", "),
            if (length(feature_cols) > 10) "..." else "", ". ",
            "Got: ", paste(head(expected_names, 10), collapse = ", "),
            if (length(expected_names) > 10) "..." else "",
            call. = FALSE
          )
        }
        # Names match -- ensure same order as model expects
        names(newdata) <- expected_names
        newdata <- newdata[, feature_cols, drop = FALSE]

        # If model was trained with cohort, always use reference level so
        # predictions are cohort-agnostic. Cohort-aware metrics split downstream.
        if (".cohort" %in% names(stored_model$model)) {
          train_cohort_levels <- levels(stored_model$model$.cohort)
          newdata$.cohort <- factor(rep(train_cohort_levels[1], nrow(newdata)),
                                    levels = train_cohort_levels)
        }
        preds <- stats::predict(stored_model, newdata = newdata, type = "response")
      }

      if (length(preds) != nrow(x_selected) || anyNA(preds)) {
        stop("Invalid predictions from stored model.")
      }
      as.numeric(preds)
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
  } else {
    # No model and no scoring function - error
    stop(
      "Panel does not have a fitted model. ",
      "Use fit_panel() to fit a model before calling evaluate_panel(), ",
      "or provide a custom scoring_fn argument.",
      call. = FALSE
    )
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

  sensitivity_point <- metric_sensitivity(
    truth,
    scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )
  specificity_point <- metric_specificity(
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

  roc_auc <- metric_auc(truth, scores, positive = positive)

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

#' Compute Confusion Matrix
#'
#' @param truth Binary response factor.
#' @param scores Numeric prediction scores.
#' @param cutoff_prob Classification cutoff.
#' @param positive Label for positive class.
#' @return A 2x2 confusion matrix with cutoff, positive, and negative attributes.
#' @keywords internal
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

  # Add cohort dummies if the model was trained with them.
  # Always zero out — predictions should be cohort-agnostic for

  # transferability. Cohort-aware metrics split by cohort downstream.
  if (!is.null(meta$cohort_info)) {
    n_dummies <- meta$cohort_info$n_dummies
    dummy_cols <- matrix(0, nrow = nrow(x_mat), ncol = n_dummies)
    colnames(dummy_cols) <- paste0(".cohort_", seq_len(n_dummies))
    x_mat <- cbind(x_mat, dummy_cols)
  }

  # Predict using lambda.min
  preds <- stats::predict(model, newx = x_mat, s = "lambda.min", type = "response")[, 1]

  preds
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
