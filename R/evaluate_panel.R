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

  .validate_probability(cutoff_prob, "cutoff_prob", bounds = "open")

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
    x_raw <- .ensure_feature_colnames(x_raw)
    if (is.null(cohort)) {
      cohort_vec <- factor(rep("cohort_01", nrow(x_raw)), levels = "cohort_01")
    } else {
      if (length(cohort) != nrow(x_raw)) {
        stop("Length of `cohort` must match the number of samples in `x`.", call. = FALSE)
      }
      cohort_vec <- factor(cohort)
    }
  }

  # Validate base features and apply the panel's transform on the fly.
  prepped <- .prepare_scoring_matrix(x_raw, base_features, feature_transform,
                                     context = "validation data")
  x_selected <- prepped$x_selected
  selected <- prepped$selected

  # Validate transformed features match what the model expects
  if (length(selected) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  # Verify feature names match model expectations. CPOP-style panels store a
  # subset of all possible pair features; in that case subset the transformed
  # matrix down to the expected columns rather than warning.
  expected_features <- panel@features
  if (!setequal(selected, expected_features)) {
    if (all(expected_features %in% selected)) {
      x_selected <- x_selected[, expected_features, drop = FALSE]
      selected <- expected_features
    } else {
      stop(
        "Transformed feature names differ from model training. ",
        "Expected: ", paste(expected_features, collapse = ", "), ". ",
        "Got: ", paste(selected, collapse = ", "), ".",
        call. = FALSE
      )
    }
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
      preds <- .predict_panel_model(
        stored_model, x_selected,
        cohort = cohort_vec, expected_features = selected
      )

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
    roc_plot <- ggplot2::ggplot(roc_curve, ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
      ggplot2::geom_path(color = "#377EB8", linewidth = 0.8) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(
        data = highlight_point,
        ggplot2::aes(x = .data$fpr, y = .data$tpr),
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
