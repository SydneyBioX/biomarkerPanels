#' Stored-Model Prediction Dispatch
#'
#' Single entry point for turning a stored panel model into a numeric vector of
#' predicted probabilities. Consolidates the prediction policy that was
#' previously duplicated across [evaluate_panel()], [evaluate_panel_by_cohort()],
#' [calibrate_panel()], and the NSGA fitness paths, so those consumers cannot
#' drift apart in how they score glmnet, glm, and Neyman-Pearson models.
#'
#' @name model_prediction
#' @keywords internal
NULL

#' Predict from a Stored Panel Model
#'
#' Dispatches on model type and returns predicted probabilities. Cohort-aware
#' models are scored cohort-agnostically (glmnet dummies zeroed, glm `.cohort`
#' set to its training reference level) because cohort-aware metrics split by
#' cohort downstream.
#'
#' @param model Fitted model: a `cv.glmnet`, an `npc` classifier, or a `glm`.
#' @param x_selected Numeric matrix of (already transformed) selected features.
#' @param cohort Optional cohort factor. Retained for signature compatibility
#'   with the fitness paths; predictions are cohort-agnostic so it is unused.
#' @param expected_features Optional character vector of the transformed feature
#'   names the model was trained on. When supplied (the evaluation paths), glm
#'   prediction validates that the training and prediction feature names match
#'   and reorders columns to the model's expectation. When `NULL` (the fitness
#'   paths, where names match by construction) the check is skipped.
#' @return Numeric vector of predicted probabilities, length `nrow(x_selected)`.
#' @keywords internal
.predict_panel_model <- function(model, x_selected, cohort = NULL,
                                 expected_features = NULL) {
  if (inherits(model, "cv.glmnet")) {
    x_mat <- as.matrix(x_selected)

    # Add cohort dummies if the model was trained with them, zeroed so that
    # prediction is cohort-agnostic.
    meta <- model$biomarkerPanels_meta
    if (!is.null(meta$cohort_info)) {
      n_dummies <- meta$cohort_info$n_dummies
      dummy_cols <- matrix(0, nrow = nrow(x_mat), ncol = n_dummies)
      colnames(dummy_cols) <- paste0(".cohort_", seq_len(n_dummies))
      x_mat <- cbind(x_mat, dummy_cols)
    }

    preds <- stats::predict(model, newx = x_mat, s = "lambda.min",
                            type = "response")[, 1]
  } else if (inherits(model, "npc")) {
    preds <- .predict_np_model(model, x_selected)
  } else if (inherits(model, "glm")) {
    # Note: as.data.frame.matrix() ignores check.names, so we sanitize names
    # explicitly to match the names used during training in .fit_final_model(),
    # which wraps in data.frame(check.names = TRUE).
    newdata <- as.data.frame(x_selected)
    names(newdata) <- make.names(names(newdata), unique = TRUE)

    model_cols <- names(model$model)
    feature_cols <- setdiff(model_cols, c(".response", ".cohort"))

    # Evaluation paths pass the panel's expected feature names: validate the
    # match and reorder to the model's expectation. CPOP-style panels train on a
    # subset of the transformed columns, so reordering picks the right columns.
    if (!is.null(expected_features)) {
      expected_names <- make.names(expected_features, unique = TRUE)
      if (!setequal(expected_names, feature_cols)) {
        stop(
          "Feature name mismatch between training and prediction data. ",
          "Expected: ", paste(utils::head(feature_cols, 10L), collapse = ", "),
          if (length(feature_cols) > 10L) "..." else "", ". ",
          "Got: ", paste(utils::head(expected_names, 10L), collapse = ", "),
          if (length(expected_names) > 10L) "..." else "",
          call. = FALSE
        )
      }
      names(newdata) <- expected_names
      newdata <- newdata[, feature_cols, drop = FALSE]
    }

    # If the model was trained with a cohort covariate, always use the reference
    # level so predictions are cohort-agnostic (matching the glmnet path).
    if (".cohort" %in% model_cols) {
      train_cohort_levels <- levels(model$model$.cohort)
      newdata$.cohort <- factor(rep(train_cohort_levels[1L], nrow(newdata)),
                                levels = train_cohort_levels)
    }

    preds <- stats::predict(model, newdata = newdata, type = "response")
  } else {
    stop("Unknown model type: ", class(model)[1L], call. = FALSE)
  }

  as.numeric(preds)
}

#' Predict from NP Classifier
#'
#' Extracts continuous probability scores from an npc model object.
#' Handles label inversion when labels were flipped for rule-out classification.
#'
#' @param model An npc object from nproc::npc().
#' @param x_selected Matrix of selected features for prediction.
#' @return Numeric vector of predicted probabilities on the original scale.
#' @keywords internal
.predict_np_model <- function(model, x_selected) {
  .check_nproc_available()

  x_mat <- as.matrix(x_selected)

  labels_flipped <- isTRUE(model$biomarkerPanels_meta$labels_flipped)

  # Get predictions from nproc - use predict method for npc objects
  np_pred <- stats::predict(model, newx = x_mat)

  # npc predict returns a list with pred.label and pred.score
  if (is.list(np_pred) && "pred.score" %in% names(np_pred)) {
    scores <- as.numeric(np_pred$pred.score)
  } else if (is.numeric(np_pred)) {
    scores <- np_pred
  } else {
    stop(
      "Unexpected prediction format from npc model. ",
      "Expected list with pred.score or numeric vector.",
      call. = FALSE
    )
  }

  # Invert if labels were flipped
  .invert_np_predictions(scores, labels_flipped)
}
