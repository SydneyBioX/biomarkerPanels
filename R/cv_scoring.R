#' Cross-Validation Scoring Functions
#'
#' Internal functions for computing cross-validated scores during NSGA
#' optimization. These functions prevent overfitting by using out-of-fold
#' predictions for fitness evaluation.
#'
#' @name cv_scoring
#' @keywords internal
NULL

#' Create Stratified Folds for Cross-Validation
#'
#' Creates class-aware fold assignments that maintain class balance within each fold.
#'
#' @param y Factor response variable.
#' @param k Number of folds.
#' @return Integer vector of fold assignments (1 to k).
#' @keywords internal
.create_stratified_folds <- function(y, k) {

  n <- length(y)
  fold_ids <- integer(n)

  # Shuffle within each class
  for (level in levels(y)) {
    idx <- which(y == level)
    shuffled <- sample(idx)
    fold_ids[shuffled] <- rep(seq_len(k), length.out = length(shuffled))
  }

  fold_ids
}

#' Compute CV-based Scores for Fitness Evaluation
#'
#' For each fold, fits a model on k-1 folds and predicts on the held-out fold.
#' Returns out-of-fold predictions for all samples.
#'
#' @param x_selected Matrix of selected features (samples x features).
#' @param truth Binary response factor.
#' @param fold_ids Integer vector assigning each sample to a fold.
#' @param cohort Optional cohort indicator.
#' @param regularized Logical; if TRUE, use glmnet regularized regression.
#'   Default is FALSE.
#' @param alpha Elastic net mixing parameter when regularized = TRUE
#'   (0 = ridge, 1 = lasso). Default is 0.
#' @param glm_design_terms Optional precomputed GLM design terms for each fold
#'   when `regularized = FALSE`.
#' @return Numeric vector of out-of-fold predicted probabilities.
#' @keywords internal
.compute_cv_scores <- function(x_selected, truth, fold_ids, cohort = NULL,
                               regularized = FALSE, alpha = 0,
                               glm_design_terms = NULL) {
  n <- length(truth)

  # Validate: must have features
  if (ncol(x_selected) == 0L) {
    stop("No features selected for cross-validation scoring.", call. = FALSE)
  }

  # Accumulator for out-of-fold predictions across all samples

  predictions <- rep(NA_real_, n)
  unique_folds <- sort(unique(fold_ids))

  for (fold in unique_folds) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)

    x_train <- x_selected[train_idx, , drop = FALSE]
    y_train <- truth[train_idx]
    x_test <- x_selected[test_idx, , drop = FALSE]

    # Validate: must have both classes in training fold
    if (length(unique(y_train)) < 2L) {
      stop("Training fold ", fold, " contains only one class. ",
           "Consider using more samples or fewer CV folds.", call. = FALSE)
    }

    # Predicted probabilities for this fold's test set
    y_pred <- if (regularized) {
      .fit_cv_fold_regularized(x_train, y_train, x_test, cohort,
                               train_idx, test_idx, fold, alpha)
    } else {
      fold_design_terms <- if (!is.null(glm_design_terms)) {
        glm_design_terms[[as.character(fold)]]
      } else {
        NULL
      }
      .fit_cv_fold_glm(x_train, y_train, x_test, cohort,
                       train_idx, test_idx, fold,
                       design_terms = fold_design_terms)
    }

    predictions[test_idx] <- y_pred
  }

  # Final validation: no NAs should remain
  if (any(is.na(predictions))) {
    stop("Cross-validation scoring produced NA values. This indicates an internal error.",
         call. = FALSE)
  }

  predictions
}

#' Fit GLM for a Single CV Fold
#'
#' @param x_train Training feature matrix.
#' @param y_train Training response factor.
#' @param x_test Test feature matrix.
#' @param cohort Full cohort vector.
#' @param train_idx Training sample indices.
#' @param test_idx Test sample indices.
#' @param fold Fold number (for error messages).
#' @param design_terms Optional precomputed GLM design terms for this fold.
#' @return Numeric vector of predictions for test fold.
#' @keywords internal
.fit_cv_fold_glm <- function(x_train, y_train, x_test, cohort,
                             train_idx, test_idx, fold,
                             design_terms = NULL) {
  tryCatch({
    y_pred <- .fit_predict_binomial_glm(
      x_train = x_train,
      truth = y_train,
      x_new = x_test,
      cohort_train = if (!is.null(cohort)) cohort[train_idx] else NULL,
      cohort_new = if (!is.null(cohort)) cohort[test_idx] else NULL,
      predict_cohort = "observed",
      design_terms = design_terms
    )

    if (anyNA(y_pred)) {
      stop("Model produced NA predictions for fold ", fold, ".", call. = FALSE)
    }

    y_pred
  }, error = function(e) {
    stop("Failed to fit model for fold ", fold, ": ", e$message, call. = FALSE)
  })
}

#' Fit Regularized Model for a Single CV Fold
#'
#' @param x_train Training feature matrix.
#' @param y_train Training response factor.
#' @param x_test Test feature matrix.
#' @param cohort Full cohort vector.
#' @param train_idx Training sample indices.
#' @param test_idx Test sample indices.
#' @param fold Fold number (for error messages).
#' @param alpha Elastic net mixing parameter.
#' @return Numeric vector of predictions for test fold.
#' @keywords internal
.fit_cv_fold_regularized <- function(x_train, y_train, x_test, cohort,
                                     train_idx, test_idx, fold, alpha) {
  tryCatch({
    x_train_mat <- as.matrix(x_train)
    y_train_vec <- as.integer(y_train) - 1L
    x_test_mat <- as.matrix(x_test)

    # Add cohort dummies if provided and has multiple levels
    if (!is.null(cohort) && length(unique(cohort[train_idx])) > 1L) {
      cohort_train <- factor(cohort[train_idx])
      cohort_test <- factor(cohort[test_idx], levels = levels(cohort_train))

      cohort_dummies_train <- stats::model.matrix(~ cohort_train - 1)[, -1, drop = FALSE]
      cohort_dummies_test <- stats::model.matrix(~ cohort_test - 1)[, -1, drop = FALSE]

      if (ncol(cohort_dummies_train) > 0L) {
        colnames(cohort_dummies_train) <- paste0(".cohort_", seq_len(ncol(cohort_dummies_train)))
        colnames(cohort_dummies_test) <- colnames(cohort_dummies_train)
        x_train_mat <- cbind(x_train_mat, cohort_dummies_train)
        x_test_mat <- cbind(x_test_mat, cohort_dummies_test)
      }
    }

    # Determine inner CV folds
    n_train <- nrow(x_train_mat)
    inner_nfolds <- min(5L, max(3L, floor(n_train / 5L)))

    fit <- glmnet::cv.glmnet(
      x = x_train_mat,
      y = y_train_vec,
      family = "binomial",
      alpha = alpha,
      nfolds = inner_nfolds,
      type.measure = "deviance"
    )

    y_pred <- stats::predict(fit, newx = x_test_mat, s = "lambda.min",
                             type = "response")[, 1]

    if (anyNA(y_pred)) {
      stop("Model produced NA predictions for fold ", fold, ".", call. = FALSE)
    }

    y_pred
  }, error = function(e) {
    stop("Failed to fit model for fold ", fold, ": ", e$message, call. = FALSE)
  })
}

#' Default Scoring Function
#'
#' Scores samples using logistic regression (or regularized regression) on
#' selected features. This function is used for in-sample scoring when
#' fitness_cv is FALSE.
#'
#' @param x_selected Matrix of selected features (samples x features).
#' @param selected_features Character vector of selected feature names.
#' @param truth Binary response factor.
#' @param cohort Optional cohort indicator.
#' @param regularized Logical; if TRUE, use glmnet regularized regression.
#'   Default is FALSE.
#' @param alpha Elastic net mixing parameter when regularized = TRUE
#'   (0 = ridge, 1 = lasso, 0.5 = elastic net). Default is 0.5.
#' @param ... Additional arguments (ignored).
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.default_scoring_fn <- function(x_selected, selected_features, truth,
                                cohort = NULL, regularized = FALSE,
                                alpha = 0.5, ...) {
  n <- length(truth)

  # Validate inputs
  if (is.null(x_selected) || ncol(x_selected) == 0L) {
    stop("No features selected for scoring.", call. = FALSE)
  }

  if (length(unique(truth)) < 2L) {
    stop("Response contains only one class. Cannot fit classification model.",
         call. = FALSE)
  }

  # Check for zero variance features
  feature_vars <- apply(x_selected, 2, stats::var, na.rm = TRUE)
  if (all(feature_vars == 0 | is.na(feature_vars))) {
    stop("All selected features have zero variance.", call. = FALSE)
  }

  predictions <- if (regularized) {
    .score_regularized(x_selected, truth, cohort, alpha, n)
  } else {
    .score_glm(x_selected, truth, cohort, n)
  }

  predictions
}

#' Score Samples Using GLM
#'
#' @param x_selected Matrix of selected features.
#' @param truth Binary response factor.
#' @param cohort Optional cohort indicator.
#' @param n Number of samples.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.score_glm <- function(x_selected, truth, cohort, n) {
  tryCatch({
    predictions <- .fit_predict_binomial_glm(
      x_train = x_selected,
      truth = truth,
      x_new = x_selected,
      cohort_train = cohort,
      cohort_new = cohort,
      predict_cohort = "observed"
    )

    if (length(predictions) != n) {
      stop("Model predictions have incorrect length (expected ", n,
           ", got ", length(predictions), ").", call. = FALSE)
    }
    if (anyNA(predictions)) {
      stop("Model produced NA predictions.", call. = FALSE)
    }

    predictions
  }, error = function(e) {
    stop("Failed to fit scoring model: ", e$message, call. = FALSE)
  })
}

#' Score Samples Using Regularized Regression
#'
#' @param x_selected Matrix of selected features.
#' @param truth Binary response factor.
#' @param cohort Optional cohort indicator.
#' @param alpha Elastic net mixing parameter.
#' @param n Number of samples.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.score_regularized <- function(x_selected, truth, cohort, alpha, n) {
  tryCatch({
    x_mat <- as.matrix(x_selected)
    y_vec <- as.integer(truth) - 1L

    # Add cohort dummies if provided
    if (!is.null(cohort) && length(unique(cohort)) > 1L) {
      cohort_dummies <- stats::model.matrix(~ factor(cohort) - 1)[, -1, drop = FALSE]
      if (ncol(cohort_dummies) > 0L) {
        colnames(cohort_dummies) <- paste0(".cohort_", seq_len(ncol(cohort_dummies)))
        x_mat <- cbind(x_mat, cohort_dummies)
      }
    }

    nfolds <- min(5L, max(3L, floor(n / 5L)))

    fit <- glmnet::cv.glmnet(
      x = x_mat,
      y = y_vec,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      type.measure = "deviance"
    )

    predictions <- stats::predict(fit, newx = x_mat, s = "lambda.min",
                                   type = "response")[, 1]

    if (length(predictions) != n) {
      stop("Model predictions have incorrect length (expected ", n,
           ", got ", length(predictions), ").", call. = FALSE)
    }
    if (anyNA(predictions)) {
      stop("Model produced NA predictions.", call. = FALSE)
    }

    predictions
  }, error = function(e) {
    stop("Failed to fit scoring model: ", e$message, call. = FALSE)
  })
}
