#' Final Model Fitting Functions
#'
#' Internal functions for fitting the final model after NSGA-II optimization
#' has selected the best features. These models are stored in the
#' BiomarkerPanelResult object for use in evaluation.
#'
#' @name model_fitting
#' @keywords internal
NULL

#' Fit Final Model for Storage
#'
#' Fits an unregularized logistic regression model on the selected features.
#'
#' @param x_selected Matrix of selected features.
#' @param truth Binary response factor.
#' @param cohort Optional cohort factor.
#' @return GLM model object, or NULL on failure.
#' @keywords internal
.fit_final_model <- function(x_selected, truth, cohort = NULL) {
  if (is.null(x_selected) || ncol(x_selected) == 0L) {
    return(NULL)
  }

  tryCatch({
    df <- data.frame(
      .response = as.integer(truth) - 1L,
      as.data.frame(x_selected, check.names = TRUE)
    )
    # Add cohort as a covariate if provided and has multiple levels
    if (!is.null(cohort) && length(unique(cohort)) > 1L) {
      df$.cohort <- factor(cohort)
    }
    fit <- stats::glm(.response ~ ., data = df, family = stats::binomial())
    fit
  }, error = function(e) {
    stop(
      "Failed to fit final logistic regression model: ",
      conditionMessage(e),
      "\nThis typically indicates perfect separation, multicollinearity, or insufficient data. ",
      "Consider reducing the number of features or using regularized = TRUE.",
      call. = FALSE
    )
  })
}

#' Fit Final Model with Cross-Validation to Reduce Training Contamination
#'
#' Uses k-fold cross-validation to fit the final model. Each fold trains on
#' the training subset and coefficients are averaged across folds. This
#' reduces training data contamination compared to fitting on all training data.
#'
#' @param x_selected Matrix of selected features.
#' @param truth Binary response factor.
#' @param cohort Optional cohort factor.
#' @param cv_folds Number of CV folds.
#' @return GLM model with averaged coefficients, or NULL on failure.
#' @keywords internal
.fit_final_model_cv <- function(x_selected, truth, cohort = NULL, cv_folds = 5L) {
  if (is.null(x_selected) || ncol(x_selected) == 0L) {
    return(NULL)
  }

  n <- nrow(x_selected)
  if (n < cv_folds * 2L) {
    # Not enough data for CV, fall back to standard fitting
    return(.fit_final_model(x_selected, truth, cohort))
  }

  tryCatch({
    # Prepare full data frame
    df <- data.frame(
      .response = as.integer(truth) - 1L,
      as.data.frame(x_selected, check.names = TRUE)
    )
    if (!is.null(cohort) && length(unique(cohort)) > 1L) {
      df$.cohort <- factor(cohort)
    }

    # Create stratified folds to maintain class balance
    fold_ids <- .create_stratified_folds(truth, cv_folds)

    # Collect coefficients from each fold
    all_coefs <- list()
    for (fold in seq_len(cv_folds)) {
      train_idx <- fold_ids != fold
      train_df <- df[train_idx, , drop = FALSE]

      fold_fit <- tryCatch(
        stats::glm(.response ~ ., data = train_df, family = stats::binomial()),
        error = function(e) NULL
      )

      if (!is.null(fold_fit)) {
        all_coefs[[length(all_coefs) + 1L]] <- stats::coef(fold_fit)
      }
    }

    if (length(all_coefs) == 0L) {
      stop(
        "All cross-validation folds failed to fit. ",
        "This indicates a systematic problem with the data or model specification. ",
        "Check for zero-variance features, perfect separation, or insufficient samples.",
        call. = FALSE
      )
    }
    if (length(all_coefs) < cv_folds) {
      warning(
        sprintf("%d of %d CV folds failed. Results may be less reliable.",
                cv_folds - length(all_coefs), cv_folds),
        call. = FALSE
      )
    }

    # Average coefficients across folds
    coef_names <- names(all_coefs[[1]])
    avg_coefs <- vapply(coef_names, function(nm) {
      vals <- vapply(all_coefs, function(cf) {
        if (nm %in% names(cf)) cf[[nm]] else NA_real_
      }, numeric(1))
      mean(vals, na.rm = TRUE)
    }, numeric(1))

    # Fit final model on all data to get the model object, then replace coefficients
    full_fit <- stats::glm(.response ~ ., data = df, family = stats::binomial())
    full_fit$coefficients <- avg_coefs
    full_fit$cv_folds_used <- length(all_coefs)
    full_fit$cv_averaged <- TRUE

    full_fit
  }, error = function(e) {
    stop(
      "CV model fitting failed: ", conditionMessage(e),
      "\nCheck for data issues such as zero-variance features or insufficient samples.",
      call. = FALSE
    )
  })
}

#' Fit Final Regularized Model for Storage
#'
#' Fits a regularized logistic regression model (elastic net via glmnet) on the
#' selected features. The fitted cv.glmnet object is stored in the
#' BiomarkerPanelResult for use in evaluation.
#'
#' @param x_selected Matrix of selected features.
#' @param truth Binary response factor.
#' @param cohort Optional cohort factor.
#' @param alpha Elastic net mixing parameter.
#' @return cv.glmnet model object, or NULL on failure.
#' @keywords internal
.fit_final_model_regularized <- function(x_selected, truth, cohort = NULL,
                                          alpha = 0.5) {
  if (is.null(x_selected) || ncol(x_selected) == 0L) {
    return(NULL)
  }

  n <- nrow(x_selected)
  x_mat <- as.matrix(x_selected)
  y_vec <- as.integer(truth) - 1L

  # Add cohort dummies if provided
  cohort_info <- NULL
  if (!is.null(cohort) && length(unique(cohort)) > 1L) {
    cohort_factor <- factor(cohort)
    cohort_dummies <- stats::model.matrix(~ cohort_factor - 1)[, -1, drop = FALSE]
    if (ncol(cohort_dummies) > 0L) {
      colnames(cohort_dummies) <- paste0(".cohort_", seq_len(ncol(cohort_dummies)))
      x_mat <- cbind(x_mat, cohort_dummies)
      cohort_info <- list(
        levels = levels(cohort_factor),
        n_dummies = ncol(cohort_dummies)
      )
    }
  }

  tryCatch({
    nfolds <- min(10L, max(3L, floor(n / 5L)))

    fit <- glmnet::cv.glmnet(
      x = x_mat,
      y = y_vec,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      type.measure = "deviance"
    )

    # Store additional metadata for prediction
    fit$biomarkerPanels_meta <- list(
      feature_names = colnames(x_selected),
      cohort_info = cohort_info,
      alpha = alpha
    )

    fit
  }, error = function(e) {
    stop(
      "Failed to fit regularized model: ", conditionMessage(e),
      "\nEnsure glmnet is installed and data contains sufficient variation.",
      call. = FALSE
    )
  })
}
