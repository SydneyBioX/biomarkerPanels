#' Final Model Fitting Functions
#'
#' Internal functions for fitting the final model after NSGA optimization
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
    warning(
      "Too few samples (n=", n, ") for ", cv_folds, "-fold CV. ",
      "Falling back to standard (non-CV) model fitting.",
      call. = FALSE
    )
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
        error = function(e) {
          warning(
            sprintf("CV fold %d/%d failed: %s", fold, cv_folds, conditionMessage(e)),
            call. = FALSE
          )
          NULL
        }
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

    non_finite <- names(avg_coefs)[!is.finite(avg_coefs)]
    if (length(non_finite) > 0L) {
      stop(
        "CV averaging produced non-finite coefficients for: ",
        paste(non_finite, collapse = ", "),
        ". This typically means all CV folds failed for these terms.",
        call. = FALSE
      )
    }

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
    fit <- .fit_binomial_glmnet(x_mat, y_vec, alpha)

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

#' Fit a Panel Model from Optimization Results
#'
#' Takes an `OptimizationResult` from [optimize_panel()] and fits a model on
#' a selected Pareto solution. Returns a `BiomarkerPanelResult` with the
#' fitted model ready for evaluation.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @param solution_id Integer specifying which Pareto solution to use. If
#'   `NULL` (default), automatically selects the solution with the best
#'   value on the first objective.
#' @param features Character vector of features to use. If provided, overrides
#'   the features from `solution_id`. Must be a subset of the feature pool.
#' @param x Optional training feature matrix. If `NULL`, uses the stored
#'   aggregated training data from the optimization result.
#' @param y Optional training response vector. If `NULL`, uses the stored
#'   response from the optimization result.
#' @param cohort Optional cohort factor. If `NULL`, uses the stored cohort
#'   from the optimization result.
#' @param regularized Logical; if `TRUE` (default), use regularized regression
#'   (elastic net via glmnet). If `FALSE`, use standard logistic regression.
#' @param regularized_alpha Elastic net mixing parameter when `regularized = TRUE`.
#'   `alpha = 1` is lasso, `alpha = 0` is ridge. Default is 0.5 (elastic net).
#' @param cv Logical; if `TRUE`, use cross-validation when fitting the final
#'   model (only applies when `regularized = FALSE`). Default is `FALSE`.
#' @param cv_folds Number of cross-validation folds when `cv = TRUE`.
#'   Default is 5.
#' @return A `BiomarkerPanelResult` with the fitted model.
#' @export
#' @seealso [optimize_panel()], [evaluate_panel()]
#' @examples
#' \dontrun{
#' # Run optimization
#' opt <- optimize_panel(x, y, objectives = define_objectives())
#'
#' # Inspect solutions
#' summarize_solutions(opt)
#'
#' # Fit model on best solution (auto-selected)
#' panel <- fit_panel(opt)
#'
#' # Or specify a particular solution
#' panel <- fit_panel(opt, solution_id = 3)
#'
#' # Evaluate on held-out data
#' evaluate_panel(panel, x_test, y_test)
#' }
fit_panel <- function(optimization_result,
                      solution_id = NULL,
                      features = NULL,
                      x = NULL,
                      y = NULL,
                      cohort = NULL,
                      regularized = TRUE,
                      regularized_alpha = 0.5,
                      cv = FALSE,
                      cv_folds = 5L) {

  .validate_probability(regularized_alpha, "regularized_alpha", bounds = "closed")
  .validate_positive_integer(cv_folds, "cv_folds", min = 2L)

  sol <- .resolve_panel_solution(optimization_result, solution_id, features,
                                 x, y, cohort)
  selected_base_features <- sol$base_features
  selected_features <- sol$features
  selected_solution_id <- sol$solution_id
  solution_metrics <- sol$metrics
  x_raw <- sol$x_raw
  x_selected <- sol$x_selected
  truth <- sol$truth
  cohort_vec <- sol$cohort

  # Fit the model
  if (regularized) {
    model <- .fit_final_model_regularized(
      x_selected, truth, cohort_vec,
      alpha = regularized_alpha
    )
  } else if (cv) {
    model <- .fit_final_model_cv(
      x_selected, truth, cohort_vec,
      cv_folds = cv_folds
    )
  } else {
    model <- .fit_final_model(x_selected, truth, cohort_vec)
  }

  if (is.null(model)) {
    stop("Model fitting failed. Check the data and parameters.", call. = FALSE)
  }

  # Build metrics from solution or recompute
  if (!is.null(solution_metrics)) {
    metrics <- solution_metrics
  } else {
    # For custom features, set metrics to NA (user should evaluate)
    metrics <- numeric(0)
  }

  # Build objectives data frame for compatibility
  objective_df <- data.frame(
    solution_id = if (is.na(selected_solution_id)) 1L else selected_solution_id,
    stringsAsFactors = FALSE
  )
  if (length(metrics) > 0L) {
    for (nm in names(metrics)) {
      objective_df[[nm]] <- metrics[[nm]]
    }
  }
  objective_df$features <- I(list(selected_features))

  # Compute training data signature
  control <- optimization_result@control
  training_data <- list(
    n = nrow(x_raw),
    p = ncol(x_raw),
    class_balance = table(truth),
    feature_pool_size = length(optimization_result@feature_pool),
    num_cohorts = if (!is.null(cohort_vec)) length(unique(cohort_vec)) else 1L
  )

  # Create BiomarkerPanelResult with both base and transformed features
  panel <- new(
    "BiomarkerPanelResult",
    base_features = selected_base_features,
    features = selected_features,
    metrics = metrics,
    objectives = objective_df,
    control = c(
      control,
      list(
        regularized = regularized,
        regularized_alpha = if (regularized) regularized_alpha else NULL,
        cv = cv,
        cv_folds = if (cv && !regularized) cv_folds else NULL,
        fitted_solution_id = selected_solution_id
      )
    ),
    training_data = training_data,
    model = model
  )

  panel
}

#' Resolve a Pareto Solution into Training-Ready Inputs
#'
#' Shared front half of [fit_panel()] and [fit_np_panel()]: chooses the
#' solution (explicit `features`, explicit `solution_id`, or auto-select on
#' the first objective), resolves training data from arguments or the stored
#' slots, validates the base features, and applies the stored feature
#' transform.
#'
#' @param optimization_result An `OptimizationResult`.
#' @param solution_id Optional solution ID; auto-selected when `NULL`.
#' @param features Optional explicit base features (overrides `solution_id`).
#' @param x,y,cohort Optional training data; defaults to the aggregated data
#'   stored on `optimization_result`.
#' @return List with `base_features`, `features` (transformed names),
#'   `solution_id` (`NA` for explicit features), `metrics` (named numeric or
#'   `NULL`), `x_raw`, `x_selected` (transform applied), `truth`, and `cohort`.
#' @keywords internal
.resolve_panel_solution <- function(optimization_result,
                                    solution_id = NULL,
                                    features = NULL,
                                    x = NULL,
                                    y = NULL,
                                    cohort = NULL) {
  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult from optimize_panel().",
         call. = FALSE)
  }

  solutions_df <- optimization_result@solutions
  if (nrow(solutions_df) == 0L) {
    stop("OptimizationResult contains no solutions.", call. = FALSE)
  }
  objective_cols <- setdiff(names(solutions_df),
                            c("solution_id", "base_features", "features"))

  # Determine which features to use
  if (!is.null(features)) {
    # User provided explicit base features
    selected_base_features <- features
    if (!all(features %in% optimization_result@feature_pool)) {
      missing <- setdiff(features, optimization_result@feature_pool)
      stop("Feature(s) not in feature pool: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    selected_solution_id <- NA_integer_
    solution_metrics <- NULL
  } else {
    # Select from Pareto solutions
    if (is.null(solution_id)) {
      # Auto-select: best on first objective
      if (length(objective_cols) == 0L) {
        stop("No objective columns found in solutions.", call. = FALSE)
      }
      first_obj <- objective_cols[1]
      direction <- optimization_result@control$objective_directions[[first_obj]] %||%
        "maximize"
      best_idx <- if (direction == "maximize") {
        which.max(solutions_df[[first_obj]])
      } else {
        which.min(solutions_df[[first_obj]])
      }
      solution_id <- solutions_df$solution_id[best_idx]
    }

    if (!solution_id %in% solutions_df$solution_id) {
      stop("solution_id ", solution_id, " not found. Valid IDs: ",
           paste(solutions_df$solution_id, collapse = ", "), call. = FALSE)
    }

    row_idx <- which(solutions_df$solution_id == solution_id)
    selected_base_features <- solutions_df$base_features[[row_idx]]
    selected_solution_id <- solution_id
    solution_metrics <- as.numeric(solutions_df[row_idx, objective_cols])
    names(solution_metrics) <- objective_cols
  }

  # Resolve training data (raw/untransformed base features)
  x_raw <- if (is.null(x)) optimization_result@aggregated_x else as.matrix(x)
  truth <- if (is.null(y)) {
    optimization_result@aggregated_y
  } else {
    ensure_binary_response(y)
  }
  cohort_vec <- if (is.null(cohort)) {
    optimization_result@aggregated_cohort
  } else {
    factor(cohort)
  }

  if (is.null(x_raw) || is.null(truth)) {
    stop("Training data not available. Provide x and y arguments.", call. = FALSE)
  }
  if (!all(selected_base_features %in% colnames(x_raw))) {
    missing <- setdiff(selected_base_features, colnames(x_raw))
    stop("Selected base feature(s) not found in training data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Apply the stored feature transform to the selected base features
  x_base_selected <- x_raw[, selected_base_features, drop = FALSE]
  feature_transform <- optimization_result@control$feature_transform %||% "none"

  if (feature_transform != "none" && length(selected_base_features) >= 2L) {
    x_selected <- .apply_feature_transform_single(x_base_selected, feature_transform)
    selected_features <- colnames(x_selected)
  } else {
    x_selected <- x_base_selected
    selected_features <- selected_base_features
  }

  list(
    base_features = selected_base_features,
    features = selected_features,
    solution_id = selected_solution_id,
    metrics = solution_metrics,
    x_raw = x_raw,
    x_selected = x_selected,
    truth = truth,
    cohort = cohort_vec
  )
}
