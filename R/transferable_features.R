#' Transferable Feature Selection
#'
#' Functions for selecting features based on ridge regression transferability
#' across cohorts, emphasizing biomarkers with consistent and strong effects.
#'
#' @name transferable_features
NULL

#' Select Ridge-Stable Transferable Features
#'
#' Fits ridge-penalised logistic regression models for each cohort and selects
#' features whose coefficients are both large in magnitude and consistent across
#' cohorts. This emphasises biomarkers that deliver transferable signal.
#'
#' @param x_list A matrix-like object, `SummarizedExperiment`, or list of such
#'   objects. Rows represent samples and columns represent features.
#' @param y_list A binary response aligned with `x_list`. Must be a factor with
#'   levels `c("No", "Yes")`. Provide a list when `x_list` is a list.
#' @param n_features Maximum number of features to return (default `50`).
#' @param lambda Optional penalty parameter supplied to [glmnet::glmnet()]. When
#'   `NULL`, cohort-specific values are chosen via cross-validation.
#' @param lambda_choice When `lambda` is `NULL`, choose either `"lambda_1se"`
#'   (default) or `"lambda_min"` from [glmnet::cv.glmnet()] as the operating
#'   point for coefficient extraction.
#' @param min_coefficient Minimum absolute coefficient required (applied to the
#'   cohort-wise minimum magnitude) before a feature is considered transferable.
#' @param require_sign_consistency If `TRUE`, only retain features whose ridge
#'   coefficients meet the `sign_consistency_threshold` across cohorts. Note
#'   that `FALSE` is equivalent to `sign_consistency_threshold = 0`: sign
#'   agreement is a non-negative fraction, so a zero threshold never filters.
#' @param sign_consistency_threshold Minimum fraction of cohorts that must agree
#'   on coefficient sign direction (default 1.0 requires 100% agreement). Set to
#'   values like 0.75 to allow partial sign consistency.
#' @param standardize Logical; passed to [glmnet::glmnet()] to control feature
#'   standardisation prior to fitting. Defaults to `TRUE`.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index to
#'   extract prior to modelling.
#' @return A list containing the selected feature identifiers, per-feature
#'   scoring metadata, coefficient matrix, and the lambda value used for each
#'   cohort.
#' @details Features are aligned across cohorts via intersection of shared
#'   column names. For more flexible alignment strategies (e.g., majority
#'   voting or median imputation), use the `feature_alignment` parameter
#'   in [optimize_panel()].
#' @export
#'
#' @importFrom glmnet cv.glmnet glmnet
select_transferable_features <- function(x_list,
                                         y_list,
                                         n_features = 50L,
                                         lambda = NULL,
                                         lambda_choice = c("lambda_1se", "lambda_min"),
                                         min_coefficient = 0,
                                         require_sign_consistency = TRUE,
                                         sign_consistency_threshold = 1.0,
                                         standardize = TRUE,
                                         assay = NULL) {
  n_features <- .validate_positive_integer(n_features, "n_features")
  lambda_choice <- match.arg(lambda_choice)
  if (!is.numeric(min_coefficient) || length(min_coefficient) != 1L) {
    stop("`min_coefficient` must be a numeric scalar.", call. = FALSE)
  }

  prepared <- .prepare_selection_inputs(x_list, y_list, assay = assay)
  matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  feature_names <- prepared$feature_names

  if (!length(feature_names)) {
    stop(
      "No common features found across cohorts for ridge regression. ",
      "Check that cohorts have overlapping feature names.",
      call. = FALSE
    )
  }

  lambda_vec <- .resolve_lambda_vector(lambda, length(matrices))

  ridge_results <- lapply(seq_along(matrices), function(i) {
    x_mat <- matrices[[i]]
    y_vec <- responses[[i]]

    if (is.null(lambda_vec)) {
      cv_fit <- glmnet::cv.glmnet(
        x = x_mat,
        y = y_vec,
        alpha = 0,
        family = "binomial",
        standardize = standardize
      )
      chosen_lambda <- if (lambda_choice == "lambda_min") {
        cv_fit$lambda.min
      } else {
        cv_fit$lambda.1se
      }

      coef_vec <- .extract_glmnet_coefficients(cv_fit, chosen_lambda, feature_names)
      list(coefficients = coef_vec, lambda = chosen_lambda)
    } else {
      fit <- glmnet::glmnet(
        x = x_mat,
        y = y_vec,
        alpha = 0,
        family = "binomial",
        standardize = standardize,
        lambda = lambda_vec[[i]]
      )
      coef_vec <- .extract_glmnet_coefficients(fit, lambda_vec[[i]], feature_names)
      list(coefficients = coef_vec, lambda = lambda_vec[[i]])
    }
  })

  coefficient_matrix <- do.call(cbind, lapply(ridge_results, `[[`, "coefficients"))
  rownames(coefficient_matrix) <- feature_names
  colnames(coefficient_matrix) <- cohort_names

  lambda_used <- vapply(ridge_results, `[[`, numeric(1), "lambda")
  names(lambda_used) <- cohort_names

  combined_scores <- .score_transferable_features(
    coefficient_matrix,
    min_coefficient = min_coefficient,
    require_sign_consistency = require_sign_consistency,
    sign_consistency_threshold = sign_consistency_threshold
  )

  settings <- list(
    lambda_choice = if (is.null(lambda_vec)) lambda_choice else "fixed",
    min_coefficient = min_coefficient,
    require_sign_consistency = require_sign_consistency,
    sign_consistency_threshold = sign_consistency_threshold,
    standardize = standardize
  )

  if (!nrow(combined_scores)) {
    warning(
      "No transferable features passed the sign-consistency and minimum-",
      "coefficient filters; returning an empty feature set. Consider relaxing ",
      "`min_coefficient`, `sign_consistency_threshold`, or setting ",
      "`require_sign_consistency = FALSE`.",
      call. = FALSE
    )
    return(list(
      features = character(),
      scores = combined_scores,
      coefficients = coefficient_matrix[NULL, , drop = FALSE],
      lambda = lambda_used,
      settings = settings
    ))
  }

  top_n <- min(n_features, nrow(combined_scores))
  top_indices <- seq_len(top_n)
  selected_scores <- combined_scores[top_indices, , drop = FALSE]
  selected_features <- selected_scores$feature
  selected_row_index <- selected_scores$row_index
  coefficients_subset <- coefficient_matrix[selected_row_index, , drop = FALSE]
  rownames(coefficients_subset) <- selected_features
  selected_scores$row_index <- NULL

  list(
    features = selected_features,
    scores = selected_scores,
    coefficients = coefficients_subset,
    lambda = lambda_used,
    settings = settings
  )
}

#' Resolve Lambda Parameter Vector
#'
#' Validates and expands the lambda parameter for ridge regression, handling
#' both single values and cohort-specific vectors.
#'
#' @param lambda User-supplied lambda or `NULL` for cross-validation.
#' @param num_cohorts Number of cohorts.
#' @return A numeric vector of length `num_cohorts`, or `NULL`.
#' @keywords internal
.resolve_lambda_vector <- function(lambda, num_cohorts) {
  if (is.null(lambda)) {
    return(NULL)
  }

  if (!is.numeric(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("`lambda` must be positive numeric.", call. = FALSE)
  }
  if (length(lambda) == 1L) {
    rep(lambda, num_cohorts)
  } else if (length(lambda) == num_cohorts) {
    lambda
  } else {
    stop("`lambda` must be length 1 or match the number of cohorts.", call. = FALSE)
  }
}

#' Extract Coefficients from Glmnet Model
#'
#' Extracts feature coefficients (excluding intercept) from a fitted glmnet
#' model at the specified lambda value.
#'
#' @param model A fitted glmnet or cv.glmnet object.
#' @param lambda The lambda value at which to extract coefficients.
#' @param feature_names Character vector of feature names to return.
#' @return A named numeric vector of coefficients.
#' @keywords internal
.extract_glmnet_coefficients <- function(model, lambda, feature_names) {
  coef_mat <- stats::coef(model, s = lambda)
  coef_vec <- as.numeric(coef_mat)[-1]
  names(coef_vec) <- rownames(coef_mat)[-1]
  coef_vec[feature_names]
}

#' Score Features by Transferability Across Cohorts
#'
#' Computes feature scores based on mean coefficient magnitude, consistency
#' across cohorts, and sign agreement. Uses C++ for performance.
#'
#' @param coefficient_matrix A matrix of ridge coefficients (features x cohorts).
#' @param min_coefficient Minimum absolute coefficient required.
#' @param require_sign_consistency Whether to filter by sign consistency.
#'   Redundant with `sign_consistency_threshold`: `FALSE` behaves exactly like
#'   `sign_consistency_threshold = 0`, because `sign_agreement >= 0` is always
#'   true.
#' @param sign_consistency_threshold Minimum fraction of cohorts agreeing on
#'   coefficient sign.
#' @return A data.frame with feature scores, sorted by decreasing score.
#' @keywords internal
.score_transferable_features <- function(coefficient_matrix,
                                         min_coefficient,
                                         require_sign_consistency,
                                         sign_consistency_threshold = 1.0) {
  if (!nrow(coefficient_matrix)) {
    return(data.frame())
  }

  feature_names <- rownames(coefficient_matrix)
  row_index <- seq_len(nrow(coefficient_matrix))

  # Call C++ implementation for the heavy computation
  cpp_result <- .score_transferable_features_cpp(coefficient_matrix)

  mean_abs <- cpp_result$mean_abs
  sd_coeff <- cpp_result$sd
  min_abs <- cpp_result$min_abs
  sign_agreement <- cpp_result$sign_agreement

  sign_consistent <- sign_agreement >= sign_consistency_threshold
  score <- mean_abs / (sd_coeff + 1e-6)

  keep <- is.finite(score) & !is.na(score) & (min_abs >= min_coefficient)
  # Redundant with the threshold: `require_sign_consistency = FALSE` gives the
  # same result as `sign_consistency_threshold = 0`, since sign_agreement is a
  # non-negative fraction and so always clears a zero threshold.
  if (require_sign_consistency) {
    keep <- keep & sign_consistent
  }

  score_df <- data.frame(
    feature = feature_names,
    row_index = row_index,
    mean_abs = mean_abs,
    sd = sd_coeff,
    min_abs = min_abs,
    sign_agreement = sign_agreement,
    sign_consistent = sign_consistent,
    score = score,
    stringsAsFactors = FALSE
  )

  score_df <- score_df[keep, , drop = FALSE]
  score_df <- score_df[order(score_df$score, decreasing = TRUE), , drop = FALSE]
  rownames(score_df) <- NULL
  score_df
}

