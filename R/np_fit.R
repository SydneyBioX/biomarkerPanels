#' Neyman-Pearson Panel Fitting
#'
#' Fits a biomarker panel using Neyman-Pearson (NP) penalized classification
#' via the nproc package. NP classification controls either the false positive
#' rate (FPR) or false negative rate (FNR) at a specified level alpha.
#'
#' @name np_fit
NULL

#' Fit a Neyman-Pearson Panel from Optimization Results
#'
#' Takes an `OptimizationResult` from [optimize_panel()] and fits a Neyman-Pearson
#' classifier on a selected Pareto solution. Returns a `TransferablePanelResult`
#' with the NP classifier and threshold.
#'
#' The NP paradigm controls one type of error (FPR or FNR) at level alpha with
#' probability at least 1-delta, while minimizing the other error type.
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
#' @param minimize_FPR Logical; if `TRUE` (default), the NP classifier controls
#'   the false positive rate (Type I error). If `FALSE`, labels are flipped
#'   internally so the classifier controls the false negative rate (miss rate),
#'   which is appropriate for rule-out diagnostic tests requiring high sensitivity.
#' @param alpha Target error rate to control (default 0.05). When `minimize_FPR = TRUE`,
#'   this is the maximum FPR. When `minimize_FPR = FALSE`, this is the maximum FNR
#'   (i.e., minimum sensitivity = 1 - alpha).
#' @param delta Tolerance for violating the alpha constraint (default 0.05).
#'   The error rate is guaranteed to be at most alpha with probability 1-delta.
#' @param method NP classification method passed to [nproc::npc()]. One of:
#'   "penlog" (penalized logistic, default), "logistic", "svm", "lda", "slda",
#'   "nb", "nnb", "ada", "tree", "randomforest".
#' @param split Number of random train/test splits for threshold estimation (default 1).
#' @param split.ratio Proportion of data used for training in each split (default 0.5).
#' @param n.cores Number of cores for parallel processing (default 1).
#' @param ... Additional arguments passed to [nproc::npc()].
#' @return A `TransferablePanelResult` with the fitted NP classifier.
#' @export
#' @seealso [optimize_panel()], [fit_panel()], [evaluate_panel()]
#' @examples
#' \dontrun{
#' # Run optimization targeting the high-sensitivity region
#' opt <- optimize_panel(
#'   x, y,
#'   objectives = define_objectives(
#'     metrics = c("pauc", "num_features"),
#'     params = list(pauc = list(sens_floor = 0.90))
#'   ),
#'   constraints = list(min_metric_constraint("sensitivity", 0.90))
#' )
#'
#' # Fit NP classifier for rule-out test (control FNR to achieve high sensitivity)
#' panel <- fit_np_panel(opt, minimize_FPR = FALSE, alpha = 0.05)
#'
#' # This guarantees sensitivity >= 95% with high probability
#' per_cohort_metrics(panel)
#'
#' # Evaluate on held-out data
#' evaluate_panel(panel, x_test, y_test)
#' }
fit_np_panel <- function(optimization_result,
                         solution_id = NULL,
                         features = NULL,
                         x = NULL,
                         y = NULL,
                         cohort = NULL,
                         minimize_FPR = TRUE,
                         alpha = 0.05,
                         delta = 0.05,
                         method = "penlog",
                         split = 1L,
                         split.ratio = 0.5,
                         n.cores = 1L,
                         ...) {
  # Validate inputs
  .check_nproc_available()
  .validate_np_params(alpha, delta, method)

  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult from optimize_panel().",
      call. = FALSE
    )
  }

  solutions_df <- optimization_result@solutions
  if (nrow(solutions_df) == 0L) {
    stop("OptimizationResult contains no solutions.", call. = FALSE)
  }

  # Determine which features to use (same pattern as fit_panel)
  if (!is.null(features)) {
    # User provided explicit base features
    selected_base_features <- features
    if (!all(features %in% optimization_result@feature_pool)) {
      missing <- setdiff(features, optimization_result@feature_pool)
      stop("Feature(s) not in feature pool: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    selected_solution_id <- NA_integer_
    solution_metrics <- NULL
  } else {
    if (is.null(solution_id)) {
      # Auto-select: best on first objective
      objective_cols <- setdiff(names(solutions_df), c("solution_id", "base_features", "features"))
      if (length(objective_cols) == 0L) {
        stop("No objective columns found in solutions.", call. = FALSE)
      }
      first_obj <- objective_cols[1]
      control <- optimization_result@control
      direction <- if (!is.null(control$objective_directions)) {
        control$objective_directions[[first_obj]]
      } else {
        "maximize"
      }
      if (direction == "maximize") {
        best_idx <- which.max(solutions_df[[first_obj]])
      } else {
        best_idx <- which.min(solutions_df[[first_obj]])
      }
      solution_id <- solutions_df$solution_id[best_idx]
    }

    if (!solution_id %in% solutions_df$solution_id) {
      stop("solution_id ", solution_id, " not found. Valid IDs: ",
        paste(solutions_df$solution_id, collapse = ", "),
        call. = FALSE
      )
    }

    row_idx <- which(solutions_df$solution_id == solution_id)
    selected_base_features <- solutions_df$base_features[[row_idx]]
    selected_solution_id <- solution_id

    objective_cols <- setdiff(names(solutions_df), c("solution_id", "base_features", "features"))
    solution_metrics <- as.numeric(solutions_df[row_idx, objective_cols])
    names(solution_metrics) <- objective_cols
  }

  # Get training data (same pattern as fit_panel)
  if (is.null(x)) {
    x_mat <- optimization_result@aggregated_x
  } else {
    x_mat <- as.matrix(x)
  }

  if (is.null(y)) {
    truth <- optimization_result@aggregated_y
  } else {
    truth <- ensure_binary_response(y)
  }

  if (is.null(cohort)) {
    cohort_vec <- optimization_result@aggregated_cohort
  } else {
    cohort_vec <- factor(cohort)
  }

  # Validate data
  if (is.null(x_mat) || is.null(truth)) {
    stop("Training data not available. Provide x and y arguments.", call. = FALSE)
  }

  # Validate base features exist in training data
  if (!all(selected_base_features %in% colnames(x_mat))) {
    missing <- setdiff(selected_base_features, colnames(x_mat))
    stop("Selected base feature(s) not found in training data: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Extract selected base features
  x_base_selected <- x_mat[, selected_base_features, drop = FALSE]

  # Apply feature transform (same pattern as fit_panel)
  feature_transform <- optimization_result@control$feature_transform
  if (is.null(feature_transform)) {
    feature_transform <- "none"  # Backward compatibility
  }

  if (feature_transform != "none" && length(selected_base_features) >= 2L) {
    x_selected <- .apply_feature_transform_single(x_base_selected, feature_transform)
    selected_features <- colnames(x_selected)
  } else {
    x_selected <- x_base_selected
    selected_features <- selected_base_features
  }

  # Prepare labels for NP classification
  label_info <- .prepare_np_labels(truth, minimize_FPR)
  y_nproc <- label_info$y_nproc
  labels_flipped <- label_info$labels_flipped

  # Fit NP classifier
  np_model <- tryCatch(
    {
      nproc::npc(
        x = as.matrix(x_selected),
        y = y_nproc,
        method = method,
        alpha = alpha,
        delta = delta,
        split = split,
        split.ratio = split.ratio,
        n.cores = n.cores,
        ...
      )
    },
    error = function(e) {
      stop(
        "NP classifier fitting failed: ", conditionMessage(e),
        "\nCheck that nproc is installed correctly and data is appropriate for the method.",
        call. = FALSE
      )
    }
  )

  # Store metadata for prediction
  np_model$biomarkerPanels_meta <- list(
    feature_names = selected_features,
    labels_flipped = labels_flipped,
    minimize_FPR = minimize_FPR,
    alpha = alpha,
    delta = delta,
    method = method
  )

  # Get predictions on training data
  np_pred <- stats::predict(np_model, newx = as.matrix(x_selected))
  raw_scores <- if (is.list(np_pred) && "pred.score" %in% names(np_pred)) {
    as.numeric(np_pred$pred.score)
  } else {
    as.numeric(np_pred)
  }

  # Invert predictions if labels were flipped
  scores <- .invert_np_predictions(raw_scores, labels_flipped)

  # Extract or compute threshold.
  # An npc object has NO top-level $cutoff; the NP-calibrated cutoff lives in
  # fits[[i]]$cutoff (one per split). The cutoff applies to the *signed*
  # pred.score returned by predict.npc. For the default penlog/logistic case
  # sign == TRUE and pred.score is P(class 1) in [0, 1], so:
  #   - minimize_FPR (no flip): predict Yes iff score >= cutoff  -> threshold = cutoff
  #   - flipped:                scores = 1 - P(No); predict Yes iff scores >= 1 - cutoff
  # When sign == FALSE the score is negated (not in [0, 1]); the flip arithmetic
  # below would be invalid, so we fall back. With split > 1 predict.npc averages
  # scores across splits and there is no single scalar cutoff matching the
  # averaged score, so a scalar NP threshold is not well defined either.
  np_fit1 <- if (length(np_model$fits) >= 1L) np_model$fits[[1L]] else NULL
  np_cutoff <- np_fit1$cutoff
  np_sign <- np_fit1$sign

  # nproc returns a length-1 cutoff per (alpha) request, or NA when the sample
  # is too small to meet the (alpha, delta) guarantee.
  np_cutoff_usable <- !is.null(np_cutoff) && length(np_cutoff) == 1L &&
    is.finite(np_cutoff)

  if (np_cutoff_usable && split > 1L) {
    warning(
      "split > 1: npc averages prediction scores across splits, so there is ",
      "no single NP-calibrated cutoff. Falling back to an empirical ",
      "training-data threshold, which does NOT carry the (alpha, delta) ",
      "guarantee. Use split = 1 to retain the NP threshold.",
      call. = FALSE
    )
  } else if (!is.null(np_cutoff) && !np_cutoff_usable) {
    warning(
      "nproc could not calibrate an NP cutoff for alpha = ", alpha,
      ", delta = ", delta, " (sample likely too small). Falling back to an ",
      "empirical training-data threshold, which does NOT carry the ",
      "(alpha, delta) guarantee.",
      call. = FALSE
    )
  }

  use_np_cutoff <- np_cutoff_usable && split == 1L &&
    (!labels_flipped || isTRUE(np_sign))

  if (use_np_cutoff) {
    threshold <- if (labels_flipped) {
      1 - np_cutoff
    } else {
      np_cutoff
    }
  } else {
    # Fallback: compute threshold from target specificity (for minimize_FPR=TRUE)
    # or target sensitivity (for minimize_FPR=FALSE)
    neg_scores <- scores[truth == "No"]
    pos_scores <- scores[truth == "Yes"]
    if (minimize_FPR) {
      target_spec <- 1 - alpha
      threshold <- if (length(neg_scores) > 0) {
        stats::quantile(neg_scores, probs = target_spec, na.rm = TRUE)
      } else {
        0.5
      }
    } else {
      target_sens <- 1 - alpha
      threshold <- if (length(pos_scores) > 0) {
        stats::quantile(pos_scores, probs = 1 - target_sens, na.rm = TRUE)
      } else {
        0.5
      }
    }
  }
  threshold <- as.numeric(threshold)

  # Compute per-cohort metrics
  if (is.null(cohort_vec)) {
    cohort_vec <- factor(rep("cohort_01", nrow(x_selected)))
  }
  per_cohort <- .compute_per_cohort_metrics(scores, truth, cohort_vec, threshold)

  # Compute weighted variance
  weighted_var <- .compute_weighted_variance(per_cohort)

  # Build metrics from solution or compute from training
  if (!is.null(solution_metrics)) {
    metrics <- solution_metrics
  } else {
    # Compute basic metrics on training data
    sens <- sum(scores >= threshold & truth == "Yes") / sum(truth == "Yes")
    spec <- sum(scores < threshold & truth == "No") / sum(truth == "No")
    metrics <- c(sensitivity = sens, specificity = spec)
  }

  # Build objectives data frame
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

  # Build control list
  control <- optimization_result@control
  training_data <- list(
    n = nrow(x_mat),
    p = ncol(x_mat),
    class_balance = table(truth),
    feature_pool_size = length(optimization_result@feature_pool),
    num_cohorts = if (!is.null(cohort_vec)) length(unique(cohort_vec)) else 1L
  )

  # Create TransferablePanelResult
  new(
    "TransferablePanelResult",
    base_features = selected_base_features,
    features = selected_features,
    metrics = metrics,
    objectives = objective_df,
    control = c(
      control,
      list(
        np_model_type = "npc",
        method = method,
        minimize_FPR = minimize_FPR,
        labels_flipped = labels_flipped,
        fitted_solution_id = selected_solution_id
      )
    ),
    training_data = training_data,
    model = np_model,
    np_threshold = threshold,
    np_alpha = alpha,
    np_delta = delta,
    per_cohort_metrics = per_cohort,
    weighted_variance = weighted_var,
    validation_metrics = list(),
    partition_info = list()
  )
}
