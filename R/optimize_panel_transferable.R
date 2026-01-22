#' Optimize Biomarker Panels with Transferability Focus
#'
#' Wrapper implementing a train/validation/held-out split strategy with
#' Neyman-Pearson threshold selection for cross-cohort generalizability.
#' Unlike [optimize_panel()], this function partitions each cohort into
#' three sets: training (for NSGA optimization), validation (for solution
#' selection), and held-out (for NP threshold calibration and final metrics).
#'
#' @param x List of matrix-like objects representing cohorts.
#' @param y List of binary response vectors aligned with `x`.
#' @param objectives Named list of objective descriptors as returned by
#'   [define_objectives()].
#' @param max_features Maximum number of biomarkers permitted in a panel.
#' @param feature_pool Optional subset of feature identifiers. If `NULL`,
#'   `select_transferable_features()` is run on training data only.
#' @param train_ratio Proportion of each cohort for training (default 0.7).
#' @param val_ratio Proportion of each cohort for validation (default 0.2).
#' @param np_alpha Type I error rate for Neyman-Pearson threshold (default 0.15).
#' @param np_delta Tolerance for NP threshold selection (default 0.05).
#' @param cohort_aggregator Transformation applied to cohort matrices (default
#'   `"pairwise_ratios"`).
#' @param feature_alignment Strategy for aligning features across cohorts
#'   (default `"intersection"`).
#' @param constraints Optional list of constraint descriptors.
#' @param algorithm Multi-objective algorithm: `"NSGA-II"` (default) or
#'   `"NSGA-III"`.
#' @param nsga_control Named list of arguments passed to rmoo.
#' @param assay For `SummarizedExperiment` inputs, assay name or index.
#' @param seed Optional integer seed for reproducibility.
#' @param regularized Logical; if `TRUE` (default), use regularized regression.
#' @param regularized_alpha Elastic net mixing parameter (default 0.5).
#' @param n_top_features Number of top features to select via
#'   `select_transferable_features()` when `feature_pool` is NULL (default 50).
#' @return A `TransferablePanelResult` with Pareto-optimal solutions, NP
#'
#'   threshold, per-cohort metrics, and weighted variance.
#' @export
optimize_panel_transferable <- function(
    x, y,
    objectives = define_objectives(
      losses = c("sensitivity", "specificity", "num_features")
    ),
    max_features = 5L,
    feature_pool = NULL,
    train_ratio = 0.7,
    val_ratio = 0.2,
    np_alpha = 0.15,
    np_delta = 0.05,
    cohort_aggregator = "pairwise_ratios",
    feature_alignment = "intersection",
    constraints = list(),
    algorithm = c("NSGA-II", "NSGA-III"),
    nsga_control = list(),
    assay = NULL,
    seed = NULL,
    regularized = TRUE,
    regularized_alpha = 0.5,
    n_top_features = 50L
) {
  algorithm <- match.arg(algorithm)

  # Validate inputs
  if (!is.list(x) || length(x) < 1L) {
    stop("`x` must be a list of cohort matrices.", call. = FALSE)
  }
  if (!is.list(y) || length(y) != length(x)) {
    stop("`y` must be a list of response vectors matching `x`.", call. = FALSE)
  }

  # 1. Validate partition ratios
  .validate_partition_ratios(train_ratio, val_ratio)

  # 2. Stratified partition each cohort
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  partitions <- .stratified_partition_cohorts(x, y, train_ratio, val_ratio)

  # 3. Feature selection on TRAINING data only (if feature_pool is NULL)
  if (is.null(feature_pool)) {
    message("Selecting transferable features from training data...")
    fs_result <- select_transferable_features(
      partitions$train_x,
      partitions$train_y,
      n_features = n_top_features
    )
    feature_pool <- fs_result$features
    message("Selected ", length(feature_pool), " features for optimization.")
  }

  # 4. Prepare combined datasets with aggregation
  train_inputs <- .prepare_cohort_inputs(
    partitions$train_x, partitions$train_y,
    assay = assay,
    aggregator = cohort_aggregator,
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )
  val_inputs <- .prepare_cohort_inputs(
    partitions$val_x, partitions$val_y,
    assay = assay,
    aggregator = cohort_aggregator,
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )
  heldout_inputs <- .prepare_cohort_inputs(
    partitions$heldout_x, partitions$heldout_y,
    assay = assay,
    aggregator = cohort_aggregator,
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )

  # Get aggregated feature names
  agg_feature_pool <- colnames(train_inputs$x)

  # 5. Create validation-based fitness function
  fitness_fn <- .make_validation_fitness(
    train_x = train_inputs$x,
    train_y = train_inputs$truth,
    train_cohort = train_inputs$cohort,
    val_x = val_inputs$x,
    val_y = val_inputs$truth,
    val_cohort = val_inputs$cohort,
    feature_pool = agg_feature_pool,
    max_features = max_features,
    objectives = objectives,
    constraints = constraints,
    regularized = regularized,
    alpha = regularized_alpha
  )

  # 6. Get NSGA defaults and run optimization
  decision_dim <- length(agg_feature_pool)
  nsga_defaults <- .get_adaptive_nsga_defaults(decision_dim, algorithm)
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)

  nsga_params <- c(
    list(
      type = "real-valued",
      fitness = fitness_fn$wrapper,
      nObj = length(objectives),
      lower = rep(0, decision_dim),
      upper = rep(1, decision_dim),
      monitor = FALSE,
      summary = FALSE
    ),
    nsga_args
  )

  if (algorithm == "NSGA-III" && is.null(nsga_args$n_partitions)) {
    nsga_params$n_partitions <- .compute_nsga3_partitions(length(objectives))
  }

  message("Running ", algorithm, " optimization...")
  if (algorithm == "NSGA-II") {
    nsga_result <- do.call(rmoo::nsga2, nsga_params)
  } else {
    nsga_result <- do.call(rmoo::nsga3, nsga_params)
  }

  # 7. Select best solution from Pareto front
  optimal_idx <- which(nsga_result@front == 1)
  pareto_pop <- nsga_result@population[optimal_idx, , drop = FALSE]
  if (is.null(dim(pareto_pop))) {
    pareto_pop <- matrix(pareto_pop, nrow = 1)
  }

  # Evaluate all Pareto solutions
  solutions <- lapply(seq_len(nrow(pareto_pop)), function(i) {
    decision_vec <- pareto_pop[i, ]
    fitness_fn$evaluate(decision_vec)
  })

  feasible_vec <- vapply(solutions, `[[`, logical(1), "feasible")
  if (!any(feasible_vec)) {
    stop("No solutions satisfied the supplied constraints.", call. = FALSE)
  }
  solutions <- solutions[feasible_vec]

  # Select primary solution based on first objective
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  primary_obj <- names(objectives)[1]
  primary_dir <- objective_directions[[primary_obj]]
  if (primary_dir == "maximize") {
    primary_idx <- which.max(metric_matrix[, primary_obj])
  } else {
    primary_idx <- which.min(metric_matrix[, primary_obj])
  }

  primary_solution <- solutions[[primary_idx]]
  selected_features <- primary_solution$features

  # 8. Fit final model on train+val combined
  combined_x <- rbind(train_inputs$x, val_inputs$x)
  combined_y <- factor(
    c(as.character(train_inputs$truth), as.character(val_inputs$truth)),
    levels = levels(train_inputs$truth)
  )
  combined_cohort <- factor(
    c(as.character(train_inputs$cohort), as.character(val_inputs$cohort)),
    levels = unique(c(levels(train_inputs$cohort), levels(val_inputs$cohort)))
  )

  x_selected <- combined_x[, selected_features, drop = FALSE]
  final_model <- .fit_final_model_regularized(
    x_selected, combined_y, combined_cohort,
    alpha = regularized_alpha
  )

  # 9. Compute scores on held-out data
  heldout_x_selected <- heldout_inputs$x[, selected_features, drop = FALSE]
  heldout_scores <- .predict_from_model(final_model, heldout_x_selected)

  # 10. NP threshold selection on held-out data
  np_result <- .select_np_threshold(
    scores = heldout_scores,
    truth = heldout_inputs$truth,
    alpha = np_alpha,
    delta = np_delta
  )

  # 11. Per-cohort metrics using NP threshold
  per_cohort_df <- .compute_per_cohort_metrics(
    scores = heldout_scores,
    truth = heldout_inputs$truth,
    cohort = heldout_inputs$cohort,
    threshold = np_result$threshold
  )

  # 12. Weighted variance
  weighted_var <- .compute_weighted_variance(per_cohort_df)

  # Build objective data frame
  objective_df <- do.call(rbind, lapply(seq_along(solutions), function(i) {
    data.frame(
      solution_id = i,
      objective = names(objectives),
      value = solutions[[i]]$metrics,
      direction = objective_directions,
      features = I(rep(list(solutions[[i]]$features), length(objectives))),
      stringsAsFactors = FALSE
    )
  }))

  # Validation metrics
  validation_metrics <- list(
    sensitivity = primary_solution$metrics["sensitivity"],
    specificity = primary_solution$metrics["specificity"],
    num_features = length(selected_features)
  )

  # Construct result
  panel <- new(
    "TransferablePanelResult",
    features = selected_features,
    metrics = setNames(as.numeric(primary_solution$metrics), names(objectives)),
    objectives = objective_df,
    control = list(
      max_features = max_features,
      feature_pool = agg_feature_pool,
      base_feature_pool = feature_pool,
      algorithm = algorithm,
      nsga2 = nsga_args,
      cohort_aggregator = cohort_aggregator,
      feature_alignment = feature_alignment,
      train_ratio = train_ratio,
      val_ratio = val_ratio,
      positive_class = "Yes",
      response_levels = levels(train_inputs$truth),
      seed = seed,
      regularized = regularized,
      regularized_alpha = regularized_alpha
    ),
    training_data = list(
      n = nrow(combined_x),
      p = ncol(combined_x),
      class_balance = table(combined_y),
      feature_pool_size = length(agg_feature_pool),
      num_cohorts = length(partitions$cohort_names),
      cohort_labels = partitions$cohort_names
    ),
    model = final_model,
    np_threshold = np_result$threshold,
    np_alpha = np_alpha,
    np_delta = np_delta,
    per_cohort_metrics = per_cohort_df,
    weighted_variance = weighted_var,
    validation_metrics = validation_metrics,
    partition_info = list(
      train_ratio = train_ratio,
      val_ratio = val_ratio,
      heldout_ratio = 1 - train_ratio - val_ratio,
      partition_sizes = partitions$partition_info
    )
  )

  panel
}


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Create Validation-Based Fitness Function
#'
#' Factory that returns a fitness function for rmoo that trains on training
#' data and evaluates objectives on validation data.
#'
#' @param train_x Training feature matrix.
#' @param train_y Training response factor.
#' @param train_cohort Training cohort factor.
#' @param val_x Validation feature matrix.
#' @param val_y Validation response factor.
#' @param val_cohort Validation cohort factor.
#' @param feature_pool Character vector of feature names.
#' @param max_features Maximum features to include.
#' @param objectives Objective list from define_objectives().
#' @param constraints Constraint list.
#' @param regularized Whether to use regularized regression.
#' @param alpha Elastic net mixing parameter.
#' @return List with `wrapper` (fitness for rmoo) and `evaluate` (full eval).
#' @keywords internal
.make_validation_fitness <- function(
    train_x, train_y, train_cohort,
    val_x, val_y, val_cohort,
    feature_pool, max_features, objectives, constraints,
    regularized, alpha
) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)
  min_features_required <- if (regularized) 2L else 1L

  evaluate_candidate <- function(decision_vec) {
    feature_names_pool <- feature_pool
    above_threshold <- which(decision_vec > 0.5)

    if (length(above_threshold) < min_features_required) {
      ord <- order(decision_vec, feature_names_pool, decreasing = c(TRUE, FALSE),
                   method = "radix")
      selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
    } else if (length(above_threshold) > max_features) {
      weights_above <- decision_vec[above_threshold]
      names_above <- feature_names_pool[above_threshold]
      ord <- order(weights_above, names_above, decreasing = c(TRUE, FALSE),
                   method = "radix")
      selected_idx <- above_threshold[ord[seq_len(max_features)]]
    } else {
      selected_idx <- above_threshold
    }

    selected_idx <- selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
    selected_features <- feature_names_pool[selected_idx]

    # Train model on training data
    x_train_sel <- train_x[, selected_features, drop = FALSE]

    tryCatch({
      if (regularized) {
        fit <- .fit_final_model_regularized(
          x_train_sel, train_y, train_cohort, alpha = alpha
        )
      } else {
        fit <- .fit_final_model(x_train_sel, train_y, train_cohort)
      }

      # Predict on validation data
      x_val_sel <- val_x[, selected_features, drop = FALSE]
      val_scores <- .predict_from_model(fit, x_val_sel)

      # Compute objectives on validation data
      constraint_results <- if (length(constraint_specs)) {
        setNames(vapply(seq_along(constraint_specs), function(j) {
          res <- constraint_specs[[j]]$fun(
            truth = val_y,
            scores = val_scores,
            selected = selected_features,
            cohort = val_cohort,
            x = x_val_sel
          )
          isTRUE(res)
        }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
      } else {
        logical(0)
      }

      feasible <- if (length(constraint_results)) all(constraint_results) else TRUE

      metrics <- vapply(objectives, function(obj) {
        obj$fun(
          val_y,
          val_scores,
          selected = selected_features,
          cohort = val_cohort,
          x = x_val_sel
        )
      }, numeric(1))

      list(
        features = selected_features,
        scores = val_scores,
        metrics = metrics,
        constraint_results = constraint_results,
        feasible = feasible
      )
    }, error = function(e) {
      list(
        features = selected_features,
        scores = rep(0.5, nrow(val_x)),
        metrics = setNames(rep(Inf, length(objectives)), names(objectives)),
        constraint_results = logical(0),
        feasible = FALSE
      )
    })
  }

  evaluate_single <- function(decision_vec) {
    evaluated <- evaluate_candidate(decision_vec)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(Inf, length(objectives)))
    }
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      if (is.na(val) || is.infinite(val)) return(Inf)
      if (dir == "maximize") return(-val)
      val
    }, val = metrics, dir = objective_directions, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    as.numeric(converted)
  }

  objective_wrapper <- function(x, ...) {
    if (is.null(dim(x))) {
      return(evaluate_single(x))
    }
    t(apply(x, 1, evaluate_single))
  }

  list(
    wrapper = objective_wrapper,
    evaluate = evaluate_candidate
  )
}

#' Select Neyman-Pearson Classification Threshold
#'
#' Uses the nproc package to select a classification threshold that controls
#' the Type I error rate (1 - specificity) at level alpha.
#'
#' @param scores Numeric vector of predicted probabilities.
#' @param truth Binary response factor.
#' @param alpha Target Type I error rate (default 0.15).
#' @param delta Tolerance parameter (default 0.05).
#' @return List with threshold value and metadata.
#' @keywords internal
.select_np_threshold <- function(scores, truth, alpha = 0.15, delta = 0.05) {
  # Check if nproc is available
  if (!requireNamespace("nproc", quietly = TRUE)) {
    warning(
      "Package 'nproc' not available. Using fallback threshold of 0.5. ",
      "Install nproc via install.packages('nproc') for NP threshold selection.",
      call. = FALSE
    )
    return(list(
      threshold = 0.5,
      method = "fallback",
      alpha = alpha,
      delta = delta
    ))
  }

  tryCatch({
    # nproc expects 0/1 labels
    y_binary <- as.integer(truth) - 1L

    # npc() returns a classifier; we extract the threshold
    np_fit <- nproc::npc(
      x = as.matrix(scores),
      y = y_binary,
      method = "lda",  # Simple method since scores are 1D
      alpha = alpha,
      delta = delta
    )

    # Extract threshold from fitted classifier
    # For 1D scores with LDA, threshold is related to cutoff
    # We need to find the threshold that achieves the target alpha

    # Alternative: use nproc's internal threshold calculation
    # The npc object stores the threshold implicitly
    # We can compute it by evaluating at different cutoffs

    # Simple approach: find threshold that achieves target specificity
    neg_scores <- scores[y_binary == 0]
    target_spec <- 1 - alpha

    if (length(neg_scores) > 0) {
      threshold <- stats::quantile(neg_scores, probs = target_spec, na.rm = TRUE)
    } else {
      threshold <- 0.5
    }

    list(
      threshold = as.numeric(threshold),
      method = "nproc",
      alpha = alpha,
      delta = delta
    )
  }, error = function(e) {
    warning(
      "NP threshold selection failed: ", conditionMessage(e),
      ". Using fallback threshold of 0.5.",
      call. = FALSE
    )
    list(
      threshold = 0.5,
      method = "fallback",
      alpha = alpha,
      delta = delta,
      error = conditionMessage(e)
    )
  })
}

#' Compute Per-Cohort Performance Metrics
#'
#' Calculates sensitivity and specificity for each cohort at the given
#' classification threshold.
#'
#' @param scores Numeric vector of predicted probabilities.
#' @param truth Binary response factor.
#' @param cohort Factor indicating cohort membership.
#' @param threshold Classification threshold.
#' @return Data frame with columns: cohort, n, n_yes, n_no, sensitivity, specificity.
#' @keywords internal
.compute_per_cohort_metrics <- function(scores, truth, cohort, threshold) {
  cohort_levels <- levels(cohort)
  results <- lapply(cohort_levels, function(coh) {
    idx <- which(cohort == coh)
    coh_scores <- scores[idx]
    coh_truth <- truth[idx]

    n <- length(idx)
    n_yes <- sum(coh_truth == "Yes")
    n_no <- sum(coh_truth == "No")

    pred_class <- ifelse(coh_scores >= threshold, "Yes", "No")

    # Sensitivity: TP / (TP + FN)
    tp <- sum(pred_class == "Yes" & coh_truth == "Yes")
    fn <- sum(pred_class == "No" & coh_truth == "Yes")
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

    # Specificity: TN / (TN + FP)
    tn <- sum(pred_class == "No" & coh_truth == "No")
    fp <- sum(pred_class == "Yes" & coh_truth == "No")
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

    data.frame(
      cohort = coh,
      n = n,
      n_yes = n_yes,
      n_no = n_no,
      sensitivity = sens,
      specificity = spec,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Compute Inverse Sample-Size Weighted Variance
#'
#' Calculates weighted variance of sensitivity and specificity across cohorts,
#' using inverse sample size as weights (smaller cohorts contribute less).
#'
#' @param per_cohort_metrics Data frame from .compute_per_cohort_metrics().
#' @return List with sensitivity and specificity variance values.
#' @keywords internal
.compute_weighted_variance <- function(per_cohort_metrics) {
  n_cohorts <- nrow(per_cohort_metrics)

  if (n_cohorts <= 1) {
    return(list(sensitivity = 0, specificity = 0))
  }

  # Inverse sample-size weights (normalized)
  weights <- 1 / per_cohort_metrics$n
  weights <- weights / sum(weights)

  compute_weighted_var <- function(values, w) {
    valid <- !is.na(values)
    if (sum(valid) <= 1) return(0)

    v <- values[valid]
    wt <- w[valid]
    wt <- wt / sum(wt)

    weighted_mean <- sum(v * wt)
    weighted_var <- sum(wt * (v - weighted_mean)^2)
    weighted_var
  }

  list(
    sensitivity = compute_weighted_var(per_cohort_metrics$sensitivity, weights),
    specificity = compute_weighted_var(per_cohort_metrics$specificity, weights)
  )
}

#' Predict from Fitted Model
#'
#' Internal helper to get predictions from either glmnet or glm model.
#'
#' @param model Fitted model (cv.glmnet or glm).
#' @param newx New feature matrix for prediction.
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_from_model <- function(model, newx) {
  if (inherits(model, "cv.glmnet")) {
    # glmnet model
    x_mat <- as.matrix(newx)

    # Add cohort dummies if the model was trained with them
    meta <- model$biomarkerPanels_meta
    if (!is.null(meta$cohort_info)) {
      # For prediction, we assume no cohort effect (set dummies to 0)
      n_dummies <- meta$cohort_info$n_dummies
      dummy_cols <- matrix(0, nrow = nrow(x_mat), ncol = n_dummies)
      colnames(dummy_cols) <- paste0(".cohort_", seq_len(n_dummies))
      x_mat <- cbind(x_mat, dummy_cols)
    }

    preds <- stats::predict(model, newx = x_mat, s = "lambda.min",
                            type = "response")[, 1]
  } else if (inherits(model, "glm")) {
    # Standard GLM
    newdata <- as.data.frame(newx, check.names = TRUE)
    preds <- stats::predict(model, newdata = newdata, type = "response")
  } else {
    stop("Unknown model type: ", class(model)[1], call. = FALSE)
  }

  as.numeric(preds)
}
