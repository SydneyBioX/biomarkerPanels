#' Rotating-Validation Fitness for Transferable Panel Optimization
#'
#' Infrastructure for a fitness function that rotates the train/validation
#' split each generation. With a fixed split, NSGA evaluates thousands of
#' candidates against the same validation samples and can overfit to that
#' particular partition. Rotating the split per generation is the MOO
#' analogue of mini-batch SGD: no single split can be exploited, so the
#' selection pressure favours solutions that generalise across splits.
#'
#' Within-batch consistency is load-bearing: every candidate in a single
#' NSGA batch must be scored against the *same* split, otherwise Pareto
#' dominance is no longer meaningful. The wrapper advances the split
#' counter once per batch (i.e. once per generation) and passes the
#' resolved split down to every candidate evaluation.
#'
#' @name fitness_rotating_validation
#' @keywords internal
NULL

#' Generate K Stratified Train/Validation Splits
#'
#' Produces `n_splits` disjoint-within-split train/val partitions of a
#' pooled dataset, stratified jointly by `(cohort, class)` so every split
#' preserves class balance within each cohort. Splits are deterministic
#' given the supplied seed, which keeps sweep comparisons reproducible.
#'
#' @param y Binary response factor aligned with the pool.
#' @param cohort Cohort factor aligned with the pool.
#' @param n_splits Number of splits to generate.
#' @param val_ratio Proportion of each `(cohort, class)` stratum assigned
#'   to validation in each split.
#' @param seed Optional integer seed.
#' @return List of length `n_splits`; each element is a list with integer
#'   vectors `train` and `val` indexing rows of the pool.
#' @keywords internal
.generate_stratified_splits <- function(y, cohort, n_splits, val_ratio,
                                        seed = NULL) {
  if (!is.numeric(n_splits) || length(n_splits) != 1L || n_splits < 2L) {
    stop("`n_splits` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(val_ratio) || length(val_ratio) != 1L ||
      val_ratio <= 0 || val_ratio >= 1) {
    stop("`val_ratio` must be a single numeric value in (0, 1).",
         call. = FALSE)
  }
  if (length(y) != length(cohort)) {
    stop("`y` and `cohort` must have matching length.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  y <- factor(y)
  cohort <- factor(cohort)
  strata <- interaction(cohort, y, drop = TRUE)
  strata_levels <- levels(strata)
  n <- length(y)

  splits <- vector("list", as.integer(n_splits))
  for (k in seq_len(as.integer(n_splits))) {
    val_idx <- integer(0)
    for (lv in strata_levels) {
      stratum_idx <- which(strata == lv)
      n_stratum <- length(stratum_idx)
      if (n_stratum < 2L) {
        # Stratum too small to split; keep all rows in training
        next
      }
      n_val <- max(1L, round(n_stratum * val_ratio))
      # At least one sample must remain in training for this stratum
      n_val <- min(n_val, n_stratum - 1L)
      val_idx <- c(val_idx, sample(stratum_idx, n_val))
    }
    train_idx <- setdiff(seq_len(n), val_idx)
    splits[[k]] <- list(train = sort(train_idx), val = sort(val_idx))
  }
  splits
}

#' Create Rotating-Validation Fitness Function
#'
#' Factory returning a fitness function that rotates through precomputed
#' stratified train/val splits. Each NSGA batch advances the split counter;
#' every candidate in the batch is evaluated against the same split to
#' preserve Pareto dominance. The returned `evaluate` closure — used for
#' final Pareto-solution metrics — averages metrics across all K splits,
#' which is the honest out-of-rotation number.
#'
#' @param pool_x Pooled base-feature matrix (train + val combined).
#' @param pool_y Pooled response factor.
#' @param pool_cohort Pooled cohort factor.
#' @param splits List of `(train, val)` index pairs from
#'   [.generate_stratified_splits()].
#' @param feature_pool Character vector of base feature names.
#' @param max_features Maximum base features per candidate.
#' @param objectives Objective list from [define_objectives()].
#' @param constraints Constraint list.
#' @param regularized Whether the model is elastic-net.
#' @param alpha Elastic net mixing parameter.
#' @param feature_transform Name of the feature transform to apply.
#' @param min_features_required Minimum base features per candidate.
#' @param selection_threshold `"adaptive"` or a numeric gate in (0, 1).
#' @return List with `wrapper` (vectorised NSGA fitness) and `evaluate`
#'   (single-candidate evaluator averaged over all splits).
#' @keywords internal
.make_rotating_validation_fitness <- function(
  pool_x, pool_y, pool_cohort,
  splits,
  feature_pool, max_features, objectives, constraints,
  regularized, alpha,
  feature_transform = "none",
  min_features_required = NULL,
  selection_threshold = "adaptive"
) {
  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  constraint_specs <- .normalize_constraints(constraints)
  n_splits <- length(splits)
  if (n_splits < 2L) {
    stop("Rotating-validation fitness requires `splits` with length >= 2.",
         call. = FALSE)
  }

  if (is.null(min_features_required)) {
    min_features_required <- if (regularized) 2L else 1L
  }

  # Batch counter lives in a closure environment so each call to
  # `objective_wrapper` advances to the next split deterministically.
  batch_state <- new.env(parent = emptyenv())
  batch_state$idx <- 0L

  .select_base_indices <- function(decision_vec) {
    n_pool <- length(feature_pool)
    if (identical(selection_threshold, "adaptive")) {
      ord <- order(decision_vec, feature_pool,
                   decreasing = c(TRUE, FALSE), method = "radix")
      sorted_weights <- decision_vec[ord]
      if (n_pool > max_features) {
        top_n <- min(max_features + 5L, n_pool)
        weight_diffs <- -diff(sorted_weights[seq_len(top_n)])
        search_range <- seq(min_features_required, top_n - 1L)
        if (length(search_range) > 0L) {
          gap_idx <- which.max(weight_diffs[search_range])
          natural_cutoff <- min_features_required + gap_idx
        } else {
          natural_cutoff <- min_features_required
        }
        n_selected <- max(min_features_required,
                          min(natural_cutoff, max_features))
      } else {
        above_median <- sum(decision_vec > 0.5)
        n_selected <- max(min_features_required,
                          min(above_median, max_features))
      }
      ord[seq_len(n_selected)]
    } else {
      threshold <- as.numeric(selection_threshold)
      above <- which(decision_vec > threshold)
      if (length(above) < min_features_required) {
        ord <- order(decision_vec, feature_pool,
                     decreasing = c(TRUE, FALSE), method = "radix")
        ord[seq_len(min(min_features_required, length(ord)))]
      } else if (length(above) > max_features) {
        ord <- order(decision_vec[above], feature_pool[above],
                     decreasing = c(TRUE, FALSE), method = "radix")
        above[ord[seq_len(max_features)]]
      } else {
        above
      }
    }
  }

  .score_on_split <- function(split, x_transformed, selected_base_features) {
    train_idx <- split$train
    val_idx <- split$val
    x_train <- x_transformed[train_idx, , drop = FALSE]
    x_val <- x_transformed[val_idx, , drop = FALSE]
    y_train <- pool_y[train_idx]
    y_val <- pool_y[val_idx]
    coh_train <- droplevels(pool_cohort[train_idx])
    coh_val <- droplevels(pool_cohort[val_idx])

    if (length(unique(y_train)) < 2L) {
      return(NULL)
    }

    if (regularized) {
      fit <- .fit_final_model_regularized(
        x_train, y_train, coh_train, alpha = alpha
      )
    } else {
      fit <- .fit_final_model(x_train, y_train, coh_train)
    }
    val_scores <- .predict_from_model(fit, x_val, cohort = coh_val)

    constraint_results <- if (length(constraint_specs)) {
      setNames(vapply(seq_along(constraint_specs), function(j) {
        res <- constraint_specs[[j]]$fun(
          truth = y_val,
          scores = val_scores,
          selected = selected_base_features,
          cohort = coh_val,
          x = x_val
        )
        isTRUE(res)
      }, logical(1)), vapply(constraint_specs, `[[`, character(1), "label"))
    } else {
      logical(0)
    }
    feasible <- if (length(constraint_results)) {
      all(constraint_results)
    } else {
      TRUE
    }

    metrics <- vapply(objectives, function(obj) {
      obj$fun(
        y_val,
        val_scores,
        selected = selected_base_features,
        cohort = coh_val,
        x = x_val
      )
    }, numeric(1))

    list(
      metrics = metrics,
      constraint_results = constraint_results,
      feasible = feasible,
      scores = val_scores,
      val_idx = val_idx
    )
  }

  .prepare_candidate <- function(decision_vec) {
    selected_idx <- .select_base_indices(decision_vec)
    selected_idx <- selected_idx[
      order(decision_vec[selected_idx], decreasing = TRUE)
    ]
    selected_base_features <- feature_pool[selected_idx]
    x_base <- pool_x[, selected_base_features, drop = FALSE]
    if (feature_transform != "none" && ncol(x_base) >= 2L) {
      x_transformed <- .apply_feature_transform_single(x_base, feature_transform)
      selected_features <- colnames(x_transformed)
    } else {
      x_transformed <- x_base
      selected_features <- selected_base_features
    }
    list(
      base_features = selected_base_features,
      features = selected_features,
      x_transformed = x_transformed
    )
  }

  # Evaluation on a single split — used per generation by the NSGA wrapper.
  evaluate_on_split <- function(decision_vec, split) {
    prep <- .prepare_candidate(decision_vec)
    tryCatch({
      result <- .score_on_split(
        split, prep$x_transformed, prep$base_features
      )
      if (is.null(result)) {
        return(list(
          base_features = prep$base_features,
          features = prep$features,
          scores = rep(NA_real_, length(split$val)),
          metrics = setNames(rep(NA_real_, length(objectives)),
                             names(objectives)),
          constraint_results = logical(0),
          feasible = FALSE
        ))
      }
      list(
        base_features = prep$base_features,
        features = prep$features,
        scores = result$scores,
        metrics = result$metrics,
        constraint_results = result$constraint_results,
        feasible = result$feasible
      )
    }, error = function(e) {
      stop(
        "Rotating-validation candidate failed for base features [",
        paste(prep$base_features, collapse = ", "), "]: ",
        conditionMessage(e),
        call. = FALSE
      )
    })
  }

  # Evaluation averaged over all K splits — used for final Pareto metrics.
  evaluate_across_splits <- function(decision_vec) {
    prep <- .prepare_candidate(decision_vec)
    per_split <- lapply(splits, function(split) {
      tryCatch(
        .score_on_split(split, prep$x_transformed, prep$base_features),
        error = function(e) NULL
      )
    })
    valid <- !vapply(per_split, is.null, logical(1))
    if (!any(valid)) {
      return(list(
        base_features = prep$base_features,
        features = prep$features,
        scores = rep(NA_real_, nrow(prep$x_transformed)),
        metrics = setNames(rep(NA_real_, length(objectives)),
                           names(objectives)),
        constraint_results = logical(0),
        feasible = FALSE
      ))
    }
    per_split <- per_split[valid]
    metric_mat <- do.call(rbind, lapply(per_split, `[[`, "metrics"))
    mean_metrics <- colMeans(metric_mat, na.rm = TRUE)
    feasible_vec <- vapply(per_split, `[[`, logical(1), "feasible")

    # Per-sample scores: mean over splits where the sample was in val.
    n_pool <- nrow(prep$x_transformed)
    score_sum <- numeric(n_pool)
    score_count <- integer(n_pool)
    for (r in per_split) {
      idx <- r$val_idx
      score_sum[idx] <- score_sum[idx] + as.numeric(r$scores)
      score_count[idx] <- score_count[idx] + 1L
    }
    pooled_scores <- ifelse(score_count > 0L,
                            score_sum / pmax(score_count, 1L),
                            NA_real_)
    list(
      base_features = prep$base_features,
      features = prep$features,
      scores = pooled_scores,
      metrics = mean_metrics,
      constraint_results = logical(0),
      feasible = all(feasible_vec)
    )
  }

  .PENALTY <- 1e6

  evaluate_single_for_wrapper <- function(decision_vec, split) {
    evaluated <- evaluate_on_split(decision_vec, split)
    if (length(constraint_specs) && !evaluated$feasible) {
      return(rep(.PENALTY, length(objectives)))
    }
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      if (is.na(val) || is.infinite(val)) return(.PENALTY)
      if (dir == "maximize") return(-val)
      val
    }, val = metrics, dir = objective_directions,
    SIMPLIFY = TRUE, USE.NAMES = FALSE)
    as.numeric(converted)
  }

  objective_wrapper <- function(x, ...) {
    # Advance once per call so every candidate in this batch uses the
    # same split — critical for dominance comparisons to remain valid.
    batch_state$idx <- batch_state$idx + 1L
    split_idx <- ((batch_state$idx - 1L) %% n_splits) + 1L
    current_split <- splits[[split_idx]]

    if (is.null(dim(x))) {
      return(evaluate_single_for_wrapper(x, current_split))
    }
    t(apply(x, 1, function(row) {
      evaluate_single_for_wrapper(row, current_split)
    }))
  }

  list(
    wrapper = objective_wrapper,
    evaluate = evaluate_across_splits,
    n_splits = n_splits
  )
}
