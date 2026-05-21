#' Fitness Selection, Caching, and Fast GLM Helpers
#'
#' Internal utilities shared by standard and transferable panel optimization.
#'
#' @name fitness_cache_utils
#' @keywords internal
NULL

.validate_cache_controls <- function(cache_fitness, cache_max_entries) {
  if (!is.logical(cache_fitness) || length(cache_fitness) != 1L ||
      is.na(cache_fitness)) {
    stop("`cache_fitness` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(cache_max_entries) || length(cache_max_entries) != 1L ||
      is.na(cache_max_entries) || cache_max_entries < 0) {
    stop("`cache_max_entries` must be a non-negative number or Inf.",
         call. = FALSE)
  }
  invisible(NULL)
}

.panel_key <- function(base_features, context = NULL) {
  if (!length(base_features)) {
    key <- "<empty>"
  } else {
    features <- enc2utf8(as.character(base_features))
    key <- paste0(nchar(features, type = "bytes"), ":", features,
                  collapse = "|")
  }
  if (!is.null(context)) {
    paste0(context, "::", key)
  } else {
    key
  }
}

.new_fitness_cache <- function(max_entries = Inf) {
  .validate_cache_controls(TRUE, max_entries)
  cache <- new.env(parent = emptyenv())
  cache$data <- new.env(parent = emptyenv())
  cache$keys <- character()
  cache$max_entries <- as.numeric(max_entries)
  cache
}

.cache_get <- function(cache, key) {
  if (is.null(cache)) {
    return(NULL)
  }
  if (!exists(key, envir = cache$data, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = cache$data, inherits = FALSE)
}

.cache_set <- function(cache, key, value) {
  if (is.null(cache) || cache$max_entries <= 0) {
    return(invisible(value))
  }
  if (!exists(key, envir = cache$data, inherits = FALSE)) {
    if (is.finite(cache$max_entries)) {
      while (length(cache$keys) >= cache$max_entries && length(cache$keys) > 0L) {
        old_key <- cache$keys[[1L]]
        rm(list = old_key, envir = cache$data)
        cache$keys <- cache$keys[-1L]
      }
    }
    cache$keys <- c(cache$keys, key)
  }
  assign(key, value, envir = cache$data)
  invisible(value)
}

.select_panel_indices <- function(decision_vec, feature_pool, max_features,
                                  min_features_required,
                                  selection_threshold = "adaptive") {
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
      n_selected <- min(natural_cutoff, max_features)
      n_selected <- max(n_selected, min_features_required)
    } else {
      above_median <- sum(decision_vec > 0.5)
      n_selected <- max(min_features_required, min(above_median, max_features))
    }
    selected_idx <- ord[seq_len(n_selected)]
  } else {
    threshold <- as.numeric(selection_threshold)
    above_threshold <- which(decision_vec > threshold)

    if (length(above_threshold) < min_features_required) {
      ord <- order(decision_vec, feature_pool,
                   decreasing = c(TRUE, FALSE), method = "radix")
      selected_idx <- ord[seq_len(min(min_features_required, length(ord)))]
    } else if (length(above_threshold) > max_features) {
      weights_above <- decision_vec[above_threshold]
      names_above <- feature_pool[above_threshold]
      ord <- order(weights_above, names_above,
                   decreasing = c(TRUE, FALSE), method = "radix")
      selected_idx <- above_threshold[ord[seq_len(max_features)]]
    } else {
      selected_idx <- above_threshold
    }
  }

  selected_idx[order(decision_vec[selected_idx], decreasing = TRUE)]
}

.make_panel_selector <- function(feature_pool, max_features,
                                 min_features_required,
                                 selection_threshold = "adaptive") {
  force(feature_pool)
  force(max_features)
  force(min_features_required)
  force(selection_threshold)
  function(decision_vec) {
    selected_idx <- .select_panel_indices(
      decision_vec = decision_vec,
      feature_pool = feature_pool,
      max_features = max_features,
      min_features_required = min_features_required,
      selection_threshold = selection_threshold
    )
    base_features <- feature_pool[selected_idx]
    list(
      indices = selected_idx,
      base_features = base_features,
      key = .panel_key(base_features)
    )
  }
}

.make_panel_transformer <- function(matrices, feature_transform,
                                    cache_max_entries = Inf) {
  if (is.matrix(matrices) || is.data.frame(matrices)) {
    matrices <- list(x = matrices)
  }
  matrices <- lapply(matrices, function(mat) {
    if (is.matrix(mat)) mat else as.matrix(mat)
  })
  cache <- .new_fitness_cache(cache_max_entries)
  force(feature_transform)

  function(selected_base_features) {
    key <- .panel_key(selected_base_features)
    cached <- .cache_get(cache, key)
    if (!is.null(cached)) {
      return(cached)
    }

    base_matrices <- lapply(matrices, function(mat) {
      mat[, selected_base_features, drop = FALSE]
    })

    if (feature_transform != "none" && length(selected_base_features) >= 2L) {
      transformed <- lapply(base_matrices, .apply_feature_transform_single,
                            transform_name = feature_transform)
      selected_features <- colnames(transformed[[1L]])
    } else {
      transformed <- base_matrices
      selected_features <- selected_base_features
    }

    out <- list(
      matrices = transformed,
      features = selected_features
    )
    .cache_set(cache, key, out)
    out
  }
}

.convert_metrics_to_objectives <- function(metrics, objective_directions,
                                           penalty = 1e6) {
  converted <- mapply(function(val, dir) {
    if (is.na(val) || is.infinite(val)) {
      return(penalty)
    }
    if (dir == "maximize") {
      return(-val)
    }
    val
  }, val = metrics, dir = objective_directions,
  SIMPLIFY = TRUE, USE.NAMES = FALSE)
  as.numeric(converted)
}

.evaluate_fitness_population <- function(x, selector, evaluate_selection,
                                         n_objectives, cache = NULL,
                                         cache_fitness = TRUE,
                                         context = NULL) {
  eval_one <- function(decision_vec) {
    selection <- selector(decision_vec)
    key <- .panel_key(selection$base_features, context = context)
    if (isTRUE(cache_fitness)) {
      cached <- .cache_get(cache, key)
      if (!is.null(cached)) {
        return(cached)
      }
    }
    value <- evaluate_selection(selection)
    if (isTRUE(cache_fitness)) {
      .cache_set(cache, key, value)
    }
    value
  }

  if (is.null(dim(x))) {
    return(eval_one(x))
  }

  if (!isTRUE(cache_fitness)) {
    return(t(apply(x, 1L, eval_one)))
  }

  selections <- lapply(seq_len(nrow(x)), function(i) selector(x[i, ]))
  keys <- vapply(selections, function(selection) {
    .panel_key(selection$base_features, context = context)
  }, character(1))
  result <- matrix(NA_real_, nrow = nrow(x), ncol = n_objectives)

  for (key in unique(keys)) {
    idx <- which(keys == key)
    cached <- .cache_get(cache, key)
    if (!is.null(cached)) {
      value <- cached
    } else {
      value <- evaluate_selection(selections[[idx[[1L]]]])
      .cache_set(cache, key, value)
    }
    result[idx, ] <- matrix(rep(value, length(idx)),
                            nrow = length(idx), byrow = TRUE)
  }

  result
}

.cohort_dummy_matrix <- function(cohort, levels) {
  n <- length(cohort)
  if (length(levels) <= 1L) {
    return(NULL)
  }
  cohort_factor <- factor(cohort, levels = levels)
  dummies <- matrix(0, nrow = n, ncol = length(levels) - 1L)
  colnames(dummies) <- paste0(".cohort_", seq_len(ncol(dummies)))
  for (j in seq_along(levels[-1L])) {
    dummies[, j] <- ifelse(
      is.na(cohort_factor),
      NA_real_,
      as.numeric(cohort_factor == levels[-1L][[j]])
    )
  }
  dummies
}

.prepare_glm_design_terms <- function(cohort_train = NULL, cohort_new = NULL,
                                      predict_cohort = c("observed", "reference")) {
  predict_cohort <- match.arg(predict_cohort)
  if (is.null(cohort_train) || length(unique(cohort_train)) <= 1L) {
    return(list(train_extra = NULL, new_extra = NULL))
  }

  cohort_train <- factor(cohort_train)
  levels_train <- levels(cohort_train)
  train_extra <- .cohort_dummy_matrix(cohort_train, levels_train)

  if (is.null(cohort_new)) {
    cohort_new <- cohort_train
  }

  if (predict_cohort == "reference") {
    new_extra <- matrix(0, nrow = length(cohort_new),
                        ncol = length(levels_train) - 1L)
    colnames(new_extra) <- colnames(train_extra)
  } else {
    new_extra <- .cohort_dummy_matrix(cohort_new, levels_train)
  }

  list(
    train_extra = train_extra,
    new_extra = new_extra,
    levels = levels_train,
    predict_cohort = predict_cohort
  )
}

.prepare_cv_glm_design_terms <- function(cohort, fold_ids) {
  if (is.null(cohort)) {
    return(NULL)
  }
  unique_folds <- sort(unique(fold_ids))
  terms <- lapply(unique_folds, function(fold) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    .prepare_glm_design_terms(
      cohort_train = cohort[train_idx],
      cohort_new = cohort[test_idx],
      predict_cohort = "observed"
    )
  })
  names(terms) <- as.character(unique_folds)
  terms
}

.build_glm_design_matrix <- function(x, extra = NULL) {
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  design <- cbind("(Intercept)" = 1, x_mat)
  if (!is.null(extra) && ncol(extra) > 0L) {
    design <- cbind(design, extra)
  }
  design
}

.fit_predict_binomial_glm <- function(x_train, truth, x_new = NULL,
                                      cohort_train = NULL, cohort_new = NULL,
                                      predict_cohort = c("observed", "reference"),
                                      design_terms = NULL) {
  predict_cohort <- match.arg(predict_cohort)
  if (is.null(x_train) || ncol(x_train) == 0L) {
    stop("No features selected for GLM scoring.", call. = FALSE)
  }
  if (length(unique(truth)) < 2L) {
    stop("Response contains only one class. Cannot fit classification model.",
         call. = FALSE)
  }
  if (is.null(x_new)) {
    x_new <- x_train
  }
  if (is.null(design_terms)) {
    design_terms <- .prepare_glm_design_terms(
      cohort_train = cohort_train,
      cohort_new = if (is.null(cohort_new)) cohort_train else cohort_new,
      predict_cohort = predict_cohort
    )
  }

  x_design <- .build_glm_design_matrix(x_train, design_terms$train_extra)
  new_design <- .build_glm_design_matrix(x_new, design_terms$new_extra)
  y_vec <- as.integer(truth) - 1L

  fit <- stats::glm.fit(
    x = x_design,
    y = y_vec,
    family = stats::binomial()
  )
  coefs <- fit$coefficients
  coefs[is.na(coefs)] <- 0
  preds <- stats::binomial()$linkinv(drop(new_design %*% coefs))
  as.numeric(preds)
}
