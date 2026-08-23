#' Fitness Caching Utilities
#'
#' Cache primitives keyed by selected base-feature panel, the shared panel
#' transformer, and population-level fitness evaluation with panel-level
#' deduplication (NSGA populations contain many duplicate panels; without the
#' dedup every duplicate would refit its model). Panel selection lives in
#' `R/panel_selection.R`; the fast GLM path lives in `R/glm_fitting.R`.
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

  # Check for reference_feature attribute on any matrix
  ref_col <- NULL
  for (mat in matrices) {
    rc <- attr(mat, "reference_feature")
    if (!is.null(rc)) {
      ref_col <- rc
      break
    }
  }

  function(selected_base_features) {
    if (!is.null(ref_col) && !(ref_col %in% selected_base_features)) {
      selected_base_features <- c(ref_col, selected_base_features)
    }

    key <- .panel_key(selected_base_features)
    cached <- .cache_get(cache, key)
    if (!is.null(cached)) {
      return(cached)
    }

    base_matrices <- lapply(matrices, function(mat) {
      sub_mat <- mat[, selected_base_features, drop = FALSE]
      if (!is.null(ref_col)) {
        attr(sub_mat, "reference_feature") <- ref_col
      }
      sub_mat
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

