#' NSGA Configuration Helpers
#'
#' Internal functions for configuring NSGA-II and NSGA-III parameters,
#' generating initial populations, and normalizing constraints.
#'
#' @name nsga_config
#' @noRd
NULL

#' Normalize Constraints to Standard Format
#'
#' Processes constraint specifications into a standardized list format.
#'
#' @param constraints List of constraint functions or descriptors.
#' @return List of standardized constraint objects with fun and label.
#' @noRd
.normalize_constraints <- function(constraints) {
  if (!length(constraints)) {
    return(list())
  }
  if (!is.list(constraints)) {
    stop("`constraints` must be supplied as a list.", call. = FALSE)
  }
  lapply(seq_along(constraints), function(i) {
    entry <- constraints[[i]]
    label <- names(constraints)[i]
    if (is.function(entry)) {
      fun <- entry
    } else if (is.list(entry) && !is.null(entry$fun) && is.function(entry$fun)) {
      fun <- entry$fun
      if (!is.null(entry$label) && nzchar(entry$label)) {
        label <- entry$label
      }
    } else {
      stop(
        "`constraints` entries must be functions or lists containing a `fun` element.",
        call. = FALSE
      )
    }
    if (is.null(label) || !nzchar(label)) {
      label <- sprintf("constraint_%02d", i)
    }
    list(
      fun = fun,
      label = label
    )
  })
}

#' Get Adaptive NSGA Defaults Based on Feature Pool Size and Algorithm
#'
#' Computes appropriate NSGA-II or NSGA-III hyperparameters based on the
#' dimensionality of the optimization problem and the selected algorithm.
#' Larger feature pools require more generations and larger populations to
#' adequately explore the search space. NSGA-III typically benefits from
#' larger populations for good reference point coverage.
#'
#' @param n_features Number of features in the decision space.
#' @param algorithm Algorithm to use: `"NSGA-III"` or `"NSGA-II"`.
#' @return Named list of NSGA parameters.
#' @noRd
.get_adaptive_nsga_defaults <- function(n_features, algorithm = "NSGA-III") {
  # Base parameters that don't change with problem size
  # Note: rmoo uses SBX crossover (nc=20) and polynomial mutation (nm=0.2) by default
  base_params <- list(
    pcrossover = 0.7,
    pmutation = 0.2
  )

  # NSGA-III typically needs larger populations for good reference point coverage
  if (algorithm == "NSGA-III") {
    if (n_features <= 30) {
      adaptive_params <- list(popSize = 92, maxiter = 80)
    } else if (n_features <= 100) {
      adaptive_params <- list(popSize = 156, maxiter = 180)
    } else if (n_features <= 250) {
      adaptive_params <- list(popSize = 256, maxiter = 350)
    } else {
      warning("Feature pool size is very large. NSGA-III may take a long time to run. Consider feature selection to reduce feature pool size.")
      adaptive_params <- list(popSize = 504, maxiter = 700)
    }
  } else {
    # NSGA-II defaults (existing)
    if (n_features <= 30) {
      adaptive_params <- list(popSize = 64, maxiter = 60)
    } else if (n_features <= 100) {
      adaptive_params <- list(popSize = 128, maxiter = 150)
    } else if (n_features <= 250) {
      adaptive_params <- list(popSize = 192, maxiter = 300)
    } else {
      warning("Feature pool size is very large. NSGA-II may take a long time to run. Consider feature selection to reduce feature pool size.")
      adaptive_params <- list(popSize = 256, maxiter = 500)
    }
  }

  c(adaptive_params, base_params)
}

#' Compute NSGA-III Partition Count Based on Number of Objectives
#'
#' Determines the appropriate number of partitions for generating reference
#' points in NSGA-III based on the number of objectives. More objectives
#' require fewer partitions to avoid exponential growth in reference points.
#'
#' @param n_objectives Number of objectives in the optimization problem.
#' @return Integer number of partitions.
#' @noRd
.compute_nsga3_partitions <- function(n_objectives) {
  # Rule of thumb: partitions scale inversely with objectives
  # For 2-3 objectives: 12 partitions
  # For 4-5 objectives: 6 partitions
  # For 6+ objectives: 4 partitions
  if (n_objectives <= 3) {
    12L
  } else if (n_objectives <= 5) {
    6L
  } else {
    4L
  }
}

#' Generate Sparse Initial Population Suggestions
#'
#' Creates a diverse set of initial weight vectors that span different panel
#' sizes, from `min_features` to `max_features`. This helps NSGA explore
#' the full range of panel sizes rather than converging only toward max_features.
#'
#' @param n_features Number of features in the decision space.
#' @param n_suggestions Number of suggestion vectors to generate.
#' @param min_features Minimum number of features to include (active weights).
#' @param max_features Maximum number of features to include.
#' @param seed Optional integer seed for reproducibility.
#' @return Matrix with `n_suggestions` rows and `n_features` columns, where
#'   each row is a weight vector with varying numbers of high/low values.
#' @noRd
.generate_sparse_suggestions <- function(n_features,
                                         n_suggestions = 20L,
                                         min_features = 2L,
                                         max_features = 10L,
                                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  target_sizes <- unique(round(seq(min_features, max_features,
    length.out = n_suggestions
  )))
  n_suggestions <- length(target_sizes)
  suggestions <- matrix(0, nrow = n_suggestions, ncol = n_features)

  for (i in seq_len(n_suggestions)) {
    k <- target_sizes[i]
    on_idx <- sample(n_features, size = min(k, n_features))
    suggestions[i, on_idx] <- runif(length(on_idx), 0.7, 0.95)
    off_idx <- setdiff(seq_len(n_features), on_idx)
    suggestions[i, off_idx] <- runif(length(off_idx), 0.05, 0.3)
  }
  suggestions
}

#' Minimum Base Features Required by Fitting Mode and Transform
#'
#' With pairwise transforms, n base features yield n*(n-1)/2 ratios:
#' 2 base features give a single ratio (not enough for glmnet), while
#' 3 give 3 ratios. Regularized fitting on untransformed features needs
#' at least 2 columns; unregularized fitting works from 1.
#'
#' @param regularized Logical; whether glmnet fitting is used.
#' @param feature_transform Feature transform name.
#' @return Integer minimum number of base features.
#' @noRd
.min_features_required <- function(regularized, feature_transform) {
  if (regularized &&
      feature_transform %in% c("pairwise_ratios", "pairwise_log_ratios")) {
    3L
  } else if (regularized) {
    2L
  } else {
    1L
  }
}

#' Assemble NSGA Parameters and Run the Search
#'
#' Shared back half of every `fitness_mode` of [optimize_panel()]: builds the
#' rmoo parameter list, adds NSGA-III reference points when needed,
#' optionally seeds the RNG, seeds the initial population with sparse
#' suggestions spanning panel sizes, and dispatches to [rmoo::nsga2()] or
#' [rmoo::nsga3()].
#'
#' @param algorithm `"NSGA-II"` or `"NSGA-III"`.
#' @param fitness Fitness function passed to rmoo.
#' @param n_objectives Number of objectives.
#' @param decision_dim Length of the decision vector.
#' @param nsga_args Resolved NSGA control parameters (popSize, maxiter, ...).
#' @param monitor Monitor function or `FALSE`.
#' @param min_features,max_features Panel size bounds for the sparse
#'   suggestions.
#' @param seed Optional seed forwarded to the suggestion generator.
#' @param set_seed Logical; when `TRUE`, validate `seed` and call
#'   [set.seed()] before generating suggestions (the transferable pipeline
#'   seeds earlier, at partitioning).
#' @return The raw rmoo result object.
#' @noRd
.run_nsga <- function(algorithm, fitness, n_objectives, decision_dim,
                      nsga_args, monitor, min_features, max_features,
                      seed = NULL, set_seed = FALSE) {
  nsga_params <- c(
    list(
      type = "real-valued",
      fitness = fitness,
      nObj = n_objectives,
      lower = rep(0, decision_dim),
      upper = rep(1, decision_dim),
      monitor = monitor,
      summary = FALSE
    ),
    nsga_args
  )

  # Add reference points for NSGA-III if not user-specified
  if (algorithm == "NSGA-III" && is.null(nsga_args$n_partitions)) {
    nsga_params$n_partitions <- .compute_nsga3_partitions(n_objectives)
  }

  if (set_seed && !is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("`seed` must be a single integer value.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  # Seed the initial population with solutions spanning min to max features
  n_suggestions <- min(20L, nsga_params$popSize %/% 4L)
  if (n_suggestions >= 2L) {
    nsga_params$suggestions <- .generate_sparse_suggestions(
      n_features = decision_dim,
      n_suggestions = n_suggestions,
      min_features = min_features,
      max_features = max_features,
      seed = seed
    )
  }

  if (algorithm == "NSGA-II") {
    do.call(rmoo::nsga2, nsga_params)
  } else {
    do.call(rmoo::nsga3, nsga_params)
  }
}
