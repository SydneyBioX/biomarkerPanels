#' Optimize Biomarker Panels with NSGA-II
#'
#' Wrapper around [mco::nsga2()] that composes registered loss functions into a
#' multi-objective search for compact biomarker panels. Candidate solutions are
#' represented as weights over a feature pool; the top `max_features` features
#' define a panel whose performance is evaluated via the selected losses.
#' Inputs may be a single cohort (`matrix`, `data.frame`, or
#' `SummarizedExperiment`) or multiple cohorts supplied as lists of such objects.
#'
#' @param x Matrix-like object, `SummarizedExperiment`, or list of matrices /
#'   experiments representing one or more cohorts.
#' @param y Binary response (`factor`, `character`, or `logical`) aligned with
#'   `x`. When `x` is a list, `y` must be a list of the same length.
#' @param objectives Named list of objective descriptors as returned by
#'   [define_objectives()].
#' @param max_features Maximum number of biomarkers permitted in a panel.
#' @param feature_pool Optional subset of feature identifiers (names or integer
#'   indices) considered during optimization. Defaults to all features.
#' @param scoring_fn Function producing per-sample scores from the selected
#'   features. Signature:
#'   `function(x_selected, selected_features, truth, cohort = NULL, ...)`.
#' @param nsga_control Named list of arguments passed to [mco::nsga2()]. Defaults
#'   to `list(popsize = 64, generations = 60, cprob = 0.7, cdist = 5,
#'   mprob = 0.2, mdist = 10)`.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @return A `BiomarkerPanelResult` with Pareto-optimal solutions summarised in
#'   the `objectives` slot.
#' @export
optimize_panel <- function(x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "num_features")
                           ),
                           max_features = 5L,
                           feature_pool = NULL,
                           scoring_fn = NULL,
                           nsga_control = list(),
                           assay = NULL) {
  if (!requireNamespace("mco", quietly = TRUE)) {
    stop("The 'mco' package is required. Install it via BiocManager::install('mco').",
         call. = FALSE)
  }

  inputs <- .prepare_cohort_inputs(x, y, assay = assay)
  x_mat <- inputs$x
  truth <- inputs$truth
  cohort <- inputs$cohort

  feature_names <- colnames(x_mat)
  if (is.null(feature_names)) {
    feature_names <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    colnames(x_mat) <- feature_names
  }

  if (is.null(feature_pool)) {
    feature_pool <- feature_names
  } else {
    feature_pool <- .resolve_feature_pool(feature_pool, feature_names)
  }

  if (!length(feature_pool)) {
    stop("`feature_pool` produced zero features.", call. = FALSE)
  }

  if (max_features < 1L) {
    stop("`max_features` must be at least 1.", call. = FALSE)
  }

  if (max_features > length(feature_pool)) {
    max_features <- length(feature_pool)
  }

  x_pool <- x_mat[, feature_pool, drop = FALSE]
  decision_dim <- ncol(x_pool)

  if (decision_dim > 200) {
    warning("Optimizing over more than 200 features may be slow; consider ",
            "reducing `feature_pool` for exploration.")
  }

  if (is.null(scoring_fn)) {
    scoring_fn <- .default_scoring_fn
  }

  if (!is.function(scoring_fn)) {
    stop("`scoring_fn` must be a function.", call. = FALSE)
  }

  objective_directions <- vapply(objectives, `[[`, character(1), "direction")
  names(objective_directions) <- names(objectives)

  nsga_defaults <- list(
    popsize = 64,
    generations = 60,
    cprob = 0.7,
    cdist = 5,
    mprob = 0.2,
    mdist = 10
  )
  nsga_args <- utils::modifyList(nsga_defaults, nsga_control)

  evaluate_candidate <- function(decision_vec) {
    ord <- order(decision_vec, decreasing = TRUE)
    selected_idx <- ord[seq_len(min(max_features, length(ord)))]
    selected_features <- colnames(x_pool)[selected_idx]
    x_selected <- x_pool[, selected_features, drop = FALSE]

    score_args <- list(
      x_selected = x_selected,
      selected_features = selected_features,
      truth = truth
    )
    if (!is.null(cohort)) {
      score_args$cohort <- cohort
    }
    scores <- do.call(scoring_fn, score_args)
    if (!is.numeric(scores) || length(scores) != nrow(x_pool)) {
      stop("`scoring_fn` must return a numeric vector matching the number of samples.",
           call. = FALSE)
    }
    metrics <- vapply(objectives, function(obj) {
      obj$fun(truth, scores, selected = selected_features)
    }, numeric(1))
    list(
      features = selected_features,
      scores = scores,
      metrics = metrics
    )
  }

  objective_wrapper <- function(decision_vec) {
    evaluated <- evaluate_candidate(decision_vec)
    metrics <- evaluated$metrics
    converted <- mapply(function(val, dir) {
      if (is.na(val)) {
        return(Inf)
      }
      if (dir == "maximize") {
        return(-val)
      }
      val
    }, val = metrics, dir = objective_directions, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    as.numeric(converted)
  }

  nsga_params <- c(
    list(
      fn = objective_wrapper,
      idim = decision_dim,
      odim = length(objectives),
      lower.bounds = rep(0, decision_dim),
      upper.bounds = rep(1, decision_dim)
    ),
    nsga_args
  )

  nsga_result <- do.call(mco::nsga2, nsga_params)

  if (is.null(dim(nsga_result$par))) {
    nsga_result$par <- matrix(nsga_result$par, nrow = 1)
  }

  solutions <- lapply(seq_len(nrow(nsga_result$par)), function(i) {
    decision_vec <- nsga_result$par[i, ]
    evaluate_candidate(decision_vec)
  })

  metric_matrix <- do.call(rbind, lapply(solutions, `[[`, "metrics"))
  colnames(metric_matrix) <- names(objectives)

  select_primary <- function(idx_vec, direction) {
    if (direction == "maximize") {
      which.max(idx_vec)
    } else {
      which.min(idx_vec)
    }
  }

  primary_obj <- names(objectives)[1]
  primary_dir <- objective_directions[[primary_obj]]
  primary_idx <- select_primary(metric_matrix[, primary_obj], primary_dir)

  primary_solution <- solutions[[primary_idx]]
  primary_metrics <- primary_solution$metrics

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

  panel <- new(
    "BiomarkerPanelResult",
    features = primary_solution$features,
    metrics = setNames(as.numeric(primary_metrics), names(objectives)),
    objectives = objective_df,
    control = list(
      max_features = max_features,
      feature_pool = feature_pool,
      nsga2 = nsga_args,
      scoring_function = deparse(substitute(scoring_fn))
    ),
    training_data = list(
      n = nrow(x_mat),
      p = ncol(x_mat),
      class_balance = table(truth),
      feature_pool_size = length(feature_pool),
      num_cohorts = length(inputs$cohort_names),
      cohort_labels = inputs$cohort_names,
      cohort_counts = inputs$cohort_counts
    )
  )

  panel
}

.default_scoring_fn <- function(x_selected, selected_features, truth,
                                cohort = NULL, ...) {
  if (is.null(x_selected) || ncol(x_selected) == 0L) {
    return(rep(0, length(truth)))
  }
  if (ncol(x_selected) == 1L) {
    return(as.numeric(x_selected[, 1]))
  }
  rowMeans(x_selected)
}

.prepare_cohort_inputs <- function(x, y, assay = NULL) {
  if (is.list(x)) {
    if (!is.list(y)) {
      stop("When `x` is a list, `y` must also be a list.", call. = FALSE)
    }
    if (length(x) != length(y)) {
      stop("`x` and `y` lists must have the same length.", call. = FALSE)
    }
    cohort_names <- names(x)
    if (is.null(cohort_names) || any(cohort_names == "")) {
      cohort_names <- sprintf("cohort_%02d", seq_along(x))
    }

    matrices <- lapply(seq_along(x), function(i) {
      .extract_feature_matrix(x[[i]], assay = assay)
    })

    base_features <- colnames(matrices[[1]])
    if (is.null(base_features)) {
      base_features <- sprintf("feature_%04d", seq_len(ncol(matrices[[1]])))
      colnames(matrices[[1]]) <- base_features
    }

    matrices <- lapply(seq_along(matrices), function(i) {
      mat <- matrices[[i]]
      if (is.null(colnames(mat))) {
        colnames(mat) <- base_features
      }
      missing <- setdiff(base_features, colnames(mat))
      if (length(missing)) {
        stop("Feature(s) missing in cohort ", shQuote(cohort_names[[i]]), ": ",
             paste(missing, collapse = ", "), call. = FALSE)
      }
      mat[, base_features, drop = FALSE]
    })

    counts <- vapply(matrices, nrow, integer(1))
    response_list <- lapply(seq_along(y), function(i) {
      yi <- ensure_binary_response(y[[i]])
      if (length(yi) != counts[[i]]) {
        stop("Length of `y[[", i, "]]` must match number of samples in `x[[", i, "]]`.",
             call. = FALSE)
      }
      yi
    })

    combined_x <- do.call(rbind, matrices)
    truth <- factor(
      unlist(lapply(response_list, as.character), use.names = FALSE),
      levels = levels(response_list[[1]])
    )
    cohort <- factor(rep(cohort_names, counts), levels = cohort_names)

    list(
      x = combined_x,
      truth = truth,
      cohort = cohort,
      cohort_names = cohort_names,
      cohort_counts = setNames(as.list(counts), cohort_names)
    )
  } else {
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    cohort_names <- "cohort_01"
    list(
      x = x_mat,
      truth = truth,
      cohort = factor(rep(cohort_names, nrow(x_mat)), levels = cohort_names),
      cohort_names = cohort_names,
      cohort_counts = setNames(list(nrow(x_mat)), cohort_names)
    )
  }
}

.extract_feature_matrix <- function(x, assay = NULL) {
  if (inherits(x, "SummarizedExperiment")) {
    assays <- SummarizedExperiment::assayNames(x)
    if (is.null(assay)) {
      assay <- assays[1]
    }
    return(SummarizedExperiment::assay(x, assay))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`x` must be a matrix-like object.", call. = FALSE)
  }
  mode(x) <- "numeric"
  x
}

.resolve_feature_pool <- function(pool, feature_names) {
  if (is.numeric(pool)) {
    pool <- feature_names[pool]
  }
  if (!all(pool %in% feature_names)) {
    missing <- setdiff(pool, feature_names)
    stop("Feature(s) not found in `x`: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  unique(pool)
}
