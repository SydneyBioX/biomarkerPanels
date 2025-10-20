#' Evaluate a Biomarker Panel
#'
#' Compute objective values and additional summary statistics for a fitted panel
#' on validation data. This function extracts the selected features from the
#' validation set and computes predictions using the default scoring function
#' (row means of selected features).
#' When evaluating multiple cohorts, feature alignment is performed via the
#' simple intersection of shared column names; future releases will offer more
#' flexible strategies.
#'
#' @param panel A `BiomarkerPanelResult`.
#' @param x Validation feature matrix, `data.frame`, `SummarizedExperiment`, or a
#'   list of such objects for multi-cohort evaluation.
#' @param y Validation binary response factor (or list of factors when `x` is a
#'   list).
#' @param objectives Optional override for objectives. Defaults to sensitivity,
#'   specificity, and c_index.
#' @param cohort Optional factor identifying cohort membership for each sample
#'   when `x` represents a single cohort. Ignored when `x` is provided as a list.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @param scoring_fn Optional custom scoring function. If `NULL`, uses the
#'   default row-mean aggregation. Must have signature
#'   `function(x_selected, selected_features, truth, ...)`.
#' @return A list with `metrics` (named numeric vector) and `objectives`
#'   (data.frame with columns objective, value, direction).
#' @export
evaluate_panel <- function(panel, x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "c_index")
                           ),
                           cohort = NULL,
                           assay = NULL,
                           scoring_fn = NULL) {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  if (is.list(x)) {
    if (!is.null(cohort)) {
      warning("`cohort` argument is ignored when `x` is supplied as a list of cohorts.",
              call. = FALSE)
    }
    prepared <- .prepare_cohort_inputs(x, y, assay = assay)
    x_mat <- prepared$x
    truth <- prepared$truth
    cohort_vec <- prepared$cohort
  } else {
    x_mat <- .extract_feature_matrix(x, assay = assay)
    truth <- ensure_binary_response(y)
    if (nrow(x_mat) != length(truth)) {
      stop("`x` and `y` must have matching sample sizes.", call. = FALSE)
    }
    if (is.null(colnames(x_mat))) {
      stop("`x` must have column names in order to align with panel features.",
           call. = FALSE)
    }
    if (is.null(cohort)) {
      cohort_vec <- factor(rep("cohort_01", nrow(x_mat)), levels = "cohort_01")
    } else {
      if (length(cohort) != nrow(x_mat)) {
        stop("Length of `cohort` must match the number of samples in `x`.", call. = FALSE)
      }
      cohort_vec <- factor(cohort)
    }
  }

  selected <- panel@features
  if (length(selected) == 0L) {
    stop("Panel has no selected features.", call. = FALSE)
  }

  if (!all(selected %in% colnames(x_mat))) {
    missing <- setdiff(selected, colnames(x_mat))
    stop("Selected feature(s) not found in validation data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  x_selected <- x_mat[, selected, drop = FALSE]

  if (is.null(scoring_fn)) {
    scoring_fn <- .default_scoring_fn
  }

  if (!is.function(scoring_fn)) {
    stop("`scoring_fn` must be a function.", call. = FALSE)
  }

  score_args <- list(
    x_selected = x_selected,
    selected_features = selected,
    truth = truth,
    cohort = cohort_vec
  )
  scores <- do.call(scoring_fn, score_args)

  if (!is.numeric(scores) || length(scores) != length(truth)) {
    stop("`scoring_fn` must return a numeric vector matching the number of samples.",
         call. = FALSE)
  }

  objective_values <- vapply(objectives, function(obj) {
    obj$fun(
      truth,
      scores,
      selected = selected,
      cohort = cohort_vec,
      x = x_selected
    )
  }, numeric(1))

  list(
    metrics = objective_values,
    objectives = data.frame(
      objective = names(objectives),
      value = objective_values,
      direction = vapply(objectives, `[[`, character(1), "direction"),
      stringsAsFactors = FALSE
    )
  )
}
