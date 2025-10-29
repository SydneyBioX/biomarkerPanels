#' Evaluate a Biomarker Panel
#'
#' Compute objective values and additional summary statistics for a fitted panel
#' on validation data. This function extracts the selected features from the
#' validation set and computes predictions using the default scoring function:
#' a logistic regression fit on the selected features (with a centred/standardised
#' fallback if the model fails to converge).
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
#' @param cohort_aggregator Optional transformation applied to `x` before
#'   evaluation. Defaults to the aggregator stored in the fitted panel
#'   (currently `"pairwise_ratios"`), keeping training and validation pipelines
#'   aligned. Future work: expose richer harmonisation options.
#' @param cutoff_prob Classification probability threshold used for confusion
#'   matrix summaries and highlight point on the ROC curve. Defaults to `0.5`.
#' @param positive Label treated as the positive class when computing confusion
#'   matrices; defaults to `"Yes"`.
#' @return A list with `metrics` (named numeric vector), `objectives`
#'   (data.frame with columns objective, value, direction), `confusion`
#'   (2x2 matrix of counts), and `roc` (list containing a `curve` data.frame,
#'   `highlight` point, optional `plot`, and `evalm` helper data).
#' @export
evaluate_panel <- function(panel, x, y,
                           objectives = define_objectives(
                             losses = c("sensitivity", "specificity", "c_index")
                           ),
                           cohort = NULL,
                           assay = NULL,
                           scoring_fn = NULL,
                           cohort_aggregator = NULL,
                           cutoff_prob = 0.5,
                           positive = "Yes") {
  stopifnot(inherits(panel, "BiomarkerPanelResult"))

  stored_aggregator <- panel@control$cohort_aggregator
  if (is.null(cohort_aggregator)) {
    cohort_aggregator <- if (is.null(stored_aggregator)) "none" else stored_aggregator
  } else {
    cohort_aggregator <- match.arg(cohort_aggregator, c("pairwise_ratios", "none"))
  }

  if (is.list(x)) {
    if (!is.null(cohort)) {
      warning("`cohort` argument is ignored when `x` is supplied as a list of cohorts.",
              call. = FALSE)
    }
    prepared <- .prepare_cohort_inputs(x, y, assay = assay, aggregator = cohort_aggregator)
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
      colnames(x_mat) <- sprintf("feature_%04d", seq_len(ncol(x_mat)))
    }
    x_mat <- .apply_cohort_aggregator(list(x_mat), cohort_aggregator)[[1]]
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

  sensitivity_point <- loss_sensitivity(
    truth,
    scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )
  specificity_point <- loss_specificity(
    truth,
    scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )

  confusion <- .compute_confusion_matrix(
    truth = truth,
    scores = scores,
    cutoff_prob = cutoff_prob,
    positive = positive
  )

  roc_curve <- .compute_roc_curve(
    truth = truth,
    scores = scores,
    positive = positive
  )

  highlight_point <- data.frame(
    threshold = cutoff_prob,
    tpr = sensitivity_point,
    fpr = 1 - specificity_point,
    sensitivity = sensitivity_point,
    specificity = specificity_point,
    stringsAsFactors = FALSE
  )

  roc_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    roc_plot <- ggplot2::ggplot(roc_curve, ggplot2::aes(x = fpr, y = tpr)) +
      ggplot2::geom_path(color = "#377EB8", linewidth = 0.8) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(
        data = highlight_point,
        ggplot2::aes(x = fpr, y = tpr),
        color = "#E41A1C",
        size = 2.5
      ) +
      ggplot2::labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        title = "ROC Curve"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal()
  }

  roc_auc <- loss_c_index(truth, scores, positive = positive)

  list(
    metrics = objective_values,
    objectives = data.frame(
      objective = names(objectives),
      value = objective_values,
      direction = vapply(objectives, `[[`, character(1), "direction"),
      stringsAsFactors = FALSE
    ),
    confusion = confusion,
    roc = list(
      curve = roc_curve,
      highlight = highlight_point,
      auc = roc_auc,
      plot = roc_plot,
      evalm = data.frame(
        obs = truth,
        prob = scores,
        stringsAsFactors = FALSE
      )
    ),
    cutoff = cutoff_prob,
    scores = scores
  )
}

.compute_confusion_matrix <- function(truth, scores, cutoff_prob, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }
  negative <- setdiff(levels_truth, positive)
  if (!length(negative)) {
    stop("Unable to determine negative class for confusion matrix.", call. = FALSE)
  }
  negative <- negative[1]
  predicted_positive <- scores >= cutoff_prob
  predicted <- ifelse(predicted_positive, positive, negative)
  predicted <- factor(predicted, levels = levels_truth)
  table_res <- table(
    truth = truth,
    predicted = predicted
  )
  structure(
    as.matrix(table_res),
    cutoff = cutoff_prob,
    positive = positive,
    negative = negative
  )
}

.compute_roc_curve <- function(truth, scores, positive) {
  truth <- ensure_binary_response(truth, positive = positive)
  levels_truth <- levels(truth)
  if (!positive %in% levels_truth) {
    positive <- levels_truth[length(levels_truth)]
  }
  negative <- setdiff(levels_truth, positive)
  if (!length(negative)) {
    stop("Unable to determine negative class for ROC.", call. = FALSE)
  }
  negative <- negative[1]

  thresholds <- sort(unique(scores), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  pos_total <- sum(truth == positive)
  neg_total <- sum(truth == negative)

  compute_point <- function(thresh) {
    predicted_positive <- scores >= thresh
    tp <- sum(predicted_positive & truth == positive)
    fp <- sum(predicted_positive & truth == negative)
    tpr <- if (pos_total == 0) NA_real_ else tp / pos_total
    fpr <- if (neg_total == 0) NA_real_ else fp / neg_total
    c(tpr = tpr, fpr = fpr)
  }

  mat <- vapply(thresholds, compute_point, numeric(2))
  df <- data.frame(
    threshold = thresholds,
    tpr = mat["tpr", ],
    fpr = mat["fpr", ],
    sensitivity = mat["tpr", ],
    specificity = 1 - mat["fpr", ],
    stringsAsFactors = FALSE
  )
  df[order(df$fpr, df$tpr), , drop = FALSE]
}
