#' Difficulty-Stratified Accuracy Metric
#'
#' Computes balanced accuracy across easy and hard sample strata, rewarding
#' panels that correctly classify clinically difficult samples.
#'
#' @name metric_difficulty
NULL

#' Easy-Hard Balanced Accuracy
#'
#' Computes accuracy within each difficulty stratum (easy/hard) and returns
#' their mean. This rewards panels that perform well on hard-to-classify
#' samples rather than relying on easy cases to inflate overall accuracy.
#'
#' The `difficulty` vector must be supplied via the `params` argument of
#' [define_objectives()] or [build_objectives()]:
#'
#' ```
#' define_objectives(
#'   metrics = c("easy_hard_accuracy", "num_features"),
#'   params = list(easy_hard_accuracy = list(difficulty = sample_labels))
#' )
#' ```
#'
#' @param truth Binary outcome; coerced with [ensure_binary_response()].
#' @param scores Numeric scores or probabilities.
#' @param selected Ignored; kept for signature compatibility.
#' @param difficulty Factor or character vector of `"easy"` / `"hard"` labels,
#'   same length as `truth`. Must be supplied.
#' @param cutoff_prob Classification probability cutoff applied to `scores`.
#'   Ignored when `cutoff_strategy` is not `"fixed"`.
#' @param cutoff_strategy Strategy for computing cutoff. One of `"fixed"`
#'   (use `cutoff_prob`), `"prevalence"` (cutoff = class prevalence), or
#'   `"youden"` (optimal Youden's J). Default is `"fixed"`.
#' @param positive Label treated as the positive ("event") class.
#' @return Balanced accuracy across easy/hard strata, between 0 and 1.
#'   Returns `NA_real_` if classification is undefined in both strata.
#'
#' @note The `difficulty` vector is captured at objective-definition time.
#'   This works with [optimize_panel()] which evaluates on the full pooled
#'   training data. For [optimize_panel_transferable()], which subsets data
#'   into train/validation partitions, the `difficulty` vector would be
#'   misaligned with the subsetted `truth`---use only with the standard
#'   pipeline.
#'
#' @export
metric_easy_hard_accuracy <- function(truth, scores = NULL, selected = NULL,
                                      difficulty = NULL,
                                      cutoff_prob = 0.5,
                                      cutoff_strategy = c("fixed", "prevalence",
                                                          "youden"),
                                      positive = "Yes") {
  truth <- ensure_binary_response(truth)
  if (is.null(scores)) {
    stop("`scores` must be supplied to compute easy-hard accuracy.",
         call. = FALSE)
  }
  if (is.null(difficulty)) {
    stop("`difficulty` must be supplied. Pass it via `params` in ",
         "`define_objectives()` or `build_objectives()`.", call. = FALSE)
  }
  difficulty <- factor(difficulty)
  if (!all(levels(difficulty) %in% c("easy", "hard"))) {
    stop("`difficulty` must contain only 'easy' and 'hard' labels.",
         call. = FALSE)
  }
  if (length(difficulty) != length(truth)) {
    stop("`difficulty` must have the same length as `truth`.", call. = FALSE)
  }

  cutoff_strategy <- match.arg(cutoff_strategy)
  cutoff <- .compute_cutoff(truth, scores, cutoff_prob, cutoff_strategy,
                            positive)
  predicted <- scores >= cutoff
  correct <- truth == positive & predicted | truth != positive & !predicted

  # Compute accuracy per stratum
  strata <- levels(difficulty)
  accs <- vapply(strata, function(s) {
    idx <- difficulty == s
    if (!any(idx)) return(NA_real_)
    mean(correct[idx])
  }, numeric(1))

  valid <- !is.na(accs)
  if (!any(valid)) return(NA_real_)
  mean(accs[valid])
}
