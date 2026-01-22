#' Biomarker Panel Result Class
#'
#' Stores multi-objective optimization results including selected biomarkers,
#' performance summaries, and optimization metadata.
#'
#' @slot features Character vector of selected biomarkers.
#' @slot metrics Named numeric vector summarizing sensitivity, specificity, etc.
#' @slot objectives Data frame describing objective values per candidate solution.
#' @slot control List of optimization control parameters.
#' @slot training_data Signature of the training dataset (e.g., sample and assay info).
#' @slot model Fitted model object (e.g., glm) trained on the selected features.
#' @export
setClass(
  "BiomarkerPanelResult",
  slots = c(
    features = "character",
    metrics = "numeric",
    objectives = "data.frame",
    control = "list",
    training_data = "list",
    model = "ANY"
  )
)

#' Panel Metrics Accessor
#'
#' Retrieve the metrics summary from a `BiomarkerPanelResult` object.
#'
#' @param object A `BiomarkerPanelResult`.
#' @return Named numeric vector of metrics.
#' @export
setGeneric("panel_metrics", function(object) standardGeneric("panel_metrics"))

#' @describeIn panel_metrics Return the named metrics vector.
#' @export
setMethod(
  "panel_metrics",
  signature = "BiomarkerPanelResult",
  definition = function(object) {
    object@metrics
  }
)

#' Panel Model Accessor
#'
#' Retrieve the fitted model from a `BiomarkerPanelResult` object.
#'
#' @param object A `BiomarkerPanelResult`.
#' @return The fitted model object (e.g., glm).
#' @export
setGeneric("panel_model", function(object) standardGeneric("panel_model"))

#' @describeIn panel_model Return the fitted model.
#' @export
setMethod(

  "panel_model",
  signature = "BiomarkerPanelResult",
  definition = function(object) {
    object@model
  }
)

#' Transferable Biomarker Panel Result Class
#'
#' Extends `BiomarkerPanelResult` with additional slots for train/validation/held-out
#' split strategy and Neyman-Pearson threshold selection for cross-cohort
#' generalizability.
#'
#' @slot np_threshold Numeric threshold selected via Neyman-Pearson classification.
#' @slot np_alpha Type I error rate for NP threshold selection.
#' @slot np_delta Tolerance parameter for NP threshold selection.
#' @slot per_cohort_metrics Data frame with per-cohort performance metrics.
#' @slot weighted_variance List containing sensitivity and specificity variance
#'   (inverse sample-size weighted across cohorts).
#' @slot validation_metrics List of validation set metrics.
#' @slot partition_info List containing partition sizes and ratios.
#' @export
setClass(
  "TransferablePanelResult",
  contains = "BiomarkerPanelResult",
  slots = c(
    np_threshold = "numeric",
    np_alpha = "numeric",
    np_delta = "numeric",
    per_cohort_metrics = "data.frame",
    weighted_variance = "list",
    validation_metrics = "list",
    partition_info = "list"
  )
)

#' NP Threshold Accessor
#'
#' Retrieve the Neyman-Pearson classification threshold from a
#' `TransferablePanelResult` object.
#'
#' @param object A `TransferablePanelResult`.
#' @return Numeric threshold value.
#' @export
setGeneric("np_threshold", function(object) standardGeneric("np_threshold"))

#' @describeIn np_threshold Return the NP threshold.
#' @export
setMethod(
  "np_threshold",
  signature = "TransferablePanelResult",
  definition = function(object) {
    object@np_threshold
  }
)

#' Per-Cohort Metrics Accessor
#'
#' Retrieve per-cohort performance metrics from a `TransferablePanelResult` object.
#'
#' @param object A `TransferablePanelResult`.
#' @return Data frame with cohort, n, n_yes, n_no, sensitivity, specificity.
#' @export
setGeneric("per_cohort_metrics", function(object) standardGeneric("per_cohort_metrics"))

#' @describeIn per_cohort_metrics Return per-cohort performance data frame.
#' @export
setMethod(

  "per_cohort_metrics",
  signature = "TransferablePanelResult",
  definition = function(object) {
    object@per_cohort_metrics
  }
)

#' Weighted Variance Accessor
#'
#' Retrieve inverse sample-size weighted variance of sensitivity and specificity
#' across cohorts from a `TransferablePanelResult` object.
#'
#' @param object A `TransferablePanelResult`.
#' @return List with `sensitivity` and `specificity` variance values.
#' @export
setGeneric("weighted_variance", function(object) standardGeneric("weighted_variance"))

#' @describeIn weighted_variance Return weighted variance list.
#' @export
setMethod(
  "weighted_variance",
  signature = "TransferablePanelResult",
  definition = function(object) {
    object@weighted_variance
  }
)
