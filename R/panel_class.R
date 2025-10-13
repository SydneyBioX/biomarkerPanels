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
#' @export
setClass(
  "BiomarkerPanelResult",
  slots = c(
    features = "character",
    metrics = "numeric",
    objectives = "data.frame",
    control = "list",
    training_data = "list"
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
