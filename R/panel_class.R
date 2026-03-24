#' Optimization Result Class
#'
#' Stores multi-objective optimization results from NSGA-II/III including the
#' Pareto-optimal solutions, feature pool, and training data. This class is
#' returned by [optimize_panel()] and serves as input to [fit_panel()] for
#' model fitting.
#'
#' @slot solutions Data frame with one row per Pareto solution. Contains columns:
#'   `solution_id` (integer), `base_features` (list-column of character vectors
#'   with original feature names), `features` (list-column with transformed
#'   feature names used for model fitting), and one column per objective with
#'   numeric values.
#' @slot feature_pool Character vector of base features considered during optimization.
#' @slot control List of optimization control parameters including `feature_transform`.
#' @slot training_signature List with metadata about the training dataset.
#' @slot aggregated_x Raw (untransformed) training feature matrix containing base
#'   features. The feature transform is applied on-the-fly during [fit_panel()].
#' @slot aggregated_y Training response vector (factor).
#' @slot aggregated_cohort Cohort membership factor for training samples.
#' @export
setClass(
  "OptimizationResult",
  slots = c(
    solutions = "data.frame",
    feature_pool = "character",
    control = "list",
    training_signature = "list",
    aggregated_x = "ANY",
    aggregated_y = "ANY",
    aggregated_cohort = "ANY"
  )
)

#' Solutions Accessor
#'
#' Retrieve the Pareto-optimal solutions from an `OptimizationResult` object.
#'
#' @param object An `OptimizationResult`.
#' @return Data frame with one row per solution containing features and metrics.
#' @export
setGeneric("solutions", function(object) standardGeneric("solutions"))

#' @describeIn solutions Return the solutions data frame.
#' @export
setMethod(
  "solutions",
  signature = "OptimizationResult",
  definition = function(object) {
    object@solutions
  }
)

#' Number of Solutions
#'
#' Get the number of Pareto-optimal solutions in an `OptimizationResult`.
#'
#' @param object An `OptimizationResult`.
#' @return Integer count of solutions.
#' @export
setGeneric("n_solutions", function(object) standardGeneric("n_solutions"))

#' @describeIn n_solutions Return the number of solutions.
#' @export
setMethod(
  "n_solutions",
  signature = "OptimizationResult",
  definition = function(object) {
    nrow(object@solutions)
  }
)

#' Get Solution Features
#'
#' Extract the feature set for a specific Pareto solution.
#'
#' @param object An `OptimizationResult`.
#' @param solution_id Integer ID of the solution (1-based).
#' @return Character vector of feature names in the solution.
#' @export
setGeneric("get_solution_features", function(object, solution_id)
  standardGeneric("get_solution_features"))

#' @describeIn get_solution_features Return features for the specified solution.
#' @export
setMethod(
  "get_solution_features",
  signature = c(object = "OptimizationResult", solution_id = "numeric"),
  definition = function(object, solution_id) {
    solution_id <- as.integer(solution_id)
    if (solution_id < 1L || solution_id > nrow(object@solutions)) {
      stop("solution_id must be between 1 and ", nrow(object@solutions), call. = FALSE)
    }
    object@solutions$features[[solution_id]]
  }
)

#' Get Solution Base Features
#'
#' Extract the base (untransformed) feature set for a specific Pareto solution.
#' These are the original feature names before any transformation (e.g., pairwise
#' ratios) is applied.
#'
#' @param object An `OptimizationResult`.
#' @param solution_id Integer ID of the solution (1-based).
#' @return Character vector of base feature names in the solution.
#' @export
setGeneric("get_solution_base_features", function(object, solution_id)
  standardGeneric("get_solution_base_features"))

#' @describeIn get_solution_base_features Return base features for the specified solution.
#' @export
setMethod(
  "get_solution_base_features",
  signature = c(object = "OptimizationResult", solution_id = "numeric"),
  definition = function(object, solution_id) {
    solution_id <- as.integer(solution_id)
    if (solution_id < 1L || solution_id > nrow(object@solutions)) {
      stop("solution_id must be between 1 and ", nrow(object@solutions), call. = FALSE)
    }
    object@solutions$base_features[[solution_id]]
  }
)

#' Biomarker Panel Result Class
#'
#' Stores a fitted biomarker panel including selected biomarkers,
#' performance summaries, and the trained model. This class is
#' returned by [fit_panel()] and serves as input to [evaluate_panel()].
#' @name BiomarkerPanelResult
#' @slot base_features Character vector of original (untransformed) biomarker names.
#'   These are the features selected by the optimization algorithm before any
#'   transformation is applied.
#' @slot features Character vector of transformed biomarker names used for model
#'   fitting. For pairwise transforms, these are the ratio names (e.g., "A--B").
#' @slot metrics Named numeric vector summarizing sensitivity, specificity, etc.
#' @slot objectives Data frame describing objective values per candidate solution.
#' @slot control List of optimization control parameters including `feature_transform`.
#' @slot training_data Signature of the training dataset (e.g., sample and assay info).
#' @slot model Fitted model object (e.g., glm or cv.glmnet) trained on the
#'   transformed features.
#' @export
setClass(
  "BiomarkerPanelResult",
  slots = c(
    base_features = "character",
    features = "character",
    metrics = "numeric",
    objectives = "data.frame",
    control = "list",
    training_data = "list",
    model = "ANY"
  )
)

#' Panel Base Features Accessor
#'
#' Retrieve the original (untransformed) feature names from a `BiomarkerPanelResult`
#' object. These are the features selected by the optimization before any
#' transformation was applied.
#'
#' @param object A `BiomarkerPanelResult`.
#' @return Character vector of base feature names.
#' @export
setGeneric("panel_base_features", function(object) standardGeneric("panel_base_features"))

#' @describeIn panel_base_features Return the base feature names.
#' @export
setMethod(
  "panel_base_features",
  signature = "BiomarkerPanelResult",
  definition = function(object) {
    object@base_features
  }
)

#' Panel Features Accessor
#'
#' Retrieve the transformed feature names from a `BiomarkerPanelResult` object.
#' For pairwise transforms, these are ratio names like "GeneA--GeneB".
#'
#' @param object A `BiomarkerPanelResult`.
#' @return Character vector of transformed feature names.
#' @export
setGeneric("panel_features", function(object) standardGeneric("panel_features"))

#' @describeIn panel_features Return the transformed feature names.
#' @export
setMethod(
  "panel_features",
  signature = "BiomarkerPanelResult",
  definition = function(object) {
    object@features
  }
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
