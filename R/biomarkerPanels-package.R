#' biomarkerPanels: Multi-objective Optimization for Compact Biomarker Panels
#'
#' The biomarkerPanels package helps researchers derive small, interpretable
#' biomarker panels that satisfy competing performance criteria such as high
#' sensitivity at fixed specificity, low assay cost, and transferability across
#' cohorts. The package integrates with common Bioconductor data structures
#' including `SummarizedExperiment` to streamline analysis from expression data
#' to deployable panel definitions.
#'
#' @section Modules:
#' * Objective construction helpers (`define_objectives`).
#' * Optimization engines (`optimize_panel`).
#' * Evaluation utilities (`evaluate_panel`, `plot_pareto_front`).
#'
#' @docType package
#' @name biomarkerPanels
NULL
