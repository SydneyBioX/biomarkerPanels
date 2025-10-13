#' biomarkerPanels: Multi-objective Optimization for Compact Biomarker Panels
#'
#' The biomarkerPanels package helps researchers derive small, interpretable
#' biomarker panels that satisfy competing performance criteria such as high
#' sensitivity at fixed specificity, low assay cost, and robustness across
#' cohorts. The package integrates with Bioconductor data structures including
#' `SummarizedExperiment` to streamline analysis from expression data to
#' deployable panel definitions.
#'
#' @section Modules:
#' * Data simulation tools for benchmarking algorithms (`simulate_expression_data`).
#' * Objective construction helpers (`define_objectives`).
#' * Optimization engines (`optimize_panel`).
#' * Evaluation utilities (`evaluate_panel`, `plot_pareto_front`).
#'
#' @docType package
#' @name biomarkerPanels
NULL
