#' Simulated Gene Expression Data
#'
#' A synthetic log-scale gene expression dataset containing four cohorts,
#' designed for benchmarking multi-objective biomarker panel optimization.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{x_list}{A list of 4 numeric matrices, each with 100 rows (samples) and 250 columns (features/genes). Dimnames are provided.}
#'   \item{y_list}{A list of 4 factors representing binary outcome labels ("No", "Yes") for each cohort.}
#'   \item{metadata}{A list of simulation parameters including `seed`, `n_samples`, `n_features`, `informative_genes`, `signal_coefficients`, `signal_shift`, and `global_shift`.}
#' }
#'
#' @source Simulated via \code{simulate_expression_data(p = 250, n = 100, k = 4, seed = 20240220)}.
"fake_gene_expression"
