# Test script: Compare sequential vs parallel NSGA timing
# Purpose: Benchmark whether parallelization helps after optimize_panel refactoring
#
# Background: Previous testing found parallel execution was ~7x slower than sequential.
# The refactoring separated panel optimization from model fitting (commit 848afa5),
# but fitness evaluation still performs glmnet fitting for each candidate.
#
# rmoo handles parallel setup automatically when parallel = TRUE:
# - Auto-detects cores via parallel::detectCores()
# - Creates cluster and registers with doParallel
# - No manual registerDoParallel() needed

library(biomarkerPanels)

# Load test fixture
fixture_path <- "tests/data/fake_gene_expression.Rds"
if (!file.exists(fixture_path)) {
  stop("Test fixture not found. Run from package root directory.")
}
test_data <- readRDS(fixture_path)

# Use 3 cohorts for more realistic testing
x_list <- test_data$x_list[1:3]
y_list <- test_data$y_list[1:3]

# Select features for testing
message("Selecting features...")
top_features <- get_top_de_features(
  x_list,
  y_list,
  n_features = 20
)

objectives <- define_objectives(
  losses = c("sensitivity", "specificity", "num_features")
)

#' Run a single benchmark scenario
#' @param scenario_name Character name for the scenario
#' @param x_list List of expression matrices
#' @param y_list List of response vectors
#' @param objectives Objective definition
#' @param feature_pool Character vector of candidate features
#' @param popSize Population size for NSGA
#' @param maxiter Maximum iterations
#' @param fitness_cv Whether to use CV for fitness evaluation
#' @param regularized Whether to use regularized fitting
#' @return List with timing results
run_benchmark <- function(scenario_name, x_list, y_list, objectives,
                          feature_pool, popSize, maxiter,
                          fitness_cv = FALSE, regularized = FALSE) {
  message(sprintf("\n=== Scenario: %s ===", scenario_name))
  message(sprintf("  popSize=%d, maxiter=%d, fitness_cv=%s, regularized=%s",
                  popSize, maxiter, fitness_cv, regularized))

  common_args <- list(
    x = x_list,
    y = y_list,
    objectives = objectives,
    max_features = 5,
    feature_pool = feature_pool,
    cohort_aggregator = "none",  # Skip aggregation to isolate NSGA timing
    algorithm = "NSGA-II",
    seed = 42,
    fitness_cv = fitness_cv,
    regularized = regularized
  )

  nsga_settings <- list(popSize = popSize, maxiter = maxiter)

  # Sequential run
  message("  Running sequential...")
  t_seq <- system.time({
    result_seq <- do.call(optimize_panel, c(
      common_args,
      list(nsga_control = c(nsga_settings, list(parallel = FALSE)))
    ))
  })

  # Parallel run
  message("  Running parallel...")
  t_par <- system.time({
    result_par <- do.call(optimize_panel, c(
      common_args,
      list(nsga_control = c(nsga_settings, list(parallel = TRUE)))
    ))
  })

  # Extract results using new API
  n_sol_seq <- n_solutions(result_seq)
  n_sol_par <- n_solutions(result_par)

  # Get features from first solution for comparison
  features_seq <- get_solution_features(result_seq, 1)
  features_par <- get_solution_features(result_par, 1)

  speedup <- t_seq["elapsed"] / t_par["elapsed"]

  message(sprintf("  Sequential: %.2f sec (%d solutions)",
                  t_seq["elapsed"], n_sol_seq))
  message(sprintf("  Parallel:   %.2f sec (%d solutions)",
                  t_par["elapsed"], n_sol_par))
  message(sprintf("  Speedup:    %.2fx", speedup))

  list(
    scenario = scenario_name,
    popSize = popSize,
    maxiter = maxiter,
    fitness_cv = fitness_cv,
    regularized = regularized,
    seq_time = t_seq["elapsed"],
    par_time = t_par["elapsed"],
    speedup = speedup,
    n_sol_seq = n_sol_seq,
    n_sol_par = n_sol_par,
    features_seq = features_seq,
    features_par = features_par
  )
}

# Define benchmark scenarios
# Test multiple configurations to find when (if ever) parallelization helps:
# - Baseline: Fast fitness = overhead likely dominates
# - CV enabled: 5x more fits per candidate
# - Large pop: More parallelizable work
# - CV + Large pop: Best case for parallel benefit

scenarios <- list(
  list(name = "Baseline (fast fitness)", popSize = 64, maxiter = 30,
       fitness_cv = FALSE, regularized = FALSE),
  list(name = "CV enabled", popSize = 64, maxiter = 30,
       fitness_cv = TRUE, regularized = TRUE),
  list(name = "Large population", popSize = 128, maxiter = 30,
       fitness_cv = FALSE, regularized = FALSE),
  list(name = "CV + Large pop", popSize = 128, maxiter = 30,
       fitness_cv = TRUE, regularized = TRUE)
)

# Run all benchmarks
results <- list()
for (scenario in scenarios) {
  result <- run_benchmark(
    scenario_name = scenario$name,
    x_list = x_list,
    y_list = y_list,
    objectives = objectives,
    feature_pool = top_features,
    popSize = scenario$popSize,
    maxiter = scenario$maxiter,
    fitness_cv = scenario$fitness_cv,
    regularized = scenario$regularized
  )
  results[[length(results) + 1]] <- result
}

# Create summary table
message("\n")
message("=" |> rep(70) |> paste(collapse = ""))
message("BENCHMARK SUMMARY")
message("=" |> rep(70) |> paste(collapse = ""))

summary_df <- data.frame(
  Scenario = sapply(results, `[[`, "scenario"),
  popSize = sapply(results, `[[`, "popSize"),
  maxiter = sapply(results, `[[`, "maxiter"),
  fitness_cv = sapply(results, `[[`, "fitness_cv"),
  Sequential = sprintf("%.1fs", sapply(results, `[[`, "seq_time")),
  Parallel = sprintf("%.1fs", sapply(results, `[[`, "par_time")),
  Speedup = sprintf("%.2fx", sapply(results, `[[`, "speedup")),
  stringsAsFactors = FALSE
)

print(summary_df, row.names = FALSE)

# Determine overall recommendation
speedups <- sapply(results, `[[`, "speedup")
any_faster <- any(speedups > 1.1)
all_slower <- all(speedups < 0.9)

message("\n")
message("-" |> rep(70) |> paste(collapse = ""))
message("RECOMMENDATION")
message("-" |> rep(70) |> paste(collapse = ""))

if (all_slower) {
  message("Parallel execution is consistently SLOWER across all scenarios.")
  message("Keep current CLAUDE.md recommendation: Do NOT use parallel = TRUE")
} else if (any_faster) {
  faster_scenarios <- summary_df$Scenario[speedups > 1.1]
  message("Parallel execution is FASTER in some scenarios:")
  message(paste("  -", faster_scenarios, collapse = "\n"))
  message("\nConsider updating CLAUDE.md with conditional guidance.")
} else {
  message("Results are mixed - parallel and sequential are roughly equivalent.")
  message("Current recommendation to avoid parallel is reasonable but not critical.")
}

message("\nBenchmark complete.")
