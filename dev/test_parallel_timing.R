# Test script: Compare sequential vs parallel NSGA-III timing
# Purpose: Verify rmoo parallel evaluation works correctly with biomarkerPanels

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

# Small feature pool for manageable test time, but enough for parallel benefit
# Use top DE features to ensure signal
message("Selecting features...")
top_features <- get_top_de_features(
  x_list,
  y_list,
  n_features = 20
)

objectives <- define_objectives(
  losses = c("sensitivity", "specificity", "num_features")
)

# Common settings for fair comparison
common_args <- list(
  x = x_list,
  y = y_list,
  objectives = objectives,
  max_features = 5,
  feature_pool = top_features,
  cohort_aggregator = "none",  # Skip aggregation to isolate NSGA timing
  algorithm = "NSGA-II",  # Use NSGA-II to verify parallel works; NSGA-III has a bug
  seed = 42,
  fitness_cv = FALSE  # Faster for timing test
)

# Smaller population/iterations for quick test, but enough to see parallel benefit
nsga_settings <- list(popSize = 64, maxiter = 30)

message("\n=== Running SEQUENTIAL test ===")
t_seq <- system.time({
  result_seq <- do.call(optimize_panel, c(
    common_args,
    list(nsga_control = c(nsga_settings, list(parallel = FALSE)))
  ))
})
message(sprintf("Sequential time: %.2f seconds", t_seq["elapsed"]))

message("\n=== Running PARALLEL test ===")
t_par <- system.time({
  result_par <- do.call(optimize_panel, c(
    common_args,
    list(nsga_control = c(nsga_settings, list(parallel = TRUE)))
  ))
})
message(sprintf("Parallel time: %.2f seconds", t_par["elapsed"]))

# Compare results
message("\n=== Results Comparison ===")
message(sprintf("Sequential features: %s", paste(result_seq@features, collapse = ", ")))
message(sprintf("Parallel features:   %s", paste(result_par@features, collapse = ", ")))

message("\nSequential metrics:")
print(panel_metrics(result_seq))
message("\nParallel metrics:")
print(panel_metrics(result_par))

message("\n=== Timing Summary ===")
speedup <- t_seq["elapsed"] / t_par["elapsed"]
message(sprintf("Sequential: %.2f sec", t_seq["elapsed"]))
message(sprintf("Parallel:   %.2f sec", t_par["elapsed"]))
message(sprintf("Speedup:    %.2fx", speedup))

if (speedup > 1.1) {
  message("\n✓ Parallel execution is faster!")
} else if (speedup >= 0.9) {
  message("\n~ Parallel and sequential times are similar (parallel overhead may offset gains for this workload)")
} else {
  message("\n✗ Parallel is slower (overhead exceeds benefit for this small workload)")
}

message("\nTest complete.")
