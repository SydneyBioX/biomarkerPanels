# Helper to create a simple panel with fitted model for testing
.create_test_panel_with_model <- function(features, x_train, y_train) {
  # Fit a simple glm model
  df <- data.frame(
    .response = as.integer(y_train) - 1L,
    as.data.frame(x_train[, features, drop = FALSE], check.names = TRUE)
  )
  model <- stats::glm(.response ~ ., data = df, family = stats::binomial())

  metrics <- c(sensitivity = 0.5, specificity = 0.5, auc = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(features), nrow(objective_df)))

  new(
    "BiomarkerPanelResult",
    base_features = features,
    features = features,
    metrics = metrics,
    objectives = objective_df,
    control = list(positive_class = "Yes"),
    training_data = list(),
    model = model
  )
}

# Simple scoring function for testing (used when no model needed)
.simple_scoring_fn <- function(x_selected, selected_features, truth, cohort = NULL, ...) {
  # Use row means as simple score
  rowMeans(x_selected)
}

# Helper to create an OptimizationResult for stability and frequency testing
.make_mock_optimization_result <- function() {
  solutions <- data.frame(
    solution_id = 1:3,
    sensitivity = c(0.8, 0.75, 0.82),
    specificity = c(0.7, 0.72, 0.68),
    stringsAsFactors = FALSE
  )
  solutions$features <- I(list(
    c("A--B", "A--C"),
    c("B--C", "A--C", "B--D"),
    c("A--B", "C--D")
  ))
  solutions$base_features <- I(list(
    c("A", "B", "C"),
    c("A", "B", "C", "D"),
    c("A", "B", "C", "D")
  ))

  new(
    "OptimizationResult",
    solutions = solutions,
    feature_pool = c("A", "B", "C", "D"),
    control = list(feature_transform = "pairwise_ratios"),
    training_signature = list(),
    aggregated_x = matrix(0, nrow = 1, ncol = 1),
    aggregated_y = factor("No"),
    aggregated_cohort = factor("c1")
  )
}

# Helper to conditionally skip slow tests
skip_slow_tests <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    return(invisible())
  }
  if (isTRUE(as.logical(Sys.getenv("CI")))) {
    return(invisible())
  }
  testthat::skip("Skipping slow tests: set NOT_CRAN=true to run")
}

# Helper to conditionally skip extended tests
skip_extended_tests <- function() {
  if (identical(Sys.getenv("EXTENDED_TESTS"), "true")) {
    return(invisible())
  }
  testthat::skip("Skipping extended tests: set EXTENDED_TESTS=true to run")
}

# Helper to get path to test fixture data
fixture_path <- function(filename) {
  testthat::test_path("..", "data", filename)
}
