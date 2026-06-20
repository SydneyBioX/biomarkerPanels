test_that("select_panel_top_sensitivity chooses highest external sensitivity", {
  perf <- data.frame(
    solution_id = c("s1", "s2", "s3"),
    external_sens_at_spec90 = c(0.72, 0.81, 0.81),
    internal_sens_at_spec90 = c(0.82, 0.88, 0.83),
    panel_size = c(4L, 5L, 4L)
  )

  res <- select_panel_top_sensitivity(perf, max_loss = 0.1, panel_size_col = "panel_size")
  expect_equal(res$solution_id, "s3")
  expect_equal(res$performance_loss, res$internal_sens_at_spec90 - res$external_sens_at_spec90)
})

test_that("compute_inclusion_frequencies aggregates across panels", {
  objectives <- data.frame(
    solution_id = c(1, 1, 2, 2),
    objective = c("sens", "spec", "sens", "spec"),
    value = c(0.1, 0.2, 0.3, 0.4),
    direction = c("maximize", "maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    c("A", "B"),
    c("A", "B"),
    c("B", "C", "D"),
    c("B", "C", "D")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8, specificity = 0.75),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  freq <- compute_inclusion_frequencies(list(panel, c("C", "D", "E")))
  expect_equal(freq$feature, c("B", "C", "D", "A", "E"))
  expect_equal(freq$count, c(2L, 2L, 2L, 1L, 1L))
  expect_equal(freq$proportion, c(2/3, 2/3, 2/3, 1/3, 1/3))

  panel_features <- select_panel_inclusion_frequency(freq, panel_size = 4)
  expect_equal(panel_features, c("B", "C", "D", "A"))
})

test_that("select_panel_by_pathway integrates pathway knowledge", {
  objectives <- data.frame(
    solution_id = c(1, 1, 2, 2),
    objective = c("sens", "spec", "sens", "spec"),
    value = c(0.1, 0.2, 0.3, 0.4),
    direction = c("maximize", "maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    c("A", "B"),
    c("A", "B"),
    c("B", "C", "D"),
    c("B", "C", "D")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8, specificity = 0.75),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  freq <- compute_inclusion_frequencies(list(panel, c("C", "D", "E")))

  pathway_mapping <- data.frame(
    feature = c("A", "B", "C", "D", "E"),
    pathway = c("P1", "P1", "P2", "P3", "P3"),
    stringsAsFactors = FALSE
  )

  feature_perf <- data.frame(
    feature = c("A", "B", "C", "D", "E"),
    external_sens_at_spec90 = c(0.55, 0.72, 0.68, 0.4, 0.6),
    stringsAsFactors = FALSE
  )

  selection <- select_panel_by_pathway(
    high_sensitivity_panels = list(panel),
    inclusion_frequency = freq,
    pathway_mapping = pathway_mapping,
    relevant_pathways = c("P1", "P2"),
    feature_performance = feature_perf,
    panel_size = 4,
    min_sensitivity = 0.5
  )

  expect_equal(selection$features, c("B", "C", "A"))
  expect_true(all(selection$pathway_assignments$pathway %in% c("P1", "P2")))
})

# ============================================================================
# analyze_feature_stability tests
# ============================================================================

test_that("analyze_feature_stability returns correct structure", {
  objectives <- data.frame(
    solution_id = c(1, 1, 2, 2, 3, 3),
    objective = rep(c("sens", "spec"), 3),
    value = c(0.8, 0.7, 0.75, 0.72, 0.82, 0.68),
    direction = rep("maximize", 6),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    c("A", "B", "C"),
    c("A", "B", "C"),
    c("B", "C", "D"),
    c("B", "C", "D"),
    c("A", "C", "D", "E"),
    c("A", "C", "D", "E")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_s3_class(stability, "FeatureStabilityResult")
  expect_named(stability, c("frequencies", "jaccard_matrix", "n_solutions", "mean_jaccard"))

  expect_equal(nrow(stability$frequencies), 5L)
  expect_equal(colnames(stability$frequencies), c("feature", "count", "proportion"))
  expect_equal(stability$frequencies$feature[1], "C")
  expect_equal(stability$frequencies$count[1], 3L)

  expect_equal(dim(stability$jaccard_matrix), c(3L, 3L))
  expect_equal(rownames(stability$jaccard_matrix), c("1", "2", "3"))
  expect_equal(unname(diag(stability$jaccard_matrix)), rep(1.0, 3))
  expect_equal(stability$jaccard_matrix[1, 2], stability$jaccard_matrix[2, 1])

  expect_equal(stability$n_solutions, 3L)
  expect_true(stability$mean_jaccard >= 0 && stability$mean_jaccard <= 1)
})

test_that("analyze_feature_stability computes Jaccard correctly", {
  objectives <- data.frame(
    solution_id = c(1, 2),
    objective = c("sens", "sens"),
    value = c(0.8, 0.75),
    direction = c("maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    c("A", "B"),
    c("B", "C")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_equal(stability$jaccard_matrix[1, 2], 1/3, tolerance = 1e-10)
  expect_equal(stability$mean_jaccard, 1/3, tolerance = 1e-10)
})

test_that("analyze_feature_stability handles empty result", {
  panel <- new(
    "BiomarkerPanelResult",
    features = character(),
    metrics = numeric(),
    objectives = data.frame(),
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_s3_class(stability, "FeatureStabilityResult")
  expect_equal(nrow(stability$frequencies), 0L)
  expect_equal(dim(stability$jaccard_matrix), c(0L, 0L))
  expect_equal(stability$n_solutions, 0L)
  expect_true(is.na(stability$mean_jaccard))
})

test_that("analyze_feature_stability handles single solution", {
  objectives <- data.frame(
    solution_id = c(1, 1),
    objective = c("sens", "spec"),
    value = c(0.8, 0.7),
    direction = c("maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(c("A", "B"), c("A", "B")))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_equal(stability$n_solutions, 1L)
  expect_equal(dim(stability$jaccard_matrix), c(1L, 1L))
  expect_equal(stability$jaccard_matrix[1, 1], 1.0)
  expect_true(is.na(stability$mean_jaccard))
})

test_that("analyze_feature_stability handles empty features in solution", {
  objectives <- data.frame(
    solution_id = c(1, 2),
    objective = c("sens", "sens"),
    value = c(0.8, 0.75),
    direction = c("maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    character(),
    c("A", "B")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_equal(stability$jaccard_matrix[1, 1], 1.0)
  expect_equal(stability$jaccard_matrix[1, 2], 0.0)
  expect_equal(stability$jaccard_matrix[2, 2], 1.0)
})

# ============================================================================
# plot_feature_stability tests
# ============================================================================

test_that("plot_feature_stability returns ggplot object", {
  skip_if_not_installed("ggplot2")

  objectives <- data.frame(
    solution_id = c(1, 1, 2, 2),
    objective = rep(c("sens", "spec"), 2),
    value = c(0.8, 0.7, 0.75, 0.72),
    direction = rep("maximize", 4),
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(
    c("A", "B", "C"),
    c("A", "B", "C"),
    c("B", "C", "D"),
    c("B", "C", "D")
  ))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)
  p <- plot_feature_stability(stability)

  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_stability respects top_n parameter", {
  skip_if_not_installed("ggplot2")

  objectives <- data.frame(
    solution_id = 1,
    objective = "sens",
    value = 0.8,
    direction = "maximize",
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(LETTERS[1:15]))

  panel <- new(
    "BiomarkerPanelResult",
    features = LETTERS[1:5],
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)
  p <- plot_feature_stability(stability, top_n = 5)

  expect_s3_class(p, "ggplot")
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  expect_lte(nrow(plot_data), 5L)
})

test_that("plot_feature_stability handles empty stability result", {
  skip_if_not_installed("ggplot2")

  panel <- new(
    "BiomarkerPanelResult",
    features = character(),
    metrics = numeric(),
    objectives = data.frame(),
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)
  p <- plot_feature_stability(stability)

  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_stability validates inputs", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_feature_stability(list(a = 1)),
    "must be output from analyze_feature_stability"
  )

  objectives <- data.frame(
    solution_id = 1,
    objective = "sens",
    value = 0.8,
    direction = "maximize",
    stringsAsFactors = FALSE
  )
  objectives$features <- I(list(c("A", "B")))

  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = objectives,
    control = list(),
    training_data = list()
  )

  stability <- analyze_feature_stability(panel)

  expect_error(
    plot_feature_stability(stability, top_n = -1),
    "positive integer"
  )

  expect_error(
    plot_feature_stability(stability, highlight_threshold = 1.5),
    "between 0 and 1"
  )
})

# ============================================================================
# OptimizationResult support and feature_type tests
# ============================================================================

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

test_that("compute_inclusion_frequencies accepts OptimizationResult", {
  opt <- .make_mock_optimization_result()
  freq <- compute_inclusion_frequencies(list(opt))

  expect_true(is.data.frame(freq))
  expect_equal(colnames(freq), c("feature", "count", "proportion"))
  # "A--C" appears in solutions 1 and 2, "A--B" in 1 and 3
  expect_true("A--C" %in% freq$feature)
  expect_true("A--B" %in% freq$feature)
  expect_equal(freq$proportion[freq$feature == "A--C"], 2/3)
})

test_that("compute_inclusion_frequencies accepts unwrapped OptimizationResult", {
  opt <- .make_mock_optimization_result()
  freq <- compute_inclusion_frequencies(opt)

  expect_true(is.data.frame(freq))
  expect_true(nrow(freq) > 0L)
})

test_that("compute_inclusion_frequencies feature_type = base_features works", {
  opt <- .make_mock_optimization_result()
  freq <- compute_inclusion_frequencies(list(opt), feature_type = "base_features")

  expect_true(is.data.frame(freq))
  # Base features are A, B, C, D
  expect_true(all(freq$feature %in% c("A", "B", "C", "D")))
  # A, B, C appear in all 3 solutions
  expect_equal(freq$proportion[freq$feature == "A"], 1.0)
  # D appears in solutions 2 and 3
  expect_equal(freq$proportion[freq$feature == "D"], 2/3)
})

test_that("analyze_feature_stability accepts OptimizationResult", {
  opt <- .make_mock_optimization_result()
  stability <- analyze_feature_stability(opt)

  expect_s3_class(stability, "FeatureStabilityResult")
  expect_equal(stability$n_solutions, 3L)
  expect_equal(dim(stability$jaccard_matrix), c(3L, 3L))
  expect_true(stability$mean_jaccard >= 0 && stability$mean_jaccard <= 1)
})

test_that("analyze_feature_stability with base_features on OptimizationResult", {
  opt <- .make_mock_optimization_result()
  stability <- analyze_feature_stability(opt, feature_type = "base_features")

  expect_s3_class(stability, "FeatureStabilityResult")
  expect_true(all(stability$frequencies$feature %in% c("A", "B", "C", "D")))
  # Jaccard for base features: sol1={A,B,C}, sol2={A,B,C,D}, sol3={A,B,C,D}
  # J(1,2) = 3/4 = 0.75, J(1,3) = 3/4 = 0.75, J(2,3) = 4/4 = 1.0
  expect_equal(stability$jaccard_matrix[1, 2], 3/4, tolerance = 1e-10)
  expect_equal(stability$jaccard_matrix[2, 3], 1.0, tolerance = 1e-10)
})

test_that("analyze_feature_stability rejects base_features for BiomarkerPanelResult", {
  panel <- new(
    "BiomarkerPanelResult",
    features = c("A", "B"),
    metrics = c(sensitivity = 0.8),
    objectives = data.frame(),
    control = list(),
    training_data = list()
  )

  expect_error(
    analyze_feature_stability(panel, feature_type = "base_features"),
    "base_features.*only supported for OptimizationResult"
  )
})

test_that("analyze_feature_stability rejects invalid input", {
  expect_error(
    analyze_feature_stability(data.frame(a = 1)),
    "OptimizationResult or BiomarkerPanelResult"
  )
})
