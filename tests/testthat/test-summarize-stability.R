# Tests for mathematical stability analyses of selected features across solutions, such as Jaccard similarity matrices.
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
