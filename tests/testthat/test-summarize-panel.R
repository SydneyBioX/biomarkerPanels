# Tests for summarization functions: inclusion frequencies and frequency-based panel assembly.
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
