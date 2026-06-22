# Tests for graphical summarization functions like plot_feature_stability (ggplot2 integration).
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
