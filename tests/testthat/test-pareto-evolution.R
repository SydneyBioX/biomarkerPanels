make_evo_result <- function(history) {
  obj_names <- c("sensitivity", "specificity")
  sols <- data.frame(solution_id = 1L)
  sols$base_features <- I(list(character(0)))
  sols$features     <- I(list(character(0)))
  for (n in obj_names) sols[[n]] <- 1

  new(
    "OptimizationResult",
    solutions = sols,
    feature_pool = character(0),
    control = list(objective_directions = c(sensitivity = "maximize",
                                            specificity = "maximize")),
    training_signature = list(),
    aggregated_x = matrix(numeric(0), 0, 0),
    aggregated_y = factor(character(0), levels = c("No", "Yes")),
    aggregated_cohort = factor(character(0)),
    history = history
  )
}

make_history <- function() {
  rows <- list(
    list(gen = 1L, ind = 1L, rank = 1L, sens = 0.60, spec = 0.65, feats = c("A", "B")),
    list(gen = 1L, ind = 2L, rank = 1L, sens = 0.70, spec = 0.55, feats = c("A", "C")),
    list(gen = 1L, ind = 3L, rank = 2L, sens = 0.40, spec = 0.50, feats = c("D")),
    list(gen = 2L, ind = 1L, rank = 1L, sens = 0.75, spec = 0.70, feats = c("A", "B")),
    list(gen = 2L, ind = 2L, rank = 1L, sens = 0.80, spec = 0.65, feats = c("A", "B", "C")),
    list(gen = 2L, ind = 3L, rank = 2L, sens = 0.50, spec = 0.50, feats = c("E")),
    list(gen = 3L, ind = 1L, rank = 1L, sens = 0.85, spec = 0.80, feats = c("A", "C")),
    list(gen = 3L, ind = 2L, rank = 1L, sens = 0.90, spec = 0.75, feats = c("B", "C")),
    list(gen = 3L, ind = 3L, rank = 1L, sens = 0.78, spec = 0.85, feats = c("A", "B", "C"))
  )

  df <- data.frame(
    generation = vapply(rows, `[[`, integer(1), "gen"),
    individual = vapply(rows, `[[`, integer(1), "ind"),
    rank = vapply(rows, `[[`, integer(1), "rank"),
    is_pareto = vapply(rows, function(r) r$rank == 1L, logical(1)),
    sensitivity = vapply(rows, `[[`, numeric(1), "sens"),
    specificity = vapply(rows, `[[`, numeric(1), "spec"),
    n_features = vapply(rows, function(r) length(r$feats), integer(1))
  )
  df$base_features <- I(lapply(rows, `[[`, "feats"))
  df
}

test_that("plot_pareto_evolution returns a magick animation with one frame per generation", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(make_history())

  anim <- plot_pareto_evolution(res, width = 320, height = 200, fps = 2)
  expect_s3_class(anim, "magick-image")
  expect_equal(length(anim), 3L)
})

test_that("plot_pareto_evolution subsamples generations when given a count", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(make_history())

  anim <- plot_pareto_evolution(res, width = 320, height = 200, fps = 2,
                                generations = 2L)
  expect_equal(length(anim), 2L)
})

test_that("plot_pareto_evolution accepts an explicit generation vector", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(make_history())

  anim <- plot_pareto_evolution(res, width = 320, height = 200, fps = 2,
                                generations = c(1L, 3L))
  expect_equal(length(anim), 2L)
})

test_that("plot_pareto_evolution writes a GIF to disk when path is supplied", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(make_history())
  out <- tempfile(fileext = ".gif")
  on.exit(unlink(out), add = TRUE)

  plot_pareto_evolution(res, width = 320, height = 200, fps = 2, path = out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

test_that("plot_pareto_evolution errors when no history is recorded", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(list())
  expect_error(plot_pareto_evolution(res), "record_history = TRUE")
})

test_that("plot_pareto_evolution errors on unknown metric column", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("magick")

  res <- make_evo_result(make_history())
  expect_error(
    plot_pareto_evolution(res, x_metric = "not_a_metric"),
    "missing from nsga_history"
  )
})
