make_opt_result <- function(points, directions, history = list()) {
  obj_names <- names(directions)
  sols <- data.frame(solution_id = seq_len(nrow(points)))
  sols$base_features <- I(replicate(nrow(points), character(0), simplify = FALSE))
  sols$features     <- I(replicate(nrow(points), character(0), simplify = FALSE))
  for (j in seq_along(obj_names)) sols[[obj_names[j]]] <- points[, j]

  new(
    "OptimizationResult",
    solutions = sols,
    feature_pool = character(0),
    control = list(objective_directions = directions),
    training_signature = list(),
    aggregated_x = matrix(numeric(0), 0, 0),
    aggregated_y = factor(character(0), levels = c("No", "Yes")),
    aggregated_cohort = factor(character(0)),
    history = history
  )
}

test_that("hypervolume matches a hand-computed 2D value", {
  skip_if_not_installed("emoa")

  pts <- rbind(c(0.9, 0.7), c(0.8, 0.85), c(0.6, 0.9))
  res <- make_opt_result(pts, c(sensitivity = "maximize", specificity = "maximize"))

  hv <- hypervolume(res, ref = c(0.5, 0.5))
  expect_equal(as.numeric(hv), 0.13, tolerance = 1e-9)
  expect_equal(attr(hv, "ref"), c(sensitivity = 0.5, specificity = 0.5))
})

test_that("hypervolume handles mixed maximize/minimize directions", {
  skip_if_not_installed("emoa")

  # AUC (maximize) vs num_features (minimize). Manually computed reference.
  pts <- rbind(c(0.85, 3), c(0.90, 5), c(0.92, 7))
  res <- make_opt_result(pts, c(auc = "maximize", num_features = "minimize"))

  hv_user <- hypervolume(res, ref = c(0.80, 8))
  # In min form: pts = (-0.85,3), (-0.90,5), (-0.92,7); ref=(-0.80,8)
  # Sort by x: (-0.92,7), (-0.90,5), (-0.85,3); y best descending: 7,5,3
  # contrib1 = (ref_x - x1) * (ref_y - y1) = (-0.80-(-0.92))*(8-7) = 0.12*1 = 0.12
  # contrib2 = (ref_x - x2) * (y1 - y2) = (-0.80-(-0.90))*(7-5) = 0.10*2 = 0.20
  # contrib3 = (ref_x - x3) * (y2 - y3) = (-0.80-(-0.85))*(5-3) = 0.05*2 = 0.10
  # total = 0.42
  expect_equal(as.numeric(hv_user), 0.42, tolerance = 1e-9)
})

test_that("hypervolume strictly increases with a new non-dominated point", {
  skip_if_not_installed("emoa")

  base <- rbind(c(0.9, 0.7), c(0.8, 0.85))
  augmented <- rbind(base, c(0.6, 0.9))
  dirs <- c(sensitivity = "maximize", specificity = "maximize")

  hv_base <- hypervolume(make_opt_result(base, dirs), ref = c(0.5, 0.5))
  hv_aug  <- hypervolume(make_opt_result(augmented, dirs), ref = c(0.5, 0.5))

  expect_gt(as.numeric(hv_aug), as.numeric(hv_base))
})

test_that("hypervolume auto reference yields a finite, positive value", {
  skip_if_not_installed("emoa")

  pts <- rbind(c(0.9, 0.7), c(0.8, 0.85), c(0.6, 0.9))
  res <- make_opt_result(pts, c(sensitivity = "maximize", specificity = "maximize"))

  hv <- hypervolume(res, ref_buffer = 0.1)
  expect_true(is.finite(hv))
  expect_gt(as.numeric(hv), 0)
  ref <- attr(hv, "ref")
  expect_lt(ref[["sensitivity"]], min(pts[, 1]))
  expect_lt(ref[["specificity"]], min(pts[, 2]))
})

test_that("hypervolume errors on bad reference length", {
  skip_if_not_installed("emoa")
  pts <- rbind(c(0.9, 0.7), c(0.8, 0.85))
  res <- make_opt_result(pts, c(sensitivity = "maximize", specificity = "maximize"))
  expect_error(hypervolume(res, ref = c(0.5)), "length 2")
})

test_that("hypervolume_history walks generations and is non-decreasing for typical fronts", {
  skip_if_not_installed("emoa")

  hist <- data.frame(
    generation = c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L),
    individual = c(1L, 2L, 1L, 2L, 3L, 1L, 2L, 3L),
    rank       = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    is_pareto  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    sensitivity = c(0.7, 0.6, 0.8, 0.7, 0.6, 0.9, 0.8, 0.6),
    specificity = c(0.6, 0.7, 0.7, 0.8, 0.85, 0.7, 0.85, 0.9),
    n_features = c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L)
  )
  hist$base_features <- I(replicate(nrow(hist), character(0), simplify = FALSE))

  res <- make_opt_result(
    hist[hist$generation == 3L, c("sensitivity", "specificity")],
    c(sensitivity = "maximize", specificity = "maximize"),
    history = hist
  )

  traj <- hypervolume_history(res, ref = c(0.5, 0.5))
  expect_equal(traj$generation, c(1L, 2L, 3L))
  expect_true(all(diff(traj$hypervolume) >= -1e-9))
  expect_equal(attr(traj, "ref"), c(sensitivity = 0.5, specificity = 0.5))
})

test_that("hypervolume_history errors when no history was recorded", {
  pts <- rbind(c(0.9, 0.7))
  res <- make_opt_result(pts, c(sensitivity = "maximize", specificity = "maximize"))
  expect_error(hypervolume_history(res), "record_history = TRUE")
})
