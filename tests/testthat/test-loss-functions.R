test_that("loss_sensitivity and loss_specificity behave as expected", {
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.8, 0.2)
  expect_equal(loss_sensitivity(truth, scores), 1)
  expect_equal(loss_specificity(truth, scores), 1)
  expect_false(is.na(loss_balanced_accuracy(truth, scores)))
})

test_that("loss_c_index computes concordance", {
  truth <- factor(rep(c("No", "Yes"), each = 10), levels = c("No", "Yes"))
  scores <- c(seq(0.1, 0.9, length.out = 10), seq(0.2, 1, length.out = 10))
  expect_gt(loss_c_index(truth, scores), 0.5)
})

test_that("loss_num_features counts correctly", {
  expect_equal(loss_num_features(selected = c("g1", "g2")), 2)
  expect_equal(loss_num_features(selected = c(TRUE, FALSE, TRUE)), 2)
  expect_equal(loss_num_features(selected = NULL), 0)
})

test_that("build_objectives wraps registered losses", {
  objs <- build_objectives(c("sensitivity", "num_features"),
                           params = list(sensitivity = list(cutoff_prob = 0.7)))
  truth <- factor(c("No", "Yes", "Yes", "No"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.6, 0.3)
  selected <- c("a", "b", "c")
  expect_equal(names(objs), c("sensitivity", "num_features"))
  expect_equal(objs$sensitivity$direction, "maximize")
  expect_equal(objs$num_features$direction, "minimize")
  expect_true(is.numeric(objs$sensitivity$fun(truth, scores, selected)))
  expect_equal(objs$num_features$fun(truth, scores, selected), 3)
})

test_that("register_loss_function adds new entries", {
  custom <- function(truth, scores = NULL, selected = NULL, ...) 42
  name <- paste0("custom_", sample(1000, 1))
  register_loss_function(name, custom, direction = "maximize", overwrite = TRUE)
  objs <- build_objectives(name)
  expect_equal(objs[[name]]$fun(NULL, NULL, NULL), 42)
})

test_that("cohort-aware losses compute transfer metrics", {
  truth <- factor(
    c("No", "Yes", "Yes", "No", "Yes", "No"),
    levels = c("No", "Yes")
  )
  scores <- c(0.1, 0.9, 0.8, 0.2, 0.6, 0.4)
  cohort <- factor(c("A", "A", "B", "B", "B", "A"))
  expect_equal(
    loss_min_cohort_sensitivity(truth, scores, cohort = cohort),
    min(
      loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
      loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
    )
  )
  expect_equal(
    loss_min_cohort_specificity(truth, scores, cohort = cohort),
    min(
      loss_specificity(truth[cohort == "A"], scores[cohort == "A"]),
      loss_specificity(truth[cohort == "B"], scores[cohort == "B"])
    )
  )
  gap <- loss_cohort_sensitivity_gap(truth, scores, cohort = cohort)
  expect_gte(gap, 0)
  expect_equal(gap, max(
    loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
    loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
  ) - min(
    loss_sensitivity(truth[cohort == "A"], scores[cohort == "A"]),
    loss_sensitivity(truth[cohort == "B"], scores[cohort == "B"])
  ))

  brier <- loss_max_cohort_brier(truth, scores, cohort = cohort)
  expect_gte(brier, 0)

  x <- matrix(
    c(1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 7),
    nrow = 6, ncol = 2, byrow = FALSE
  )
  shift <- loss_max_cohort_mean_shift(truth, scores, cohort = cohort, x = x)
  expect_gte(shift, 0)

  objs <- build_objectives(c("min_cohort_sensitivity", "max_cohort_mean_shift"))
  expect_equal(
    objs$min_cohort_sensitivity$fun(truth, scores, cohort = cohort),
    loss_min_cohort_sensitivity(truth, scores, cohort = cohort)
  )
  expect_equal(
    objs$max_cohort_mean_shift$fun(truth, scores, cohort = cohort, x = x),
    loss_max_cohort_mean_shift(truth, scores, cohort = cohort, x = x)
  )
})
