# Tests for the metric registry framework and objective definition functions that interface with the optimizer.
test_that("build_objectives registers and retrieves correctly", {
  objs <- build_objectives(c("sensitivity", "specificity", "num_features"))

  expect_type(objs, "list")
  expect_equal(length(objs), 3)

  expect_equal(names(objs), c("sensitivity", "specificity", "num_features"))

  for (m in names(objs)) {
    expect_true(is.function(objs[[m]]$fun))
    expect_type(objs[[m]]$direction, "character")
    expect_true(objs[[m]]$direction %in% c("maximize", "minimize"))
  }
})

test_that("define_objectives preserves cutoff_strategy", {
  objs <- define_objectives(metrics = c("sensitivity", "specificity", "num_features"), cutoff_strategy = "prevalence")

  expect_type(objs, "list")
  # Ensure the functions returned by define_objectives encapsulate the strategy
  # We test this by checking if the formals or closure environment holds it,
  # or functionally by injecting a test case

  f_sens <- objs$sensitivity$fun
  truth <- factor(c("No", "Yes", "No", "No", "No", "No", "No", "No", "No", "No"), levels = c("No", "Yes")) # 10% prevalence
  scores <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  # With prevalence = 10%, cutoff is top 10% of scores -> >= 0.1
  # True positive has score 0.2, which is >= 0.1 -> TP=1
  val <- f_sens(truth, scores)
  expect_equal(val, 1)
})

test_that("metric registry allows addition and query of custom metrics", {
  custom_fun <- function(truth, scores, ...) {
    mean(scores)
  }

  register_metric("mean_score", custom_fun, "maximize")

  objs <- build_objectives("mean_score")
  expect_equal(names(objs), "mean_score")
  expect_equal(objs$mean_score$direction, "maximize")

  truth <- factor(c("No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.4, 0.8)
  expect_equal(objs$mean_score$fun(truth, scores), 0.6)

  # Clean up registry
  rm("mean_score", envir = biomarkerPanels:::.metric_registry)
})

test_that("define_objectives constructs proper objective functions", {
  objs <- define_objectives(
    metrics = c("sensitivity", "num_features"),
    cutoff_prob = 0.5
  )

  expect_true(is.list(objs))
  expect_true("sensitivity" %in% names(objs))
  expect_true("num_features" %in% names(objs))

  expect_equal(objs$sensitivity$direction, "maximize")
  expect_equal(objs$num_features$direction, "minimize")
})

test_that("cutoff_dependent metadata identifies metrics correctly via registry", {
  registry <- metric_registry()
  # These metrics require cutoff probability logic
  expect_true(registry$sensitivity$cutoff_dependent)
  expect_true(registry$specificity$cutoff_dependent)
  expect_true(registry$f1$cutoff_dependent)
  expect_true(registry$precision$cutoff_dependent)
  expect_true(registry$npv$cutoff_dependent)
  expect_true(registry$balanced_accuracy$cutoff_dependent)

  # These metrics operate on continuous scores
  expect_false(isTRUE(registry$auc$cutoff_dependent))
  expect_false(isTRUE(registry$min_cohort_auc$cutoff_dependent))
})



test_that("metric_easy_hard_accuracy uses difficulty from args via define_objectives wrapper", {
  difficulty_var <- c("easy", "easy", "hard", "hard")
  # Use params argument to pass difficulty
  objs <- define_objectives(
    metrics = c("easy_hard_accuracy"),
    params = list(easy_hard_accuracy = list(difficulty = difficulty_var))
  )

  # Validate the returned objective function works with the bound arguments
  truth <- factor(c("No", "Yes", "No", "Yes"), levels = c("No", "Yes"))
  scores <- c(0.1, 0.9, 0.2, 0.8)

  result <- objs$easy_hard_accuracy$fun(truth, scores)
  expect_equal(result, 1.0)
})

test_that("build_objectives handles cohort metrics properly", {
  objs <- build_objectives(c("min_cohort_auc", "cohort_auc_gap", "cohort_auc_var", "max_cohort_brier"))
  expect_equal(objs$min_cohort_auc$direction, "maximize")
  expect_equal(objs$cohort_auc_gap$direction, "minimize")
  expect_equal(objs$cohort_auc_var$direction, "minimize")
  expect_equal(objs$max_cohort_brier$direction, "minimize")
})
