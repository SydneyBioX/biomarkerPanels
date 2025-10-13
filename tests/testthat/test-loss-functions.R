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
                           params = list(sensitivity = list(threshold = 0.7)))
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
