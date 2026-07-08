# Tests for the consolidated stored-model prediction dispatcher
# (.predict_panel_model). These lock in that glm, cv.glmnet, and the two
# evaluation call-site contracts (expected_features NULL vs supplied) stay
# consistent, since the dispatcher replaced several duplicated inline paths.

make_xy <- function(seed, cols, n = 40L) {
  set.seed(seed)
  x <- matrix(rnorm(n * length(cols)), n, length(cols),
              dimnames = list(NULL, cols))
  lin <- 1.2 * x[, 1] - 0.7 * x[, 2]
  y <- factor(ifelse(runif(n) < stats::plogis(lin), "Yes", "No"),
              levels = c("No", "Yes"))
  list(x = x, y = y)
}

test_that(".predict_panel_model matches a hand-rolled glm prediction", {
  d <- make_xy(1, c("gA", "gB", "gC"))
  df <- data.frame(.response = as.integer(d$y) - 1L,
                   as.data.frame(d$x, check.names = TRUE))
  model <- stats::glm(.response ~ ., data = df, family = stats::binomial())

  new_d <- make_xy(2, c("gA", "gB", "gC"))
  ref_nd <- as.data.frame(new_d$x)
  names(ref_nd) <- make.names(names(ref_nd), unique = TRUE)
  ref <- as.numeric(stats::predict(model, newdata = ref_nd, type = "response"))

  # expected_features = NULL (fitness-path contract)
  expect_equal(.predict_panel_model(model, new_d$x), ref)
  # expected_features supplied (evaluation-path contract) — same values,
  # since the reorder/validation is output-neutral when names already match.
  expect_equal(
    .predict_panel_model(model, new_d$x, expected_features = colnames(new_d$x)),
    ref
  )
})

test_that(".predict_panel_model validates feature names when expected_features supplied", {
  d <- make_xy(1, c("gA", "gB"))
  df <- data.frame(.response = as.integer(d$y) - 1L,
                   as.data.frame(d$x, check.names = TRUE))
  model <- stats::glm(.response ~ ., data = df, family = stats::binomial())
  new_d <- make_xy(2, c("gA", "gB"))

  # A mismatching expected_features must error (evaluation contract);
  # NULL expected_features must not (fitness contract).
  expect_error(
    .predict_panel_model(model, new_d$x, expected_features = c("gA", "WRONG")),
    "Feature name mismatch"
  )
  expect_silent(.predict_panel_model(model, new_d$x))
})

test_that(".predict_panel_model handles cv.glmnet and errors on unknown models", {
  skip_if_not_installed("glmnet")
  d <- make_xy(1, c("gA", "gB", "gC"))
  set.seed(3)
  gm <- glmnet::cv.glmnet(d$x, as.integer(d$y) - 1L,
                          family = "binomial", nfolds = 5)
  gm$biomarkerPanels_meta <- list(cohort_info = NULL)
  new_d <- make_xy(2, c("gA", "gB", "gC"))

  ref <- as.numeric(stats::predict(gm, newx = new_d$x, s = "lambda.min",
                                   type = "response")[, 1])
  expect_equal(.predict_panel_model(gm, new_d$x), ref)

  expect_error(.predict_panel_model(list(a = 1), new_d$x), "Unknown model type")
})
