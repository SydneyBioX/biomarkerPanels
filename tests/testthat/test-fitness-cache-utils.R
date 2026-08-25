test_that("shared panel selector matches adaptive and fixed threshold behavior", {
  features <- c("b", "a", "c", "d", "e")
  decision <- c(0.2, 0.9, 0.8, 0.8, 0.1)

  adaptive <- .make_panel_selector(
    feature_pool = features,
    max_features = 3L,
    min_features_required = 1L,
    selection_threshold = "adaptive"
  )(decision)
  expect_equal(adaptive$base_features, c("a", "c", "d"))

  fixed <- .make_panel_selector(
    feature_pool = features,
    max_features = 2L,
    min_features_required = 1L,
    selection_threshold = 0.75
  )(decision)
  expect_equal(fixed$base_features, c("a", "c"))

  fallback <- .make_panel_selector(
    feature_pool = features,
    max_features = 2L,
    min_features_required = 2L,
    selection_threshold = 0.95
  )(decision)
  expect_equal(fallback$base_features, c("a", "c"))
})

test_that("fast binomial glm predictions match glm for no-cohort and cohort cases", {
  set.seed(1001)
  n_train <- 36L
  n_new <- 12L
  x_train <- matrix(rnorm(n_train * 3L), nrow = n_train)
  x_new <- matrix(rnorm(n_new * 3L), nrow = n_new)
  colnames(x_train) <- colnames(x_new) <- paste0("g", 1:3)
  truth <- factor(rep(c("No", "Yes"), each = n_train / 2L),
                  levels = c("No", "Yes"))

  df_train <- data.frame(
    .response = as.integer(truth) - 1L,
    as.data.frame(x_train, check.names = TRUE)
  )
  fit <- stats::glm(.response ~ ., data = df_train, family = stats::binomial())
  df_new <- data.frame(as.data.frame(x_new, check.names = TRUE))
  pred_glm <- as.numeric(stats::predict(fit, newdata = df_new, type = "response"))
  pred_fast <- .fit_predict_binomial_glm(x_train, truth, x_new)
  expect_equal(pred_fast, pred_glm, tolerance = 1e-10)

  cohort_train <- factor(rep(c("A", "B", "C"), length.out = n_train))
  cohort_new <- factor(rep(c("A", "B", "C"), length.out = n_new),
                       levels = levels(cohort_train))
  df_train$.cohort <- cohort_train
  fit_cohort <- stats::glm(.response ~ ., data = df_train,
                           family = stats::binomial())
  df_new$.cohort <- cohort_new
  pred_glm_cohort <- as.numeric(stats::predict(
    fit_cohort, newdata = df_new, type = "response"
  ))
  pred_fast_cohort <- .fit_predict_binomial_glm(
    x_train = x_train,
    truth = truth,
    x_new = x_new,
    cohort_train = cohort_train,
    cohort_new = cohort_new,
    predict_cohort = "observed"
  )
  expect_equal(pred_fast_cohort, pred_glm_cohort, tolerance = 1e-10)

  stored_fit <- .fit_final_model(x_train, truth, cohort_train)
  pred_reference <- .predict_panel_model(stored_fit, x_new, cohort = cohort_new)
  pred_fast_reference <- .fit_predict_binomial_glm(
    x_train = x_train,
    truth = truth,
    x_new = x_new,
    cohort_train = cohort_train,
    cohort_new = cohort_new,
    predict_cohort = "reference"
  )
  expect_equal(pred_fast_reference, pred_reference, tolerance = 1e-10)
})

.make_cache_test_data <- function() {
  set.seed(2002)
  n_train <- 32L
  n_val <- 16L
  p <- 5L
  train_x <- matrix(rnorm(n_train * p), nrow = n_train)
  val_x <- matrix(rnorm(n_val * p), nrow = n_val)
  colnames(train_x) <- colnames(val_x) <- paste0("g", seq_len(p))
  train_y <- factor(rep(c("No", "Yes"), each = n_train / 2L),
                    levels = c("No", "Yes"))
  val_y <- factor(rep(c("No", "Yes"), each = n_val / 2L),
                  levels = c("No", "Yes"))
  train_cohort <- factor(rep(c("A", "B"), length.out = n_train))
  val_cohort <- factor(rep(c("A", "B"), length.out = n_val))
  list(
    train_x = train_x,
    val_x = val_x,
    train_y = train_y,
    val_y = val_y,
    train_cohort = train_cohort,
    val_cohort = val_cohort
  )
}

test_that("fitness cache primitives round-trip panel keys", {
  cache <- .new_fitness_cache()
  key <- .panel_key(c("g2", "g1"))

  # Miss on an empty cache, and a NULL cache is a safe no-op.
  expect_null(.cache_get(cache, key))
  expect_null(.cache_get(NULL, key))
  expect_equal(.cache_set(NULL, key, 1), 1)

  .cache_set(cache, key, c(1, 2))
  expect_equal(.cache_get(cache, key), c(1, 2))

  # Same panel, same key; a different panel is a distinct entry.
  expect_identical(.panel_key(c("g2", "g1")), key)
  expect_null(.cache_get(cache, .panel_key(c("g1", "g2"))))

  # Context namespaces the key.
  ctx_key <- .panel_key(c("g2", "g1"), context = "split1")
  expect_false(identical(ctx_key, key))
  expect_null(.cache_get(cache, ctx_key))

  # Re-setting an existing key overwrites in place.
  .cache_set(cache, key, c(3, 4))
  expect_equal(.cache_get(cache, key), c(3, 4))
})

test_that("validation fitness assigns identical objectives to duplicate panels", {
  dat <- .make_cache_test_data()
  objectives <- define_objectives(metrics = c("sensitivity", "specificity"))
  pop <- rbind(
    c(0.9, 0.8, 0.1, 0.2, 0.3),
    c(0.9, 0.8, 0.1, 0.2, 0.3),
    c(0.1, 0.2, 0.95, 0.85, 0.3)
  )

  fit <- .make_validation_fitness(
    train_x = dat$train_x,
    train_y = dat$train_y,
    train_cohort = dat$train_cohort,
    val_x = dat$val_x,
    val_y = dat$val_y,
    val_cohort = dat$val_cohort,
    feature_pool = colnames(dat$train_x),
    max_features = 2L,
    objectives = objectives,
    constraints = list(),
    regularized = FALSE,
    alpha = 0.5,
    feature_transform = "none",
    min_features_required = 1L,
    selection_threshold = 0.5
  )

  values <- fit$wrapper(pop)
  expect_true(all(is.finite(values)))
  expect_equal(values[1, ], values[2, ], tolerance = 1e-10)
})

test_that("duplicate selected panels are evaluated once when cache is enabled", {
  dat <- .make_cache_test_data()
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  counted_objective <- list(
    mean_score = list(
      label = "mean_score",
      direction = "maximize",
      fun = function(truth, scores, selected = NULL, ...) {
        counter$n <- counter$n + 1L
        mean(scores)
      }
    )
  )
  pop <- rbind(
    c(0.9, 0.8, 0.1, 0.2, 0.3),
    c(0.9, 0.8, 0.1, 0.2, 0.3)
  )

  cached <- .make_validation_fitness(
    train_x = dat$train_x,
    train_y = dat$train_y,
    train_cohort = dat$train_cohort,
    val_x = dat$val_x,
    val_y = dat$val_y,
    val_cohort = dat$val_cohort,
    feature_pool = colnames(dat$train_x),
    max_features = 2L,
    objectives = counted_objective,
    constraints = list(),
    regularized = FALSE,
    alpha = 0.5,
    feature_transform = "none",
    min_features_required = 1L,
    selection_threshold = 0.5
  )
  cached$wrapper(pop)
  expect_equal(counter$n, 1L)
})

test_that("validation fitness cache handles unregularized pairwise ratios", {
  dat <- .make_cache_test_data()
  objectives <- define_objectives(metrics = c("auc", "num_features"))
  pop <- rbind(
    c(0.9, 0.8, 0.7, 0.1, 0.2),
    c(0.9, 0.8, 0.7, 0.1, 0.2),
    c(0.1, 0.2, 0.95, 0.85, 0.8)
  )

  fit <- .make_validation_fitness(
    train_x = dat$train_x,
    train_y = dat$train_y,
    train_cohort = dat$train_cohort,
    val_x = dat$val_x,
    val_y = dat$val_y,
    val_cohort = dat$val_cohort,
    feature_pool = colnames(dat$train_x),
    max_features = 3L,
    objectives = objectives,
    constraints = list(),
    regularized = FALSE,
    alpha = 0.5,
    feature_transform = "pairwise_ratios",
    min_features_required = 2L,
    selection_threshold = 0.5
  )

  values <- fit$wrapper(pop)
  expect_true(all(is.finite(values)))
  expect_equal(values[1, ], values[2, ], tolerance = 1e-10)
  # Repeat evaluation is served from the cache and must be identical.
  expect_equal(fit$wrapper(pop), values, tolerance = 1e-10)
})

test_that("rotating validation cache keys are split-specific and batch-consistent", {
  set.seed(3003)
  pool_x <- matrix(rnorm(24L * 3L), nrow = 24L)
  colnames(pool_x) <- paste0("g", 1:3)
  pool_y <- factor(
    c("No", "Yes", "No", "No", "Yes", "No", "No", "No",
      "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No",
      "Yes", "Yes", "No", "Yes", "No", "Yes", "No", "No"),
    levels = c("No", "Yes")
  )
  pool_cohort <- factor(rep(c("A", "B"), length.out = 24L))
  splits <- list(
    list(train = 9:24, val = 1:8),
    list(train = 1:16, val = 17:24)
  )
  split_objective <- list(
    val_yes = list(
      label = "val_yes",
      direction = "maximize",
      fun = function(truth, scores, selected = NULL, ...) {
        sum(truth == "Yes")
      }
    )
  )
  fit <- .make_rotating_validation_fitness(
    pool_x = pool_x,
    pool_y = pool_y,
    pool_cohort = pool_cohort,
    splits = splits,
    feature_pool = colnames(pool_x),
    max_features = 1L,
    objectives = split_objective,
    constraints = list(),
    regularized = FALSE,
    alpha = 0.5,
    feature_transform = "none",
    min_features_required = 1L,
    selection_threshold = 0.5
  )

  pop <- rbind(
    c(0.9, 0.1, 0.2),
    c(0.1, 0.9, 0.2)
  )
  first <- fit$wrapper(pop)
  second <- fit$wrapper(pop)

  expect_equal(first[, 1], c(-2, -2))
  expect_equal(second[, 1], c(-4, -4))
})
