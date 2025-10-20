make_ratio_fixture <- function() {
  make_cohort <- function(shift, seed) {
    n <- 8L
    p <- 4L
    set.seed(seed)
    x <- matrix(rnorm(n * p, sd = 0.05), nrow = n, ncol = p)
    colnames(x) <- paste0("gene", seq_len(p))
    y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
    x[, 1] <- x[, 1] + rep(c(0, shift), each = n / 2)
    x[, 2] <- x[, 2] + 0.2
    x[, 3] <- x[, 3] + rep(c(0, shift / 4), each = n / 2)
    baseline <- seq(-0.2, 0.2, length.out = n)
    x[, 4] <- x[, 4] + baseline
    list(x = x, y = y)
  }

  list(
    cohorts = list(
      make_cohort(2.5, 1),
      make_cohort(2.0, 2)
    )
  )
}

test_that("get_top_de_features ranks shared signals", {
  fixture <- make_ratio_fixture()

  selected <- get_top_de_features(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_features = 3,
    combination_method = "Stouffer"
  )

  expect_type(selected, "character")
  expect_lte(length(selected), 3)
  expect_true("gene1" %in% selected)
})

test_that("select_features_for_ratios returns stable and informative sets", {
  fixture <- make_ratio_fixture()

  ratios <- select_features_for_ratios(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_stable = 2,
    n_informative = 3,
    stability_method = "precision_weighted",
    combination_method = "Stouffer"
  )

  expect_named(ratios, c("stable", "informative"))
  expect_length(ratios$stable, 2)
  expect_length(ratios$informative, 3)
  expect_true("gene2" %in% ratios$stable)
  expect_true("gene1" %in% ratios$informative)

  alt <- select_features_for_ratios(
    lapply(fixture$cohorts, `[[`, "x"),
    lapply(fixture$cohorts, `[[`, "y"),
    n_stable = 2,
    n_informative = 2,
    stability_method = "inverse_t_se"
  )
  expect_length(alt$stable, 2)
  expect_length(alt$informative, 2)

  expect_error(
    select_features_for_ratios(
      lapply(fixture$cohorts, `[[`, "x"),
      fixture$cohorts[[1]]$y,
      n_stable = 2,
      n_informative = 2
    ),
    "`x_list` and `y_list` must have the same length"
  )

  expect_error(
    select_features_for_ratios(
      lapply(fixture$cohorts, `[[`, "x"),
      lapply(fixture$cohorts, `[[`, "y"),
      stability_method = "unknown-method"
    ),
    "should be one of"
  )
})

test_that("select_transferable_features prioritises consistent ridge coefficients", {
  make_cohort <- function(seed, delta) {
    set.seed(seed)
    n <- 80L
    x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(x) <- c("geneA", "geneB", "geneC")
    linear <- 1.5 * x[, "geneA"] - 0.5 * x[, "geneB"] + delta
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohort1 <- make_cohort(101, -0.2)
  cohort2 <- make_cohort(202, 0.2)

  set.seed(999)
  transferable <- select_transferable_features(
    list(cohort1$x, cohort2$x),
    list(cohort1$y, cohort2$y),
    n_features = 2,
    lambda_choice = "lambda_min"
  )

  expect_type(transferable$features, "character")
  expect_true("geneA" %in% transferable$features)
  expect_true(all(transferable$scores$sign_consistent))
  expect_true(all(transferable$scores$min_abs >= 0))
  expect_equal(colnames(transferable$coefficients), c("cohort_01", "cohort_02"))

  expect_error(
    select_transferable_features(
      list(cohort1$x, cohort2$x),
      cohort1$y
    ),
    "`x_list` and `y_list` must have the same length"
  )
})

test_that("select_transferable_features intersects cohort feature sets", {
  make_dataset <- function(seed, cols) {
    set.seed(seed)
    n <- 60L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 1.2 * x[, "gene_common1"] - 0.8 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  ds1 <- make_dataset(111, c("gene_common1", "gene_common2", "gene_unique1"))
  ds2 <- make_dataset(222, c("gene_unique2", "gene_common2", "gene_common1"))

  result <- select_transferable_features(
    list(ds1$x, ds2$x),
    list(ds1$y, ds2$y),
    n_features = 2,
    lambda = 0.1,
    min_coefficient = 0.01,
    standardize = FALSE
  )

  expect_true(length(result$features) <= 2)
  expect_true(all(result$features %in% c("gene_common1", "gene_common2")))
  expect_equal(colnames(result$coefficients), c("cohort_01", "cohort_02"))
  expect_true(all(rownames(result$coefficients) %in% c("gene_common1", "gene_common2")))
})
