# Tests for optimizing biomarker panels across multiple independent validation cohorts (multi-objective settings).
test_that("optimize_panel handles multiple cohorts", {
  skip_slow_tests()
  skip_if_not(file.exists(fixture_path("fake_gene_expression.Rds")))
  fixture <- readRDS(fixture_path("fake_gene_expression.Rds"))

  res <- optimize_panel(
    x = fixture$x_list,
    y = fixture$y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 4,
    feature_pool = colnames(fixture$x_list[[1]])[seq_len(12)],
    feature_transform = "none",
    nsga_control = list(popSize = 16, maxiter = 15)
  )

  expect_s4_class(res, "OptimizationResult")
  expect_true(n_solutions(res) >= 1)

  # Check feature count constraint
  sols <- solutions(res)
  for (i in seq_len(nrow(sols))) {
    expect_true(length(sols$features[[i]]) <= 4)
  }
})

test_that("optimize_panel intersects feature sets across cohorts", {
  skip_slow_tests()
  make_cohort <- function(seed, cols) {
    set.seed(seed)
    n <- 200L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 0.9 * x[, "gene_common1"] - 0.7 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  c1 <- make_cohort(1, c("gene_common1", "gene_common2", "gene_unique1"))
  c2 <- make_cohort(2, c("gene_common2", "gene_common1", "gene_unique2"))

  res <- optimize_panel(
    x = list(c1$x, c2$x),
    y = list(c1$y, c2$y),
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "none",
    nsga_control = list(popSize = 12, maxiter = 8)
  )

  expect_s4_class(res, "OptimizationResult")

  # All features in solutions should be from the intersection
  sols <- solutions(res)
  all_features <- unique(unlist(sols$features))
  expect_true(all(all_features %in% c("gene_common1", "gene_common2")))
})

test_that("feature_alignment = 'majority' keeps features in >= 50% cohorts", {
  skip_slow_tests()
  # Create 4 cohorts with varying feature overlap
  make_cohort <- function(seed, features) {
    set.seed(seed)
    n <- 200L
    x <- matrix(rnorm(n * length(features)), nrow = n, ncol = length(features))
    colnames(x) <- features
    linear <- 0.8 * x[, "common1"] + 0.5 * x[, "common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohorts <- list(
    make_cohort(1, c("common1", "common2", "partial1", "partial2", "unique1")),
    make_cohort(2, c("common1", "common2", "partial1", "partial2", "unique2")),
    make_cohort(3, c("common1", "common2", "partial1", "unique3")),
    make_cohort(4, c("common1", "common2", "unique4"))
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  expect_warning({
    res <- optimize_panel(
      x = x_list,
      y = y_list,
      objectives = define_objectives(metrics = c("sensitivity", "specificity")),
      max_features = 2,
      feature_transform = "none",
      feature_alignment = "majority",
      nsga_control = list(popSize = 12, maxiter = 8)
    )
  }, regexp = "Imputing")

  expect_s4_class(res, "OptimizationResult")
  expect_equal(res@control$feature_alignment, "majority")

  pool_features <- res@control$feature_pool
  expect_true(all(c("common1", "common2") %in% pool_features))
  expect_true("partial1" %in% pool_features)  # 3/4 = 75%
  expect_true("partial2" %in% pool_features)  # 2/4 = 50%
  expect_false(any(c("unique1", "unique2", "unique3", "unique4") %in% pool_features))
})

test_that("feature_alignment = 'intersection' is default and drops partial features", {
  skip_slow_tests()
  make_cohort <- function(seed, features) {
    set.seed(seed)
    n <- 200L
    x <- matrix(rnorm(n * length(features)), nrow = n, ncol = length(features))
    colnames(x) <- features
    y <- factor(sample(c("No", "Yes"), n, replace = TRUE), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  cohorts <- list(
    make_cohort(1, c("geneA", "geneB", "geneC")),
    make_cohort(2, c("geneA", "geneB", "geneD"))
  )

  x_list <- lapply(cohorts, `[[`, "x")
  y_list <- lapply(cohorts, `[[`, "y")

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 2,
    feature_transform = "none",
    feature_alignment = "intersection",
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  expect_equal(res@control$feature_alignment, "intersection")
  expect_equal(sort(res@control$feature_pool), c("geneA", "geneB"))
})
