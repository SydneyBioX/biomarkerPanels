test_that("optimize_panel returns a BiomarkerPanelResult", {
  set.seed(123)
  sim <- simulate_expression_data(p = 30, n = 25, k = 1, seed = 42)
  x <- sim$x_list[[1]]
  y <- sim$y_list[[1]]

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity", "num_features")),
    max_features = 3,
    feature_pool = colnames(x)[seq_len(8)],
    cohort_aggregator = "none",
    nsga_control = list(popsize = 20, generations = 20)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(length(res@features) <= 3)
  expect_named(res@metrics, c("sensitivity", "specificity", "num_features"))
  expect_s3_class(res@objectives, "data.frame")
  expect_true(all(res@objectives$objective %in% c("sensitivity", "specificity", "num_features")))
  expect_equal(res@training_data$num_cohorts, 1L)
  expect_equal(res@training_data$cohort_labels, "cohort_01")
})

test_that("optimize_panel handles multiple cohorts", {
  skip_if_not(file.exists(fixture_path("fake_gene_expression.Rds")))
  fixture <- readRDS(fixture_path("fake_gene_expression.Rds"))

  res <- optimize_panel(
    x = fixture$x_list,
    y = fixture$y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity", "num_features")),
    max_features = 4,
    feature_pool = colnames(fixture$x_list[[1]])[seq_len(12)],
    cohort_aggregator = "none",
    nsga_control = list(popsize = 16, generations = 15)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(length(res@features) <= 4)
  expect_equal(res@training_data$num_cohorts, length(fixture$x_list))
  expect_true(all(names(res@training_data$cohort_counts) %in% names(fixture$x_list)))
})

test_that("optimize_panel intersects feature sets across cohorts", {
  make_cohort <- function(seed, cols) {
    set.seed(seed)
    n <- 40L
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
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "none",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(all(res@features %in% c("gene_common1", "gene_common2")))
})

test_that("min_metric_constraint builds feasible constraint", {
  constraint <- min_metric_constraint("sensitivity", threshold = 0.9)
  expect_type(constraint$label, "character")
  expect_true(is.function(constraint$fun))
  expect_equal(constraint$loss, "sensitivity")
  expect_equal(constraint$direction, "maximize")
})

test_that("optimize_panel enforces minimum metric constraints", {
  set.seed(5001)
  n <- 60L
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("gene_common1", "gene_common2")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", "gene_common1"] <- x[y == "Yes", "gene_common1"] + 3
  x[y == "No", "gene_common1"] <- x[y == "No", "gene_common1"] - 3

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("specificity")),
    max_features = 1,
    feature_pool = colnames(x),
    constraints = list(min_metric_constraint("sensitivity", threshold = 0.9)),
    cohort_aggregator = "none",
    nsga_control = list(popsize = 16, generations = 15)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true("min_sensitivity_0.9" %in% res@control$constraints)
  eval_res <- evaluate_panel(res, x, y)
  expect_gte(eval_res$metrics["sensitivity"], 0.9)
})

test_that("optimize_panel errors when constraints infeasible", {
  set.seed(5002)
  n <- 40L
  x <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(x) <- c("gene_common1", "gene_common2")
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

  expect_error(
    optimize_panel(
      x = x,
      y = y,
      objectives = define_objectives(losses = c("specificity")),
      max_features = 1,
      feature_pool = colnames(x),
      constraints = list(min_metric_constraint("sensitivity", threshold = 1.01)),
      cohort_aggregator = "none",
      nsga_control = list(popsize = 12, generations = 10)
    ),
    "No solutions satisfied the supplied constraints"
  )
})

test_that("pairwise cohort aggregator produces contrast features", {
  set.seed(1234)
  make_mat <- function(seed) {
    set.seed(seed)
    m <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(m) <- c("A", "B", "C")
    m
  }
  x_list <- list(make_mat(1), make_mat(2))
  y_list <- lapply(1:2, function(i) factor(rep(c("No", "Yes"), length.out = 10), levels = c("No", "Yes")))

  res <- optimize_panel(
    x = x_list,
    y = y_list,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 2,
    cohort_aggregator = "pairwise_ratios",
    nsga_control = list(popsize = 12, generations = 8)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(res@control$cohort_aggregator, "pairwise_ratios")
  expect_true(all(grepl("--", res@control$feature_pool)))
})

test_that("feature_pool accepts base features with pairwise aggregator", {
  set.seed(777)
  x <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(x) <- c("GeneA", "GeneB", "GeneC")
  y <- factor(sample(c("No", "Yes"), 20, replace = TRUE), levels = c("No", "Yes"))

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(losses = c("sensitivity", "specificity")),
    max_features = 1,
    feature_pool = c("GeneA", "GeneC"),
    cohort_aggregator = "pairwise_ratios",
    nsga_control = list(popsize = 12, generations = 5)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_equal(sort(res@control$base_feature_pool), sort(c("GeneA", "GeneC")))
  expect_true(all(grepl("--", res@features, fixed = TRUE)))
  components <- unique(unlist(strsplit(res@features, "--", fixed = TRUE)))
  expect_true(all(components %in% c("GeneA", "GeneC")))
})
