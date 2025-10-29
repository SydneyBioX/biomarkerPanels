test_that("evaluate_panel scores held-out multi-cohort data", {
  panel_features <- c("gene_common1", "gene_common2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, c_index = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list()
  )

  make_validation <- function(seed, cols) {
    set.seed(seed)
    n <- 30L
    x <- matrix(rnorm(n * length(cols)), nrow = n, ncol = length(cols))
    colnames(x) <- cols
    linear <- 1.1 * x[, "gene_common1"] - 0.6 * x[, "gene_common2"]
    prob <- stats::plogis(linear)
    y <- factor(ifelse(runif(n) < prob, "Yes", "No"), levels = c("No", "Yes"))
    list(x = x, y = y)
  }

  v1 <- make_validation(101, c("gene_common1", "gene_common2", "gene_extra1"))
  v2 <- make_validation(202, c("gene_common2", "gene_common1", "gene_extra2"))

  eval_res <- evaluate_panel(
    panel = panel,
    x = list(v1$x, v2$x),
    y = list(v1$y, v2$y)
  )

  expect_named(eval_res$metrics, c("sensitivity", "specificity", "c_index"))
  expect_true(all(is.finite(eval_res$metrics)))
  expect_true(is.matrix(eval_res$confusion))
  expect_equal(dim(eval_res$confusion), c(2, 2))
  expect_true(all(c("No", "Yes") %in% rownames(eval_res$confusion)))
  expect_true(all(c("No", "Yes") %in% colnames(eval_res$confusion)))
  expect_true(is.list(eval_res$roc))
  expect_true(all(c("curve", "highlight", "auc", "plot", "evalm") %in% names(eval_res$roc)))
  expect_s3_class(eval_res$roc$curve, "data.frame")
  expect_true(all(c("threshold", "tpr", "fpr", "sensitivity", "specificity") %in% colnames(eval_res$roc$curve)))
  expect_s3_class(eval_res$roc$highlight, "data.frame")
  expect_equal(nrow(eval_res$roc$highlight), 1L)
  expect_true(is.null(eval_res$roc$plot) || inherits(eval_res$roc$plot, "ggplot"))
  expect_s3_class(eval_res$roc$evalm, "data.frame")
  expect_true(all(c("obs", "prob") %in% colnames(eval_res$roc$evalm)))
})

test_that("evaluate_panel requires column names for single-cohort data", {
  panel_features <- c("gene_common1", "gene_common2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, c_index = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list()
  )

  x <- matrix(rnorm(20), nrow = 10)
  y <- factor(rep(c("No", "Yes"), each = 5), levels = c("No", "Yes"))

  expect_error(
    evaluate_panel(panel, x, y),
    "Selected feature\\(s\\) not found"
  )
})

test_that("evaluate_panel accepts custom cutoff for confusion and ROC highlight", {
  panel_features <- c("g1", "g2")
  metrics <- c(sensitivity = 0.5, specificity = 0.5, c_index = 0.5)
  objective_df <- data.frame(
    solution_id = 1L,
    objective = names(metrics),
    value = metrics,
    direction = c("maximize", "maximize", "maximize"),
    stringsAsFactors = FALSE
  )
  objective_df$features <- I(rep(list(panel_features), nrow(objective_df)))
  panel <- new(
    "BiomarkerPanelResult",
    features = panel_features,
    metrics = metrics,
    objectives = objective_df,
    control = list(),
    training_data = list()
  )

  set.seed(123)
  x <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("g1", "g2")))
  linear <- 1.2 * x[, "g1"] - 0.4 * x[, "g2"]
  prob <- stats::plogis(linear)
  y <- factor(ifelse(runif(20) < prob, "Yes", "No"), levels = c("No", "Yes"))

  eval_res <- evaluate_panel(
    panel = panel,
    x = x,
    y = y,
    cutoff_prob = 0.3
  )

  highlight <- eval_res$roc$highlight
  expect_equal(highlight$threshold, 0.3)
  expect_equal(highlight$sensitivity, loss_sensitivity(y, eval_res$scores, cutoff_prob = 0.3))
  expect_equal(highlight$specificity, loss_specificity(y, eval_res$scores, cutoff_prob = 0.3))
})
