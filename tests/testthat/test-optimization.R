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
    nsga_control = list(popsize = 16, generations = 15)
  )

  expect_s4_class(res, "BiomarkerPanelResult")
  expect_true(length(res@features) <= 4)
  expect_equal(res@training_data$num_cohorts, length(fixture$x_list))
  expect_true(all(names(res@training_data$cohort_counts) %in% names(fixture$x_list)))
})
