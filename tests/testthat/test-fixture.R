test_that("fixture dataset loads", {
  path <- fixture_path("fake_gene_expression.Rds")
  skip_if_not(file.exists(path))
  obj <- readRDS(path)
  expect_true(is.list(obj))
  expect_length(obj$x_list, 4)
})
