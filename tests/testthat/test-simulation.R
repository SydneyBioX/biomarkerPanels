test_that("simulate_expression_data yields expected structure", {
  data <- simulate_expression_data(p = 50, n = 10, k = 2, seed = 1)
  expect_named(data, c("x_list", "y_list", "metadata"))
  expect_length(data$x_list, 2)
  expect_length(data$y_list, 2)
  expect_equal(dim(data$x_list[[1]]), c(10, 50))
  expect_s3_class(data$y_list[[1]], "factor")
})
