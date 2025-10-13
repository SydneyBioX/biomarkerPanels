test_that("ensure_binary_response enforces levels", {
  fac <- factor(c("No", "Yes", "No"), levels = c("No", "Yes"))
  expect_identical(ensure_binary_response(fac), fac)

  chr <- c("No", "Yes", "No")
  expect_identical(levels(ensure_binary_response(chr)), c("No", "Yes"))

  expect_error(ensure_binary_response(factor(c("A", "B"))),
               "must be a factor")
})
