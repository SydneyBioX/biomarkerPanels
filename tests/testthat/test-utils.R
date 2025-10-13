test_that("ensure_binary_response coerces and orders levels", {
  fac <- factor(c("No", "Yes", "No"), levels = c("No", "Yes"))
  expect_identical(ensure_binary_response(fac), fac)

  chr <- c("No", "Yes", "No")
  expect_identical(levels(ensure_binary_response(chr)), c("No", "Yes"))

  expect_error(
    ensure_binary_response(c("No", "Maybe", "Yes")),
    "exactly two outcome values"
  )

  expect_warning(
    coerced <- ensure_binary_response(factor(c("A", "B"))),
    "Coercing"
  )
  expect_identical(levels(coerced), c("No", "Yes"))
  expect_identical(as.character(coerced), c("No", "Yes"))

  custom <- ensure_binary_response(
    c("control", "FGR", "FGR"),
    positive = "FGR",
    negative = "control"
  )
  expect_identical(levels(custom), c("No", "Yes"))
  expect_identical(as.character(custom), c("No", "Yes", "Yes"))
})
