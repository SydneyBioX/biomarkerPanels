make_prep_fixture <- function(features_1 = c("geneA", "geneB", "geneC"),
                              features_2 = c("geneC", "geneB", "geneD")) {
  set.seed(42)
  mk <- function(feats) {
    n <- 6L
    m <- matrix(rnorm(n * length(feats)), nrow = n)
    colnames(m) <- feats
    m
  }
  y <- factor(rep(c("No", "Yes"), each = 3L), levels = c("No", "Yes"))
  list(
    x = list(cohort_a = mk(features_1), cohort_b = mk(features_2)),
    y = list(y, y)
  )
}

test_that(".prepare_selection_inputs intersects and preserves cohort-1 order", {
  fx <- make_prep_fixture()

  prepared <- .prepare_selection_inputs(fx$x, fx$y)

  # geneA / geneD are cohort-specific; geneB, geneC are shared. Order follows
  # cohort 1's columns, not alphabetical order.
  expect_identical(prepared$feature_names, c("geneB", "geneC"))
  expect_identical(prepared$cohort_names, c("cohort_a", "cohort_b"))
  expect_true(all(vapply(prepared$matrices, ncol, integer(1)) == 2L))
  expect_identical(colnames(prepared$matrices[[2]]), c("geneB", "geneC"))
})

test_that(".prepare_selection_inputs honours feature_order = 'sorted'", {
  fx <- make_prep_fixture(
    features_1 = c("geneC", "geneB", "geneA"),
    features_2 = c("geneA", "geneB", "geneC")
  )

  prepared <- .prepare_selection_inputs(fx$x, fx$y, feature_order = "sorted")

  expect_identical(prepared$feature_names, c("geneA", "geneB", "geneC"))
})

test_that(".prepare_selection_inputs encodes responses as requested", {
  fx <- make_prep_fixture()

  int_prep <- .prepare_selection_inputs(fx$x, fx$y, response = "integer")
  expect_type(int_prep$responses[[1]], "integer")
  expect_identical(sort(unique(int_prep$responses[[1]])), c(0L, 1L))

  fct_prep <- .prepare_selection_inputs(fx$x, fx$y, response = "factor")
  expect_s3_class(fct_prep$responses[[1]], "factor")
  expect_identical(levels(fct_prep$responses[[1]]), c("No", "Yes"))
})

test_that(".prepare_selection_inputs union keeps ragged matrices", {
  fx <- make_prep_fixture()

  prepared <- .prepare_selection_inputs(fx$x, fx$y, align = "union")

  # Union takes every feature, ordered by first appearance across cohorts.
  expect_identical(prepared$feature_names,
                   c("geneA", "geneB", "geneC", "geneD"))
  # Matrices are deliberately left un-subset under union alignment.
  expect_identical(colnames(prepared$matrices[[1]]),
                   c("geneA", "geneB", "geneC"))
  expect_identical(colnames(prepared$matrices[[2]]),
                   c("geneC", "geneB", "geneD"))
})

test_that(".prepare_selection_inputs reports per-cohort counts on empty overlap", {
  fx <- make_prep_fixture(
    features_1 = c("geneA", "geneB"),
    features_2 = c("geneY", "geneZ")
  )

  expect_error(
    .prepare_selection_inputs(fx$x, fx$y),
    "No shared features were found across cohorts"
  )
  expect_error(
    .prepare_selection_inputs(fx$x, fx$y),
    "cohort_a: 2, cohort_b: 2"
  )
})

test_that(".prepare_selection_inputs validates lengths and sample counts", {
  fx <- make_prep_fixture()

  expect_error(
    .prepare_selection_inputs(fx$x, fx$y[1]),
    "must have the same length"
  )

  bad_y <- fx$y
  bad_y[[2]] <- bad_y[[2]][1:2]
  expect_error(
    .prepare_selection_inputs(fx$x, bad_y),
    "must match the length"
  )
})

test_that(".prepare_selection_inputs accepts a single non-list cohort", {
  fx <- make_prep_fixture()

  prepared <- .prepare_selection_inputs(fx$x[[1]], fx$y[[1]])

  expect_length(prepared$matrices, 1L)
  expect_identical(prepared$feature_names, c("geneA", "geneB", "geneC"))
  expect_identical(prepared$cohort_names, "cohort_01")
})
