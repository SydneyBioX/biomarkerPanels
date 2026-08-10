# Fixture: pooled matrix with known signal / batch / null gene structure.
# Signal genes track the outcome only, batch genes track the cohort only.
make_batch_fixture <- function() {
  make_cohort <- function(seed, batch_shift) {
    set.seed(seed)
    n <- 30L
    p_signal <- 5L
    p_batch <- 5L
    p_null <- 10L
    p <- p_signal + p_batch + p_null

    x <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
    colnames(x) <- c(
      paste0("signal_", seq_len(p_signal)),
      paste0("batch_", seq_len(p_batch)),
      paste0("null_", seq_len(p_null))
    )
    y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))

    for (j in seq_len(p_signal)) {
      x[, j] <- x[, j] + rep(c(-1, 1), each = n / 2) * (1.5 + 0.2 * j)
    }
    for (j in seq_len(p_batch)) {
      x[, p_signal + j] <- x[, p_signal + j] + batch_shift * (2.5 + 0.3 * j)
    }

    list(x = x, y = y)
  }

  c1 <- make_cohort(101, -1.5)
  c2 <- make_cohort(202, 0.0)
  c3 <- make_cohort(303, 1.5)

  list(
    x = rbind(c1$x, c2$x, c3$x),
    y = factor(c(as.character(c1$y), as.character(c2$y), as.character(c3$y)),
      levels = c("No", "Yes")
    ),
    cohort = factor(rep(c("A", "B", "C"), each = nrow(c1$x))),
    signal_genes = paste0("signal_", 1:5),
    batch_genes = paste0("batch_", 1:5),
    null_genes = paste0("null_", 1:10)
  )
}

select_batch <- function(fx, n_features = 5L, ...) {
  select_batch_associated_features(
    x = fx$x, cohort = fx$cohort, y = fx$y,
    n_features = n_features, n_pcs = 10L, ...
  )
}

test_that("select_batch_associated_features returns the documented structure", {
  fx <- make_batch_fixture()
  result <- select_batch(fx)

  expect_named(result, c(
    "features", "biology_epsilon", "diagnostics", "problematic_cohort",
    "one_vs_rest", "pure_batch_pcs", "gene_scores", "pca"
  ))
  expect_length(result$features, 5L)
  expect_s3_class(result$pca, "prcomp")
  expect_length(select_batch(fx, n_features = 500L)$features, ncol(fx$x))
})

test_that("selection recovers batch genes and noise cannot inflate the ratio", {
  fx <- make_batch_fixture()
  result <- select_batch(fx)

  # Regression: with a fixed tiny epsilon, null genes whose batch AND biology
  # scores were both numerical noise produced enormous ratios and displaced
  # real batch genes from the top of the ranking.
  expect_setequal(result$features, fx$batch_genes)

  ratios <- stats::setNames(result$gene_scores$ratio, result$gene_scores$gene)
  expect_gt(min(ratios[fx$batch_genes]), max(ratios[fx$null_genes]))

  # Epsilon is anchored to the biology-score distribution, not a constant.
  expect_equal(
    result$biology_epsilon,
    unname(stats::quantile(result$gene_scores$biology_score, 0.5))
  )
})

test_that("biology_penalty_quantile = 1 reduces ranking to batch_score alone", {
  fx <- make_batch_fixture()
  gs <- select_batch(fx)$gene_scores

  expect_identical(
    select_batch(fx, biology_penalty_quantile = 1.0)$features,
    gs$gene[order(-gs$batch_score)][1:5]
  )
})

test_that("y = NULL falls back to a positive epsilon and batch-only ranking", {
  fx <- make_batch_fixture()
  result <- select_batch_associated_features(
    x = fx$x, cohort = fx$cohort, y = NULL, n_features = 5L, n_pcs = 10L
  )

  # Every biology score is zero, so the quantile degenerates; the fallback must
  # keep `ratio` finite (not 0/0 = NaN) and monotone in batch_score.
  expect_true(all(result$gene_scores$biology_score == 0))
  expect_gt(result$biology_epsilon, 0)
  expect_true(all(is.finite(result$gene_scores$ratio)))
  gs <- result$gene_scores
  expect_identical(result$features, gs$gene[order(-gs$batch_score)][1:5])
})

test_that("diagnostics identify pure-batch PCs and the separated cohort", {
  fx <- make_batch_fixture()
  result <- select_batch(fx)

  # The strict criterion should fire on clean data, not one of the fallbacks.
  expect_true(length(result$pure_batch_pcs) > 0L)
  pure <- result$diagnostics[result$diagnostics$is_pure_batch, ]
  expect_true(all(pure$batch_pvalue < 0.10))
  expect_true(all(pure$biology_pvalue > 0.20))

  expect_equal(nrow(result$one_vs_rest), nlevels(fx$cohort))
  # Cohort B has no batch shift, so it should not be the extreme one.
  expect_false(identical(result$problematic_cohort, "B"))
})

test_that("input validation rejects malformed arguments", {
  fx <- make_batch_fixture()

  expect_error(
    select_batch_associated_features(x = "not a matrix", cohort = fx$cohort),
    "must be a matrix or data.frame"
  )

  unnamed <- fx$x
  colnames(unnamed) <- NULL
  expect_error(
    select_batch_associated_features(x = unnamed, cohort = fx$cohort),
    "must have column names"
  )

  expect_error(
    select_batch_associated_features(x = fx$x, cohort = fx$cohort[-1]),
    "same length as the number of rows"
  )

  expect_error(
    select_batch_associated_features(x = fx$x, cohort = fx$cohort, y = fx$y[-1]),
    "same length as the number of rows"
  )

  expect_error(
    select_batch(fx, biology_penalty_quantile = 1.5),
    "biology_penalty_quantile"
  )
})

test_that("select_denominator_features is deprecated but still works", {
  fx <- make_batch_fixture()

  expect_warning(
    old <- select_denominator_features(
      fx$x, fx$cohort, fx$y, n_denominators = 5L, n_pcs = 10L
    ),
    "deprecated"
  )

  # Old name keeps its argument order and exposes the result under both the
  # legacy and current element names.
  expect_identical(old$denominators, select_batch(fx)$features)
  expect_identical(old$features, old$denominators)
})
