# Tests for compute_diversity() — annotation-free co-expression diversity.

.coexpr_modules       <- biomarkerPanels:::.coexpr_modules
.coexpr_diversity     <- biomarkerPanels:::.coexpr_diversity
.center_within_cohort <- biomarkerPanels:::.center_within_cohort

# Build a small OptimizationResult with two internally-correlated, mutually
# independent feature groups (A1..A4 and B1..B4) plus some noise features.
# Optionally plant a between-cohort mean shift to exercise within-cohort
# centering.
make_opt <- function(n = 80L, cohort_shift = 0, seed = 1L) {
  set.seed(seed)
  half <- n %/% 2L
  cohort <- factor(rep(c("c1", "c2"), c(half, n - half)))

  zA <- stats::rnorm(n)
  zB <- stats::rnorm(n)
  noise <- function() stats::rnorm(n, sd = 0.05)

  mat <- cbind(
    A1 = zA + noise(), A2 = zA + noise(), A3 = zA + noise(), A4 = zA + noise(),
    B1 = zB + noise(), B2 = zB + noise(), B3 = zB + noise(), B4 = zB + noise(),
    N1 = stats::rnorm(n), N2 = stats::rnorm(n)
  )

  # Plant a cohort batch shift on the A block: it inflates pooled correlation
  # of A features with anything sharing the shift, but vanishes after
  # within-cohort centering.
  if (cohort_shift != 0) {
    shift <- ifelse(cohort == "c2", cohort_shift, 0)
    mat[, c("A1", "A2", "A3", "A4")] <- mat[, c("A1", "A2", "A3", "A4")] + shift
    mat[, c("N1", "N2")] <- mat[, c("N1", "N2")] + shift
  }

  solutions <- data.frame(solution_id = 1:3, stringsAsFactors = FALSE)
  solutions$base_features <- I(list(
    c("A1", "A2", "A3", "A4"),          # redundant: one module
    c("A1", "B1", "N1", "N2"),          # spread across modules
    c("A1", "A2", "B1", "B2")           # in between
  ))
  solutions$features <- solutions$base_features

  new("OptimizationResult",
      solutions = solutions,
      feature_pool = colnames(mat),
      control = list(),
      training_signature = list(),
      aggregated_x = mat,
      aggregated_y = factor(rep(c("No", "Yes"), length.out = n),
                            levels = c("No", "Yes")),
      aggregated_cohort = cohort)
}

test_that("compute_diversity appends one column, preserves rows and order", {
  opt <- make_opt()
  summary_df <- summarize_solutions(opt)
  out <- compute_diversity(summary_df, opt, n_modules = 4L)

  expect_s3_class(out, "data.frame")
  expect_true("coexpr_diversity" %in% names(out))
  expect_equal(nrow(out), nrow(summary_df))
  expect_identical(out$solution_id, summary_df$solution_id)
  # existing columns untouched
  expect_identical(out[names(summary_df)], summary_df)
})

test_that("scores lie in (0, 1] and a redundant panel scores below a spread one", {
  opt <- make_opt()
  out <- compute_diversity(summarize_solutions(opt), opt, n_modules = 4L)
  d <- out$coexpr_diversity

  expect_true(all(d > 0 & d <= 1))
  # solution 1 = all-A (one module) should be the least diverse;
  # solution 2 = A/B/N1/N2 spread should be the most diverse.
  expect_lt(d[1], d[2])
  expect_equal(d[1], 0.25)  # 1 distinct module / 4 features
})

test_that("OptimizationResult can be passed directly", {
  opt <- make_opt()
  direct <- compute_diversity(opt, n_modules = 4L)
  viadf  <- compute_diversity(summarize_solutions(opt), opt, n_modules = 4L)
  expect_equal(direct$coexpr_diversity, viadf$coexpr_diversity)
})

test_that("n_modules larger than the pool is capped with a warning", {
  opt <- make_opt()
  expect_warning(
    out <- compute_diversity(summarize_solutions(opt), opt, n_modules = 999L),
    "capping"
  )
  expect_true(all(is.finite(out$coexpr_diversity)))
})

test_that("within-cohort centering removes a planted batch shift", {
  # A1 (signal block) and N1 (noise) are independent, but a large between-cohort
  # shift applied to both makes them look strongly co-expressed in the pooled
  # correlation. Within-cohort centering should strip that spurious correlation.
  opt <- make_opt(cohort_shift = 8, seed = 7L)
  mat <- opt@aggregated_x

  cor_naive <- stats::cor(mat[, "A1"], mat[, "N1"])
  centered <- .center_within_cohort(mat, opt@aggregated_cohort)
  cor_centered <- stats::cor(centered[, "A1"], centered[, "N1"])

  expect_gt(abs(cor_naive), 0.7)      # batch shift fakes co-expression
  expect_lt(abs(cor_centered), 0.3)   # centering restores independence

  # And the recovered modules keep A separate from N under centering.
  modules <- .coexpr_modules(mat, opt@aggregated_cohort, n_modules = 4L,
                             linkage = "ward.D2", cor_method = "pearson")
  expect_false(modules[["A1"]] == modules[["N1"]])
})

test_that("constant features are dropped with a warning", {
  opt <- make_opt()
  opt@aggregated_x[, "N2"] <- 1  # zero variance
  expect_warning(
    compute_diversity(summarize_solutions(opt), opt, n_modules = 4L),
    "constant feature"
  )
})

test_that("name collision is rejected", {
  opt <- make_opt()
  summary_df <- summarize_solutions(opt)
  summary_df$coexpr_diversity <- 0
  expect_error(
    compute_diversity(summary_df, opt),
    "already exists"
  )
})

test_that("missing training matrix yields NA with a warning", {
  opt <- make_opt()
  opt@aggregated_x <- matrix(numeric(0), 0, 0)
  expect_warning(
    out <- compute_diversity(summarize_solutions(opt), opt),
    "no training matrix"
  )
  expect_true(all(is.na(out$coexpr_diversity)))
})

test_that(".coexpr_diversity handles unmapped / empty panels", {
  modules <- c(A = 1L, B = 2L, C = 2L)
  expect_equal(.coexpr_diversity(c("A", "B"), modules), 1.0)   # 2 modules / 2
  expect_equal(.coexpr_diversity(c("B", "C"), modules), 0.5)   # 1 module / 2
  expect_true(is.na(.coexpr_diversity(c("Z", "Y"), modules)))  # none map
  expect_true(is.na(.coexpr_diversity(character(0), modules)))
})
