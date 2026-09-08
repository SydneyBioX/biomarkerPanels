# .build_pareto_solutions_df(): post-processing of the raw NSGA population.

setClass("FakeNsgaResult", representation(population = "matrix", front = "numeric"))

# Genome = weight per pool feature; a panel is the set of features with
# weight > 0.5. Metrics are a deterministic function of the panel so identical
# panels always re-evaluate to identical rows.
fake_pool <- c("A", "B", "C", "D")
fake_evaluate <- function(decision_vec) {
  feats <- fake_pool[decision_vec > 0.5]
  list(
    feasible = length(feats) > 0L,
    metrics = c(auc = 0.5 + 0.1 * length(feats), num_features = length(feats)),
    base_features = feats,
    features = feats
  )
}
fake_objectives <- list(auc = list(), num_features = list())
fake_directions <- c(auc = "maximize", num_features = "minimize")

test_that("duplicate panels on the rank-1 front collapse to one solution", {
  # Four distinct genomes; the first three all decode to {A, B} (ordering of
  # weights differs, so the decoded feature order differs too). The fourth
  # decodes to {A}, which is non-dominated (fewer features).
  pop <- rbind(
    c(0.9, 0.8, 0.1, 0.2),
    c(0.7, 0.95, 0.3, 0.0),
    c(0.6, 0.6, 0.4, 0.1),
    c(0.9, 0.1, 0.1, 0.1)
  )
  res <- new("FakeNsgaResult", population = pop, front = rep(1, 4))

  df <- .build_pareto_solutions_df(res, fake_evaluate, fake_objectives, fake_directions)

  expect_equal(nrow(df), 2L)
  expect_equal(df$solution_id, 1:2)
  keys <- vapply(df$base_features, function(f) paste(sort(f), collapse = "|"), "")
  expect_false(any(duplicated(keys)))
  expect_setequal(keys, c("A|B", "A"))
})

test_that("dedupe ignores feature order within a panel", {
  # Same set, opposite decode order (feature order in base_features follows
  # genome order in the fake evaluator only through fake_pool, so build the
  # reversed panel explicitly via a custom evaluator).
  evaluate_ordered <- function(decision_vec) {
    feats <- if (decision_vec[1] > 0.5) c("A", "B") else c("B", "A")
    list(feasible = TRUE, metrics = c(auc = 0.7, num_features = 2),
         base_features = feats, features = feats)
  }
  pop <- rbind(c(1, 0), c(0, 1))
  res <- new("FakeNsgaResult", population = pop, front = c(1, 1))

  df <- .build_pareto_solutions_df(res, evaluate_ordered, fake_objectives, fake_directions)
  expect_equal(nrow(df), 1L)
})

test_that("distinct non-dominated panels are all retained", {
  pop <- rbind(
    c(0.9, 0.1, 0.1, 0.1),  # {A}:       auc 0.6, 1 feature
    c(0.9, 0.9, 0.1, 0.1),  # {A, B}:    auc 0.7, 2 features
    c(0.9, 0.9, 0.9, 0.1)   # {A, B, C}: auc 0.8, 3 features
  )
  res <- new("FakeNsgaResult", population = pop, front = rep(1, 3))
  df <- .build_pareto_solutions_df(res, fake_evaluate, fake_objectives, fake_directions)
  expect_equal(nrow(df), 3L)
  expect_equal(df$num_features, c(1, 2, 3))
})
