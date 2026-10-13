# Core tests for optimize_panel, validating basic structural outputs, bounds checking, and handling of base scenarios.
test_that("Pareto front contains no dominated solutions", {
  skip_slow_tests()
  set.seed(777)
  n <- 40L
  p <- 10L
  x <- matrix(rnorm(n * p), nrow = n)
  colnames(x) <- paste0("gene_", seq_len(p))
  y <- factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes"))
  x[y == "Yes", 1:3] <- x[y == "Yes", 1:3] + 2

  res <- optimize_panel(
    x = x,
    y = y,
    objectives = define_objectives(metrics = c("sensitivity", "specificity", "num_features")),
    max_features = 6,
    feature_pool = colnames(x),
    feature_transform = "none",
    nsga_control = list(popSize = 20, maxiter = 20)
  )

  sols <- solutions(res)
  n_sol <- nrow(sols)
  if (n_sol < 2) skip("Only one Pareto solution; dominance check trivially passes")

  obj_cols <- setdiff(names(sols), c("solution_id", "base_features", "features"))
  directions <- res@control$objective_directions[obj_cols]

  # For each pair of solutions, verify neither dominates the other on ALL objectives
  for (i in seq_len(n_sol - 1)) {
    for (j in (i + 1):n_sol) {
      vals_i <- as.numeric(sols[i, obj_cols])
      vals_j <- as.numeric(sols[j, obj_cols])

      # Compare with direction: "maximize" means higher is better
      i_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vi > vj else vi < vj
      }, vals_i, vals_j, directions)

      i_equal_or_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vi >= vj else vi <= vj
      }, vals_i, vals_j, directions)

      j_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vj > vi else vj < vi
      }, vals_i, vals_j, directions)

      j_equal_or_better <- mapply(function(vi, vj, dir) {
        if (dir == "maximize") vj >= vi else vj <= vi
      }, vals_i, vals_j, directions)

      # i dominates j if i is >= on all and > on at least one
      i_dominates_j <- all(i_equal_or_better) && any(i_better)
      j_dominates_i <- all(j_equal_or_better) && any(j_better)

      expect_false(
        i_dominates_j,
        info = sprintf("Solution %d dominates solution %d", i, j)
      )
      expect_false(
        j_dominates_i,
        info = sprintf("Solution %d dominates solution %d", j, i)
      )
    }
  }
})

# ============================================================================
# fitness_mode dispatch and partitioning
# ============================================================================

test_that("fitness_mode = 'loco' with no partitioning stores no held-out data", {
  skip_slow_tests()
  sim <- simulate_expression_data(p = 20, n = 120, k = 2, seed = 11)

  res <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(sim$x_list[[1]])[seq_len(6)],
    feature_transform = "none",
    fitness_mode = "loco",
    train_ratio = 1,
    seed = 7,
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  expect_null(res@control$heldout_x)
  expect_null(res@control$partition_info)
  # Nothing was split off, so the result carries every input sample.
  expect_equal(nrow(res@aggregated_x), sum(vapply(sim$x_list, nrow, integer(1))))
})

test_that("within_cohort_val partitioning stores held-out data and partition info", {
  skip_slow_tests()
  sim <- simulate_expression_data(p = 20, n = 120, k = 2, seed = 12)

  res <- optimize_panel(
    x = sim$x_list,
    y = sim$y_list,
    objectives = define_objectives(metrics = c("sensitivity", "specificity")),
    max_features = 3,
    feature_pool = colnames(sim$x_list[[1]])[seq_len(6)],
    feature_transform = "none",
    fitness_mode = "within_cohort_val",
    train_ratio = 0.7,
    val_ratio = 0.2,
    seed = 7,
    nsga_control = list(popSize = 12, maxiter = 5)
  )

  n_total <- sum(vapply(sim$x_list, nrow, integer(1)))
  expect_equal(length(res@control$heldout_y), nrow(res@control$heldout_x))
  expect_equal(res@control$partition_info$heldout_ratio, 0.1)
  # aggregated data is train + val; held-out is disjoint from it.
  expect_equal(nrow(res@aggregated_x) + nrow(res@control$heldout_x), n_total)
})
