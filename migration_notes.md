# mco → rmoo Migration Notes

## Executive Summary

This document maps all `mco` package usage in biomarkerPanels and provides the translation to `rmoo`. The migration is straightforward with one main structural change: mco returns a list, while rmoo returns an S4 object.

---

## Files Requiring Changes

### 1. `R/optimize_panel.R` (Primary)

| Line(s) | Current Code | Migration Action |
|---------|-------------|------------------|
| 3 | Roxygen: `[mco::nsga2()]` | Change to `[rmoo::nsga2()]` |
| 47 | Roxygen: `[mco::nsga2()]` | Change to `[rmoo::nsga2()]` |
| 48-49 | Roxygen param defaults | Update parameter names in docs |
| 98-101 | `requireNamespace("mco", ...)` | Change to `requireNamespace("rmoo", ...)` |
| 370-379 | Parameter construction with mco names | Translate parameter names |
| 389 | `do.call(mco::nsga2, nsga_params)` | Change to `do.call(rmoo::nsga2, nsga_params)` |
| 391-393 | `nsga_result$par` extraction | Change to `nsga_result@population` with front filtering |
| 395-398 | Iterating over `nsga_result$par` rows | Filter for `nsga_result@front == 1` first |

**Note**: No `@importFrom mco` tags exist (confirmed via grep).

### 2. `R/cohort_preparation.R`

| Line(s) | Current Code | Migration Action |
|---------|-------------|------------------|
| 361-380 | `.get_adaptive_nsga_defaults()` | Rename parameters: `popsize`→`popSize`, `generations`→`maxiter`, `cprob`→`pcrossover`, `mprob`→`pmutation` |

### 3. `DESCRIPTION`

| Line | Current | Migration Action |
|------|---------|------------------|
| 24 | `mco,` | Replace with `rmoo,` |

### 4. `tests/testthat/test-optimization.R`

| Line(s) | Current Code | Migration Action |
|---------|-------------|------------------|
| 15, 39, 70, etc. | `nsga_control = list(popsize = ...)` | Update parameter names |
| 325-326 | `res@control$nsga2$popsize` assertions | Update parameter names in assertions |

### 5. Auto-Generated Files (Regenerate After Changes)

These files contain mco references but are auto-generated:

| File | Action |
|------|--------|
| `man/optimize_panel.Rd` | Regenerate via `devtools::document()` |
| `doc/*.html` | Regenerate via `devtools::build_vignettes()` |

---

## Parameter Mapping Table

| mco Parameter | rmoo Parameter | Notes |
|---------------|----------------|-------|
| `fn` | `fitness` | Same signature: `function(x)` returning numeric vector |
| `idim` | (implicit) | Inferred from `length(lower)` |
| `odim` | `nObj` | Number of objectives |
| `lower.bounds` | `lower` | Vector of lower bounds |
| `upper.bounds` | `upper` | Vector of upper bounds |
| `popsize` | `popSize` | Population size (note: camelCase) |
| `generations` | `maxiter` | Number of generations |
| `cprob` | `pcrossover` | Crossover probability (default: 0.8) |
| `mprob` | `pmutation` | Mutation probability (default: 0.1) |
| `cdist` | - | rmoo uses SBX crossover with `nc=20` by default |
| `mdist` | - | rmoo uses polynomial mutation with `nm=0.2` by default |
| (none) | `type` | **Required**: `"real-valued"`, `"binary"`, or `"permutation"` |
| (none) | `monitor` | Set to `FALSE` for non-interactive use |
| (none) | `summary` | Set to `FALSE` to disable generation summaries |

---

## Result Structure Changes

### mco::nsga2 returns a list:
```r
nsga_result$par    # Matrix of Pareto-optimal decision variables (all rank-1)
nsga_result$value  # Matrix of corresponding objective values
```

### rmoo::nsga2 returns an S4 object:
```r
nsga_result@population  # Matrix of ALL individuals (not just Pareto-optimal)
nsga_result@fitness     # Matrix of objective values for ALL individuals
nsga_result@front       # Vector of Pareto front ranks (1 = optimal, 2 = next front, etc.)
```

### Critical Difference:
**mco returns only Pareto-optimal solutions in `$par`**, but **rmoo returns the full population**. You must filter by `@front == 1` to get Pareto-optimal solutions:

```r
# OLD (mco)
pareto_solutions <- nsga_result$par

# NEW (rmoo)
optimal_idx <- which(nsga_result@front == 1)
pareto_solutions <- nsga_result@population[optimal_idx, , drop = FALSE]
```

---

## Detailed Code Changes

### optimize_panel.R Line 370-379: Parameter Construction

**BEFORE (mco):**
```r
nsga_params <- c(
  list(
    fn = objective_wrapper,
    idim = decision_dim,
    odim = length(objectives),
    lower.bounds = rep(0, decision_dim),
    upper.bounds = rep(1, decision_dim)
  ),
  nsga_args
)
```

**AFTER (rmoo):**
```r
nsga_params <- c(
  list(
    type = "real-valued",
    fitness = objective_wrapper,
    nObj = length(objectives),
    lower = rep(0, decision_dim),
    upper = rep(1, decision_dim),
    monitor = FALSE,
    summary = FALSE
  ),
  nsga_args
)
```

### optimize_panel.R Line 389: Function Call

**BEFORE:**
```r
nsga_result <- do.call(mco::nsga2, nsga_params)
```

**AFTER:**
```r
nsga_result <- do.call(rmoo::nsga2, nsga_params)
```

### optimize_panel.R Lines 391-398: Result Extraction

**BEFORE (mco):**
```r
if (is.null(dim(nsga_result$par))) {
  nsga_result$par <- matrix(nsga_result$par, nrow = 1)
}

solutions <- lapply(seq_len(nrow(nsga_result$par)), function(i) {
  decision_vec <- nsga_result$par[i, ]
  evaluate_candidate(decision_vec)
})
```

**AFTER (rmoo):**
```r
# Filter to Pareto-optimal solutions (front rank == 1)
optimal_idx <- which(nsga_result@front == 1)
pareto_pop <- nsga_result@population[optimal_idx, , drop = FALSE]

# Handle edge case: single solution returns vector instead of matrix
if (is.null(dim(pareto_pop))) {
  pareto_pop <- matrix(pareto_pop, nrow = 1)
}

solutions <- lapply(seq_len(nrow(pareto_pop)), function(i) {
  decision_vec <- pareto_pop[i, ]
  evaluate_candidate(decision_vec)
})
```

### cohort_preparation.R: Adaptive Defaults

**BEFORE (mco parameter names):**
```r
.get_adaptive_nsga_defaults <- function(n_features) {
  base_params <- list(
    cprob = 0.7,
    cdist = 5,
    mprob = 0.2,
    mdist = 10
  )

  if (n_features <= 30) {
    adaptive_params <- list(popsize = 64, generations = 60)
  } else if (n_features <= 100) {
    adaptive_params <- list(popsize = 128, generations = 150)
  } else {
    adaptive_params <- list(popsize = 200, generations = 300)
  }

  c(adaptive_params, base_params)
}
```

**AFTER (rmoo parameter names):**
```r
.get_adaptive_nsga_defaults <- function(n_features) {
  base_params <- list(
    pcrossover = 0.7,
    pmutation = 0.2
  )
  # Note: rmoo uses SBX crossover (nc=20) and polynomial mutation (nm=0.2) by default
  # cdist/mdist equivalents are set via custom crossover/mutation functions if needed

  if (n_features <= 30) {
    adaptive_params <- list(popSize = 64, maxiter = 60)
  } else if (n_features <= 100) {
    adaptive_params <- list(popSize = 128, maxiter = 150)
  } else {
    adaptive_params <- list(popSize = 200, maxiter = 300)
  }

  c(adaptive_params, base_params)
}
```

---

## Handling Distribution Indices (cdist, mdist)

mco uses `cdist` and `mdist` for crossover and mutation distribution indices. rmoo handles these differently:

- **SBX Crossover** (rmoo default for real-valued): Uses `nc` parameter in `nsgareal_sbxCrossover(object, parents, nc = 20)`
- **Polynomial Mutation** (rmoo default): Uses `nm` parameter in `nsgareal_polMutation(object, parent, nm = 0.20)`

If custom distribution indices are needed, use `nsgaControl()`:
```r
# To customize crossover distribution index (equivalent to mco's cdist=5):
crossover_fn <- function(object, parents) {
  nsgareal_sbxCrossover(object, parents, nc = 5)
}
nsga2(..., crossover = crossover_fn)
```

For most use cases, rmoo's defaults (`nc=20`, `nm=0.2`) provide good performance.

---

## Test File Updates

All test files using `nsga_control` need parameter name updates:

```r
# BEFORE
nsga_control = list(popsize = 20, generations = 20)

# AFTER
nsga_control = list(popSize = 20, maxiter = 20)
```

Assertions checking stored parameters:
```r
# BEFORE
expect_equal(res@control$nsga2$popsize, 128)
expect_equal(res@control$nsga2$generations, 150)

# AFTER
expect_equal(res@control$nsga2$popSize, 128)
expect_equal(res@control$nsga2$maxiter, 150)
```

---

## DESCRIPTION Update

```diff
Imports:
    methods,
    stats,
    utils,
    Matrix,
    SummarizedExperiment,
    S4Vectors,
-   mco,
+   rmoo,
    limma,
    glmnet,
    pROC
```

---

## Optional Enhancements (Phase 3)

### 1. Algorithm Selection Parameter

rmoo supports NSGA-II, NSGA-III, and R-NSGA-II. Consider adding an `algorithm` parameter:

```r
optimize_panel <- function(..., algorithm = c("NSGA-II", "NSGA-III")) {
  algorithm <- match.arg(algorithm)

  if (algorithm == "NSGA-II") {
    result <- rmoo::nsga2(...)
  } else if (algorithm == "NSGA-III") {
    # NSGA-III requires reference points
    result <- rmoo::nsga3(..., n_partitions = 12)
  }
}
```

### 2. Parallel Processing

rmoo supports parallel fitness evaluation (via the GA package infrastructure):
```r
# Add parallel = TRUE for parallel fitness evaluation
# Requires doParallel or similar parallel backend
rmoo::nsga2(..., parallel = TRUE)
```

### 3. Visualization Methods

rmoo provides built-in plotting:
```r
plot(result)                    # Scatter plot (auto-selects 2D/3D/pairwise)
plot(result, type = "pcp")      # Parallel coordinates
plot(result, type = "heatmap")  # Heatmap
plot(result, type = "polar")    # Polar coordinates
```

---

## Verification Plan

1. **Unit Tests**: Run existing test suite after migration
   ```bash
   R -e "devtools::load_all(); testthat::test_file('tests/testthat/test-optimization.R')"
   ```

2. **Comparison Test**: Run same optimization with both packages and compare Pareto fronts:
   - Same fitness function should produce comparable (not identical due to stochastic nature) Pareto fronts
   - Verify objective values are in same range

3. **Integration Test**: Run full `optimize_panel()` workflow on test data
   ```bash
   R -e "devtools::load_all(); testthat::test_file('tests/testthat/test-optimization.R')"
   ```

4. **Package Check**: Ensure clean build
   ```bash
   R CMD check .
   ```

---

## Migration Checklist

- [ ] Update `DESCRIPTION` (mco → rmoo in Imports)
- [ ] Update `optimize_panel()` roxygen docs (lines 3, 47-49: mco → rmoo references)
- [ ] Update `optimize_panel()` requireNamespace check (lines 98-101: mco → rmoo)
- [ ] Update `optimize_panel()` parameter construction (lines 370-379)
- [ ] Update `optimize_panel()` function call (line 389: mco::nsga2 → rmoo::nsga2)
- [ ] Update `optimize_panel()` result extraction (lines 391-398: add front filtering)
- [ ] Update `.get_adaptive_nsga_defaults()` parameter names
- [ ] Update test files with new parameter names
- [ ] Run `devtools::document()` to regenerate Rd files
- [ ] Run test suite
- [ ] Run `R CMD check`
- [ ] (Optional) Add algorithm selection parameter
- [ ] (Optional) Add parallel processing support
- [ ] (Optional) Add visualization wrapper functions
