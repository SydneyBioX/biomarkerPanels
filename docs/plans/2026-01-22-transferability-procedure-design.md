# Transferability Procedure Design

## Overview

Add a transferability-aware optimization procedure via `optimize_panel_transferable()`. This wrapper implements a train/validation/held-out split strategy with Neyman-Pearson threshold selection (via the `nproc` package) to maximize cross-cohort generalizability.

## Requirements

1. Partition each cohort into train / validation / held-out splits (stratified by class)
2. Perform feature selection on training data only (no data leakage)
3. Train NSGA optimization on concatenated training data
4. Evaluate MOO objectives on validation data (pooled + per-cohort metrics)
5. Select classification threshold on held-out data using nproc's NP classification
6. Return extended result with NP threshold and per-cohort metrics

## Dependencies

- **nproc**: Neyman-Pearson classification for threshold selection
  - CRAN: https://cran.r-project.org/web/packages/nproc/index.html
  - Add to Suggests in DESCRIPTION

## API Design

```r
optimize_panel_transferable(
  x, y,                                # List of matrices and responses
  objectives = define_objectives(...),
  max_features = 5L,
  feature_pool = NULL,                 # If NULL, feature selection done internally on training data
  cohort_aggregator = "pairwise_ratios",
  feature_alignment = "intersection",

  # Partition ratios (must sum to 1.0)
  train_ratio = 0.7,
  val_ratio = 0.2,
  # heldout_ratio derived: 1 - train_ratio - val_ratio

  # Neyman-Pearson threshold selection
  np_alpha = 0.15,                     # Type I error bound (1 - desired specificity floor)
  np_delta = 0.05,                     # Tolerance for type I error violation probability

  constraints = list(),
  algorithm = "NSGA-III",
  nsga_control = list(),
  seed = NULL,
  regularized = TRUE,
  regularized_alpha = 0.5,
  ...
)
```

Returns: `TransferablePanelResult` (extends `BiomarkerPanelResult`)

### Parameter Validation

- `train_ratio + val_ratio` must be < 1.0 (remainder becomes held-out)
- `train_ratio` must be >= 0.5 (need sufficient training data)
- `val_ratio` must be >= 0.1 (need sufficient validation data)
- Derived `heldout_ratio = 1 - train_ratio - val_ratio` must be >= 0.05

## Data Flow

### Step 1: Stratified Partitioning

For each cohort independently:
- Split "Yes" class indices by ratio (train/val/heldout)
- Split "No" class indices by ratio (train/val/heldout)
- Combine to get train_i, val_i, heldout_i (each maintains class balance)

Warn if any partition has <5 samples per class.

### Step 2: Feature Selection (Training Data Only)

**Critical: Feature selection must use only training partition to prevent data leakage.**

If `feature_pool = NULL`:
1. Run `select_transferable_features()` or `get_top_de_features()` on `train_combined` only
2. Use resulting features as the candidate pool for NSGA

If `feature_pool` is provided:
- Assume user has handled leakage concerns externally
- Use provided features directly

### Step 3: Concatenate by Role

```
train_combined   = rbind(train_1, train_2, ..., train_k)
val_combined     = rbind(val_1, val_2, ..., val_k)
heldout_combined = rbind(heldout_1, heldout_2, ..., heldout_k)
```

Track cohort membership via factor vectors: `train_cohort`, `val_cohort`, `heldout_cohort`.

### Step 4: Validation-Based Fitness Evaluation

For each NSGA candidate:
1. Select features (weight > 0.5, up to max_features)
2. Train model on `train_combined`
3. Predict on `val_combined`
4. Compute primary objectives: pooled metrics on concatenated validation
5. Compute secondary objectives: cohort gap (max - min metric)
6. Return objective values

No internal CV during fitness - the explicit train/val split handles generalization.

### Step 5: Final Model Fitting

After NSGA selects optimal feature set:
1. Combine train + val data (forming ~90% of original)
2. Fit final regularized model on combined set
3. Generate probability predictions using `type = "response"` (scores in [0,1])

### Step 6: Neyman-Pearson Threshold Selection

Use the `nproc` package for principled threshold selection:

```r
# Generate predictions on held-out (pooled across cohorts)
heldout_scores <- predict(final_model, heldout_combined, type = "response")

# Use nproc to find NP-optimal threshold
# - Bounds type I error (FPR) at alpha with probability >= 1 - delta
# - Minimizes type II error (FNR) subject to that bound
np_fit <- nproc::npc(
  x = NULL,
  y = heldout_y,
  method = "custom",
  score = heldout_scores,
  alpha = np_alpha,
  delta = np_delta,
  split = 0
)

np_threshold <- np_fit$cutoff
```

### Step 7: Per-Cohort Metrics and Weighted Variance

Compute metrics at the NP threshold for each cohort:

```r
per_cohort_metrics <- lapply(unique(heldout_cohort), function(c) {
  idx <- heldout_cohort == c
  pred_class <- ifelse(heldout_scores[idx] >= np_threshold, "Yes", "No")
  truth <- heldout_y[idx]
  list(
    cohort = c,
    n = sum(idx),
    sensitivity = mean(pred_class[truth == "Yes"] == "Yes"),
    specificity = mean(pred_class[truth == "No"] == "No")
  )
})

# Weighted variance using inverse sample size
weights <- 1 / sapply(per_cohort_metrics, `[[`, "n")
weights <- weights / sum(weights)

sens_values <- sapply(per_cohort_metrics, `[[`, "sensitivity")
spec_values <- sapply(per_cohort_metrics, `[[`, "specificity")

weighted_var <- list(

sensitivity = sum(weights * (sens_values - weighted.mean(sens_values, weights))^2),
  specificity = sum(weights * (spec_values - weighted.mean(spec_values, weights))^2)
)
```

## Result Class

```r
setClass(
  "TransferablePanelResult",
  contains = "BiomarkerPanelResult",
  slots = c(
    np_threshold = "numeric",
    np_alpha = "numeric",
    np_delta = "numeric",
    per_cohort_metrics = "data.frame",
    weighted_variance = "list",        # list(sensitivity = x, specificity = y)
    validation_metrics = "list",
    partition_info = "list"            # Stores indices for reproducibility
  )
)
```

Backward compatible: `panel_metrics()`, `panel_model()`, `evaluate_panel()` work unchanged.

New accessors:
- `np_threshold(result)` - Returns the NP classification threshold
- `per_cohort_metrics(result)` - Returns data.frame of per-cohort performance
- `weighted_variance(result)` - Returns list of weighted variance by metric

## File Structure

### New Files

| File | Purpose |
|------|---------|
| `R/optimize_panel_transferable.R` | Main wrapper + validation fitness + NP threshold |
| `tests/testthat/test-transferability.R` | Unit tests |

### Modified Files

| File | Changes |
|------|---------|
| `R/cohort_preparation.R` | Add `.stratified_partition_cohorts()` |
| `R/panel_class.R` | Add `TransferablePanelResult` class + accessors |
| `DESCRIPTION` | Add nproc to Suggests |
| `NAMESPACE` | Export new function and class |

## Implementation Sequence

1. Add nproc to DESCRIPTION Suggests
2. Add `.stratified_partition_cohorts()` to `R/cohort_preparation.R`
3. Add `TransferablePanelResult` class to `R/panel_class.R`
4. Create `R/optimize_panel_transferable.R` with:
   - `.validate_partition_ratios()` helper
   - `.make_validation_fitness()` factory
   - `.select_np_threshold()` wrapper around nproc::npc
   - `.compute_per_cohort_metrics()` helper
   - Main `optimize_panel_transferable()` wrapper
5. Add tests in `tests/testthat/test-transferability.R`
6. Run `devtools::document()`

## Edge Cases

1. **Small cohorts**: Warn if any partition has <5 samples per class; error if <2
2. **nproc failure**: If nproc::npc fails, fall back to threshold = 0.5 with warning
3. **Single cohort**: Weighted variance is 0; NP threshold still valid
4. **Invalid ratios**: Error if ratios don't sum to 1.0 or violate minimums
5. **No nproc installed**: Error with informative message to install from CRAN

## Test Cases

1. Partitioning preserves class balance within tolerance
2. Feature selection uses only training data (data leakage test)
3. Basic wrapper returns TransferablePanelResult
4. NP threshold is within [0,1] (probability scores)
5. Per-cohort metrics populated for all cohorts
6. Weighted variance computed correctly
7. evaluate_panel() works with TransferablePanelResult
8. Invalid partition ratios rejected with clear error
9. Ratios not summing to 1.0 rejected
10. Single cohort handled gracefully
11. Warning emitted for small partitions
12. nproc integration works with custom scores

## Trade-offs

**What nproc provides:**
- Principled, statistically grounded threshold selection
- Guarantees on type I error control with probability >= 1 - delta
- No arbitrary grid search parameters

**Limitations:**
- Does not directly optimize for cross-cohort variance
- Type I error control (specificity) rather than type II (sensitivity)

**Mitigation:** Per-cohort metrics and weighted variance are computed and stored for transparency. Users can assess transferability from variance values and adjust `np_alpha` if needed.
