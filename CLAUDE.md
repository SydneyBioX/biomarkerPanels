# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Install package dependencies
R -e "devtools::install_deps()"

# Build and install package (includes C++ compilation)
R CMD INSTALL .

# Run all tests
R -e "devtools::test()"

# Run a single test file
R -e "testthat::test_file('tests/testthat/test-optimization.R')"

# Check package (CRAN-style)
R CMD check .

# Generate documentation from roxygen2 comments
R -e "devtools::document()"
```

## Architecture

### Core Pipeline

This R package implements multi-objective optimization (NSGA-II via `mco`) for discovering compact biomarker panels with constrained performance targets.

**Data Flow:**
1. **Feature Selection** (`R/feature_selection.R`): Filters candidate features via limma DE analysis or ridge regression transferability scoring
2. **Cohort Aggregation** (`R/utils.R`): Applies pairwise difference transformation (`pairwise_col_diff`) to dampen batch effects across cohorts
3. **Optimization** (`R/optimization.R`): `optimize_panel()` runs NSGA-II search over feature weights, evaluating candidates via registered loss functions
4. **Evaluation** (`R/evaluation.R`): `evaluate_panel()` computes validation metrics with ROC diagnostics

### Key Abstractions

**Loss Registry** (`R/loss_functions.R`):
- All objectives (sensitivity, specificity, c_index, cohort gaps, etc.) are registered in `.loss_registry` environment
- Each loss has a standard signature: `function(truth, scores, selected, cohort, x, ...)`
- Use `register_loss_function()` to add custom objectives

**BiomarkerPanelResult** (`R/panel_class.R`):
- S4 class storing selected features, metrics, Pareto-optimal solutions, and optimization settings
- Access metrics via `panel_metrics(result)`

**Cohort Aggregator**:
- Default `"pairwise_ratios"` generates all pairwise column differences via Rcpp (`src/pairwise.cpp`)
- This transformation makes features more transferable across cohorts with different baselines

### Multi-Cohort Support

When `x` is a list of matrices (one per cohort):
- Features are aligned via simple intersection of column names
- Truth/cohort vectors are concatenated for joint optimization
- Cohort-aware losses (e.g., `loss_min_cohort_sensitivity`) evaluate worst-case performance

### Constraints

Use `min_metric_constraint()` to enforce hard thresholds (e.g., minimum sensitivity) during optimization. Infeasible candidates receive infinite objective values.

## Testing

Test fixtures are loaded from `tests/data/` via `fixture_path()` helper in `tests/testthat/helper-data.R`.
