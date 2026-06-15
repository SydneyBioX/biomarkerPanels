# biomarkerPanels <img src="man/figures/moo_hexsticker.png" alt="biomarkerPanels hex sticker" align="right" width="150"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/SydneyBioX/biomarkerPanels/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/SydneyBioX/biomarkerPanels/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[Documentation](https://sydneybiox.github.io/biomarkerPanels/) | [Issues](https://github.com/SydneyBioX/biomarkerPanels/issues) | [License](LICENSE)

Multi-objective optimization for discovering compact biomarker panels from high-dimensional data.

> **Note**
> Contributions are welcome! Please report any issues. You may also contribute by opening a pull request.

## Key Features

- **Multi-objective optimization**: NSGA-II/NSGA-III genetic algorithms balance sensitivity, specificity, AUC, panel size, and custom objectives simultaneously
- **Hard constraints**: Enforce minimum performance thresholds (e.g., sensitivity >= 0.90) during optimization
- **Multi-cohort support**: Built-in transferability objectives penalize cross-site performance gaps
- **Batch effect mitigation**: Pairwise feature ratios dampen distributional shifts across cohorts
- **Extensible metric registry**: Add custom objectives with `register_metric()`
- **Bioconductor integration**: Native support for `SummarizedExperiment` inputs
- **Pareto front selection**: Choose final panels by sensitivity, feature frequency, or pathway enrichment

## 🔨 Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("SydneyBioX/biomarkerPanels")
```

**Requirements**: R >= 4.4, plus dependencies (rmoo, glmnet, limma, pROC, SummarizedExperiment)

## 🚀 Quick Start

### Single cohort optimization

```r
library(biomarkerPanels)

# Define objectives to optimize
objectives <- define_objectives(
  metrics = c("sensitivity", "specificity", "num_features")
)

# Run NSGA-II optimization (returns Pareto front, no model)
opt_result <- optimize_panel(
  x = train_matrix,
  y = train_response,
  objectives = objectives,
  max_features = 10,
  seed = 42
)

# Inspect Pareto-optimal solutions
summarize_solutions(opt_result)
#>   solution_id n_features sensitivity specificity num_features
#> 1           1          4       0.912       0.847            4
#> 2           2          6       0.934       0.821            6
#> ...

# Fit model on chosen solution (or auto-select best)
panel <- fit_panel(opt_result, solution_id = 1)
# OR: auto-select best on first objective
panel <- fit_panel(opt_result)

# Evaluate on held-out data (requires fitted model)
eval <- evaluate_panel(panel, x = test_matrix, y = test_response)
eval$metrics
#>   sensitivity   specificity           auc
#>         0.912         0.847         0.923
```

### Multi-cohort optimization

For multi-site studies, pass data as named lists. Use cohort-aware objectives like `min_cohort_auc` to ensure the panel generalizes:

```r
result <- optimize_panel(
  x = list(site_A = mat1, site_B = mat2, site_C = mat3),
  y = list(site_A = y1, site_B = y2, site_C = y3),
  objectives = define_objectives(
    metrics = c("sensitivity", "min_cohort_auc", "num_features")
  ),
  max_features = 8
)
```

### Feature pre-filtering

Reduce the search space before optimization:

```r
# Via differential expression
top_de <- get_top_de_features(x, y, n = 50)

# Via cross-cohort transferability
transferable <- select_transferable_features(x_list, y_list, n = 50)

result <- optimize_panel(x, y, feature_pool = top_de, ...)
```

## 📋 Main Functions

| Function | Description |
|----------|-------------|
| `optimize_panel()` | Run NSGA-II/III, returns `OptimizationResult` with Pareto front |
| `summarize_solutions()` | Inspect Pareto solutions with metrics and feature counts |
| `fit_panel()` | Fit model on selected solution, returns `BiomarkerPanelResult` |
| `evaluate_panel()` | Validate panel performance on held-out data (requires fitted model) |
| `define_objectives()` | Configure optimization objectives |
| `min_metric_constraint()` | Add hard performance constraints |
| `select_panel_top_sensitivity()` | Select solution from Pareto front by sensitivity |
| `select_panel_inclusion_frequency()` | Select solution by feature frequency across solutions |
| `get_top_de_features()` | Pre-filter features via differential expression |
| `select_transferable_features()` | Pre-filter features by cross-cohort stability |
| `metric_registry()` | View all available objective functions |

## 📊 Available Objectives

| Objective | Description | Direction |
|-----------|-------------|-----------|
| `sensitivity` | True positive rate | maximize |
| `specificity` | True negative rate | maximize |
| `auc` | Area under ROC curve | maximize |
| `pauc` | Partial AUC (high-sensitivity region) | maximize |
| `num_features` | Panel size | minimize |
| `min_cohort_auc` | Worst-case AUC across cohorts | maximize |
| `cohort_auc_gap` | Max AUC difference between cohorts | minimize |
| `cohort_auc_var` | Variance of per-cohort AUC values | minimize |

See `metric_registry()` for the complete list.

## 🙋 FAQ

**Q: How do I optimize for rule-out screening (high sensitivity)?**

Use `define_ruleout_objectives()` which enforces a sensitivity constraint and optimizes partial AUC in the high-sensitivity region:

```r
objectives <- define_ruleout_objectives(min_sensitivity = 0.95)
```

**Q: Should I use NSGA-II or NSGA-III?**

NSGA-II (default) works well for 2-3 objectives. NSGA-III provides better diversity for many-objective problems (4+).

**Q: How do I add a custom objective?**

```r
register_metric(
  name = "my_metric",
  fun = function(truth, scores, selected, ...) { ... },
  direction = "maximize"
)
```

## License

GPL-3

## Issues

Please report bugs and feature requests via [GitHub Issues](https://github.com/SydneyBioX/biomarkerPanels/issues).

## How to Cite

If you use biomarkerPanels in your research, please cite:

```bibtex
@software{biomarkerPanels,
  author = {Robertson, Harry},
  title = {biomarkerPanels: Multi-objective Optimization for Biomarker Panel Discovery},
  url = {https://github.com/SydneyBioX/biomarkerPanels}
}
```
