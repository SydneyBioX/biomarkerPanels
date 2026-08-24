# Developer Guide: biomarkerPanels Architecture

## High-Level Overview

The `biomarkerPanels` package is designed to derive small, interpretable biomarker panels from complex gene expression data (or other omics datasets). The architecture is built around a multi-objective optimization pipeline that balances competing criteria like sensitivity, specificity, assay cost (panel size), and transferability across multiple cohorts. 

The package is organized into several distinct modules that handle the data lifecycle from raw expression matrices to fully evaluated, deployable panels.

```text
                      +-----------------------------+
                      |   1. Cohort Data Inputs     |
                      | (cohort_inputs.R, data_...) |
                      +-----------------------------+
                                     |
                                     v
                      +-----------------------------+
                      |  2. Feature Pre-selection   |
                      | (discriminative, RUV, etc.) |
                      +-----------------------------+
                                     |
                                     v
+------------------+  +-----------------------------+  +--------------------+
|  Metrics & Objs  |->|  3. Core Optimization       |<-| Feature Transforms |
| (objectives.R)   |  | (NSGA-II / III)             |  | (pairwise, norm)   |
+------------------+  | (optimize_panel.R)          |  +--------------------+
                      +-----------------------------+
                                     |
                                     v
                      +-----------------------------+
                      | 4. Model Fit & Calibration  |
                      | (model_fitting.R, np_fit.R) |
                      +-----------------------------+
                                     |
                                     v
+------------------+  +-----------------------------+  +--------------------+
| Summarization    |<-| 5. Evaluation & Prediction  |->|   Visualization    |
| (panel_summ...)  |  | (evaluate_panel.R, etc.)    |  | (visualization.R)  |
+------------------+  +-----------------------------+  +--------------------+
```

## Module Descriptions

### 1. Inputs & Data Preparation
- **Files:** `cohort_inputs.R`, `feature_alignment.R`, `data_partitioning.R`
- **Purpose:** Harmonizes incoming expression data (which may be provided as lists of matrices, data.frames, or `SummarizedExperiment` objects) across multiple cohorts. It aligns features, resolves differences in variable availability (e.g., intersecting features), and provides safe partitioning for training/held-out splits.

### 2. Feature Pre-selection
- **Files:** `discriminative_features.R`, `ruv_features.R`, `transferable_features.R`, `denominator_features.R`
- **Purpose:** Large omics datasets are computationally expensive to optimize directly. Pre-selection filters the feature space using various statistical methods (e.g., Wilcoxon, RUV, ridge transferability) to form a highly-qualified "feature pool" that is small enough for NSGA-II/III to efficiently search.

### 3. Core Optimization (NSGA-II / NSGA-III)
- **Files:** `optimize_panel.R`, `optimize_panel_transferable.R`, `objectives.R`, `metric_registry.R`, `feature_transforms.R`, `fitness_cache_utils.R`, `cv_scoring.R`
- **Purpose:** The central engine of the package. It uses the `rmoo` package's NSGA implementations to search for pareto-optimal biomarker subsets. 
- **Dependencies:** 
  - **Metrics & Objectives:** Users define what to optimize (e.g., maximize sensitivity, minimize panel size) via `define_objectives()`.
  - **Feature Transforms:** Raw values can be transformed on-the-fly (e.g., `pairwise_ratios` or `reference_norm`) before evaluation to reduce batch effects.
  - **Fitness Caching:** Fast lookup of previously scored panels ensures optimization runs efficiently.

### 4. Model Fitting & Calibration
- **Files:** `model_fitting.R`, `np_fit.R`, `threshold_calibration.R`
- **Purpose:** Given a specific panel of selected features, this module fits a diagnostic model. It supports standard Logistic Regression/Regularized models as well as Neyman-Pearson (`np_fit.R`) bounds for rigorous control of Type-I errors (e.g., ensuring high specificity constraints).

### 5. Evaluation & Prediction
- **Files:** `evaluate_panel.R`, `evaluate_cohorts.R`, `evaluate_pareto.R`, `evaluate_sensitivity.R`
- **Purpose:** Handles out-of-sample validation on held-out test data. It evaluates single panels or the entire Pareto front and extracts comprehensive operational metrics (AUC, sensitivity, specificity, accuracy).

### 6. Summarization & Visualization
- **Files:** `panel_summarization.R`, `panel_diversity.R`, `visualization.R`, `hypervolume.R`
- **Purpose:** Translates complex optimization outputs into actionable insights. `panel_diversity.R` measures feature sharing across the Pareto front, while `visualization.R` provides high-level plotting tools to inspect pareto fronts, feature importance, and performance trade-offs.

## Unit Testing Strategy & Refactoring Plan

The `biomarkerPanels` package maintains a robust unit testing suite. However, to keep tests maintainable, performant, and easy to navigate, we utilize smaller, categorically split scripts paired with centralized helper functions.

### Testing Architecture

1. **Categorical Splitting**: Test files are named using the `test-{module}-{feature}.R` convention.
   - For example, `test-optimization.R` is split into `test-optimize-core.R`, `test-optimize-aggregators.R`, `test-optimize-cohorts.R`, and `test-optimize-algorithms.R`.
   - `test-metric-functions.R` is split into `test-metrics-base.R`, `test-metrics-cohort.R`, and `test-metrics-objectives.R`.
2. **Centralized Helpers (DRY)**: Repeated data setup and mock generation are centralized in `tests/testthat/helper-data.R` (and other `helper-*.R` files). These are automatically loaded by `testthat` for all tests.
3. **Decoupled Integration Tests**: Tests validating S4 data structures or simple threshold logic use mock objects instead of running the full stochastic optimization loops (e.g., `nsga3`).

### Execution Workflow for Refactoring

When refactoring or breaking up long test scripts, the following workflow must be strictly followed to ensure stability:

1. **One-by-One Execution**: We will work with each test script one by one. We do not attempt to split or rewrite all files simultaneously.
2. **Immediate Validation**: At the end of editing or splitting a single script, we immediately revise and review that the solution actually works by running the `devtools::test()` function (or `devtools::test(filter='pattern')`).
3. **Immediate Remediation**: If there are any problems, errors, or test failures with the script that was just edited, we will fix those errors immediately before moving on to the next test script.
