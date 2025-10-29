# biomarkerPanels <img src="src/moo_hexsticker.png" alt="biomarkerPanels hex sticker" align="right" width="150"/>

Multi-objective optimization for discovering biomarker panels
with constrained performance targets. This package is designed for integration
with Bioconductor ecosystems and supports high-dimensional datasets.

## Key Features
- Multi-objective search balancing sensitivity, specificity, and custom costs.
- Support for gene-expression inputs via `SummarizedExperiment`.
- Panel evaluation utilities including ROC curves and calibration diagnostics.
- Transferability diagnostics baked into the loss registry (minimum cohort
  sensitivity/specificity, cohort gaps, calibration, and distribution-shift
  penalties).
- Hard constraints on optimisation via metric thresholds (e.g., enforce minimum
  sensitivity before maximising specificity).
- Extensible architecture compatible with Rcpp acceleration.

## Repository Guide
- `R/`: package source modules (simulation, objectives, optimization engines).
- `data-raw/`: scripts used to prepare package datasets and fixtures.
- `inst/extdata/`: example datasets bundled with the package.
- `tests/`: testthat suite and reusable fixtures.
- `vignettes/`: R Markdown tutorials for end-to-end workflows.
- `src/`: C++ code for certain components.

## Getting Started
Development installation will use:

```r
# install.packages("remotes")
remotes::install_github("SydneyBioX/moo")
```