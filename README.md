# biomarkerPanels

![biomarkerPanels hex sticker](src/moo_hexsticker.png)

Multi-objective optimization workflows for discovering small biomarker panels
with constrained performance targets. The package is designed for integration
with Bioconductor ecosystems and supports datasets where the number of features
greatly exceeds the sample size.

## Key Features (planned)
- Multi-objective search balancing sensitivity, specificity, and custom costs.
- Support for gene-expression inputs via `SummarizedExperiment`.
- Panel evaluation utilities including ROC curves and calibration diagnostics.
- Extensible architecture compatible with Rcpp acceleration.

## Repository Guide
- `R/`: package source modules (simulation, objectives, optimization engines).
- `data-raw/`: scripts used to prepare package datasets and fixtures.
- `inst/extdata/`: example datasets bundled with the package.
- `tests/`: testthat suite and reusable fixtures.
- `vignettes/`: R Markdown tutorials for end-to-end workflows.
- `src/`: optional C++ code for high-performance components.

## Getting Started
Development installation will use:

```r
# install.packages("BiocManager")
BiocManager::install("your-org/biomarkerPanels")
```

## Contributing
We welcome contributions. Please open an issue to discuss new objectives,
optimizers, or dataset support before submitting pull requests.
