# biomarkerPanels News

## 0.0.0.9000
- Initialize package skeleton with module stubs and documentation scaffolding.
- Provide simulated gene-expression fixture for testing and demonstration.
- Preserve feature identifiers within `select_transferable_features()` outputs,
  ensuring downstream consumers (e.g., optimizer coefficient matrices) handle
  duplicate assay names safely.
- Update `optimize_panel()` to respect feature pools when using the pairwise
  cohort aggregator, mapping aggregated identifiers back to their constituent
  features to avoid missing-feature errors.
- Enhance `evaluate_panel()` with confusion-matrix summaries and ROC diagnostics
  (including optional `ggplot2` visualisations and `evalm`-ready data) so
  validation runs expose the operating point used for reported metrics.
