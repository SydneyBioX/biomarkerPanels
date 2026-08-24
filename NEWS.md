# biomarkerPanels 0.1.0

First versioned milestone: a behaviour-preserving cleanup of the package
surface and internals.

## Breaking changes

* 19 exports removed. Dead exports (zero callers) were deleted outright —
  `select_panel_top_sensitivity()`, `select_panel_by_pathway()`,
  `metric_conditional_score_shift()` and `metric_easy_hard_accuracy()` — and
  15 registry-only `metric_*` functions are now internal: reach them through their registry names in
  `define_objectives()`. `metric_specificity_at_sensitivity()` remains
  exported.

## Improvements

* `metric_specificity_at_sensitivity()` now returns 0 (the worst achievable
  value) on degenerate inputs instead of erroring or returning `NA`, matching
  `metric_sensitivity_at_specificity()`, so it is safe as an NSGA objective.
* Feature-selection front-ends share a common input-preparation path;
  metrics/evaluation internals deduplicated; fitness evaluation consolidated
  into a shared scaffold and NSGA launch path.
* File reorganisation for one-theme-per-file layout (see
  `dev/cleanup-decisions-2026-08.md` for the survey behind this release).
