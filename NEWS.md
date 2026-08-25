# biomarkerPanels 0.2.0.9000

## Bug fixes

* A `data.frame` passed as `x` to `evaluate_panel()`, `evaluate_pareto_front()`
  or the cohort helpers was mistaken for a multi-cohort list (`is.list()` is
  `TRUE` for data.frames), failing with "When `x` is a list, `y` must also be
  a list." All cohort-list dispatch now goes through `.is_cohort_list()`, which
  excludes data.frames and SummarizedExperiments; `optimize_panel_transferable()`
  likewise rejects a bare data.frame instead of splitting its columns into
  cohorts.

# biomarkerPanels 0.2.0

Executes the 11 resolved cleanup decisions
(`dev/cleanup-decisions-2026-08-RESOLVED.md`). Unlike 0.1.0, several of these
are deliberate behaviour changes.

## Breaking changes

* CPOP subsystem retired: `select_cpop_features()`, `fit_cpop_panel()` and
  `colmeans_penalty()` are gone. If a CPOP comparison is revisited, use the
  standalone CPOP package.
* Deprecated shims removed: `get_top_de_features()` (use `select_de_features()`,
  same arguments) and `select_denominator_features()` (use
  `select_batch_associated_features()`; note `n_denominators` is now
  `n_features` and the result field `$denominators` is now `$features`).
* `optimize_panel(scoring_fn = )` removed — it was only honoured with
  `fitness_cv = FALSE` and `regularized = FALSE`. `evaluate_panel(scoring_fn = )`
  is unchanged and remains the supported no-model scoring path.
* Fitness cache knobs `cache_fitness` and `cache_max_entries` removed from
  `optimize_panel()` and `optimize_panel_transferable()`; the fitness cache is
  now always on, with no eviction (also removes an O(n^2) key-tracking cost).
* `select_transferable_features(cv_control = )` removed; the internal
  `cv.glmnet` call now uses a fixed argument list.
* `TransferablePanelResult` loses its always-empty `validation_metrics` slot.
  **Saved `TransferablePanelResult` RDS objects from older versions may fail
  class validity checks on load** — re-create them with `calibrate_panel()`.
* `fit_panel()`'s inner cross-validation now caps folds at 5 (was 10), matching
  the NSGA scoring paths. `lambda.min` can land elsewhere, so refit final-model
  coefficients and `panel_metrics()` may shift slightly; NSGA optimization
  output is unchanged. Re-run downstream analyses that report fitted-model
  numbers.

## Improvements

* All three binomial `cv.glmnet` sites now share one internal helper
  (`.fit_binomial_glmnet()`), so the fold rule cannot drift again.
* `fitness_mode = "loco"` now has dedicated test coverage mirroring the
  rotating-validation tests.
* Internal helpers now use `@noRd`: `man/` drops from 212 to 74 pages, leaving
  only exported topics. Their roxygen prose remains in the source files.
* `dev/test_parallel_timing.R` retired; its timing figures predated the shared
  fitness scaffold and no longer described the current code.
* The getting-started vignette gains the live rule-out recipe
  (`specificity_at_sensitivity` + `pauc` + `min_cohort_auc` + `num_features`)
  as a worked example.
* Documented that `require_sign_consistency = FALSE` is equivalent to
  `sign_consistency_threshold = 0` in `select_transferable_features()`.

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
