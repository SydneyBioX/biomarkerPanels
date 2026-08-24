# Codebase cleanup 2026-08-23 — what was done, and the decisions left open

Branch: `worktree-codebase-cleanup` (5 commits on top of `22a975e`).
Produced by a four-agent survey (feature selection, optimization/fitness,
metrics/evaluation, repo-wide dead-code audit) followed by verified execution.
Every deletion was re-verified by grep across this repo,
`biomarkerPanels-analysis`, AND `summer-students-26/Summer-Biomarker-2026`
(the teaching repo is a second live consumer — 16 files load this package).

## Executed (behavior-preserving unless noted)

1. **Deleted dead exports**: `select_panel_top_sensitivity()`,
   `select_panel_by_pathway()` (zero callers anywhere),
   `metric_easy_hard_accuracy` (whole `metric_difficulty.R`; its own @note said
   it is incompatible with `optimize_panel_transferable()`),
   `metric_conditional_score_shift` + `.w1` (test-only), and the
   `.predict_from_model` alias. NOTE: `select_panel_inclusion_frequency()` was
   KEPT — one survey agent claimed it was dead, but it has 2 live call sites in
   `biomarkerPanels-analysis/0_other/Chloe-*.qmd`.
2. **Un-exported 15 metric functions** (kept as internals + docs): downstream
   only ever reaches metrics via registry strings in `define_objectives()`.
   `metric_specificity_at_sensitivity` stays exported (direct call in
   `5_fetalgrowth/dev_ruleout_panel_phase1.R`). Tests alias internals via
   `tests/testthat/helper-metrics.R`.
3. **Deduplication**: identical SE/matrix branches in `.prepare_cohort_inputs`;
   `.resolve_panel_solution()` shared by `fit_panel`/`fit_np_panel` (~110-line
   verbatim copy); `.build_roc()` (6 hand-rolled `pROC::roc` sites);
   `.confusion_counts()` (5 classification metrics); `.make_fitness_scaffold()`
   (4x fitness skeleton); `.run_nsga()`, `.validate_selection_threshold()`,
   `.min_features_required()`, `.combine_partitions()` (2-3x blocks in the two
   entry points); RUV `.pooled_pvalue_controls()` shared by initial + iterated
   control selection.
4. **Bug fixes** (small behavior changes, all in the safe direction):
   - RUV: a full RUV-4 fit + variance adjustment was computed and immediately
     discarded when `iterate_controls = TRUE` (pure waste; results unchanged).
   - `fit_np_panel`/`.compute_per_cohort_metrics` now return `NA` (not `NaN`)
     for empty classes, via the metric functions.
   - `evaluate_panel_by_cohort(metrics=...)` now accepts ANY registered metric
     (was a hard-coded sensitivity/specificity/auc switch).
   - `control$scoring_function` no longer stores the deparsed body of the
     default scorer.
   - LOCO single-class folds now warn instead of silently dropping a cohort
     from the metrics (mirrors the rotating-validation guard).
5. **File moves/merges**: `de_features_helpers.R` -> `de_features.R`;
   `feature_stability.R` -> `stable_features.R` and
   `feature_stability_analysis.R` -> `pareto_feature_stability.R`;
   `fitness_cache_utils.R` split into `panel_selection.R` + `glm_fitting.R` +
   cache-only remainder; `utils.R` dissolved (`pairwise_col_diff` ->
   `feature_transforms.R`, `.filter_dominated` -> `nsga_results.R`);
   `.as_cohort_list` -> `cohort_inputs.R`; `.predict_np_model` ->
   `model_prediction.R`; `metric_cohort_transferability.R` merged into
   `metric_cohort.R`.
6. **Hygiene**: phantom `define_ruleout_objectives()` removed from the
   `fit_np_panel` example (it never existed; see decision 6 below); stale
   `doc/*.html` deleted; `inst/benchmarks` and `.claude` added to
   `.Rbuildignore` (benchmarks were shipping untracked AND call a function
   that no longer exists); vignette `scoring_fn` section corrected (see
   decision 5).

## Decisions still open (deliberately NOT done)

1. **CPOP subsystem (~700 lines + 137 test lines)** serves ONE external call
   site (`select_cpop_features`); `fit_cpop_panel` and `colmeans_penalty` have
   zero callers anywhere. Options: keep as-is (published algorithm, standalone
   value); un-export `fit_cpop_panel`/`colmeans_penalty`; or retire the module.
2. **`fitness_mode = "loco"` and `"within_cohort_rotating"`** (~665 lines
   combined with their optimize_panel_transferable branches). Research
   leftovers from the FGR sweep; `loco` has NO package test and one notebook
   caller. This is the CV_STRATEGY_BRIEF.md methods question — retire or add a
   test, but decide there. The ~400-line fold-plan unification of all four
   fitness strategies should wait until this is settled.
3. **`scoring_fn` is silently ignored under default arguments.** It only takes
   effect when `fitness_cv = FALSE` AND `regularized = FALSE`. Either thread it
   through `.compute_cv_scores()` or deprecate it. (The vignette now documents
   the restriction honestly.)
4. **Inner-CV budget mismatch**: NSGA candidate scoring uses
   `nfolds = min(5, ...)` (`cv_scoring.R`) while the final `fit_panel()` model
   uses `min(10, ...)` (`model_fitting.R`). Probably unintentional; unifying
   changes optimization numbers, so it needs a deliberate call (and re-running
   any affected comparisons). Related: the three `cv.glmnet` sites could share
   a `.fit_binomial_glmnet()` helper once the nfolds question is settled.
5. **`cache_max_entries`** eviction is unreachable (never set finite anywhere;
   `cache_fitness = FALSE` is test-only). Removing both knobs simplifies five
   signatures (~30 lines) but is an API change.
6. **`define_ruleout_objectives()`**: implement (~20-line wrapper the analysis
   repo hand-rolls in at least two places) or leave removed. CLAUDE.md still
   advertises it — see CLAUDE.md edits below.
7. **Deprecated shims**: `select_denominator_features` must STAY (12 downstream
   call sites; new name has zero adoption — migrating is a downstream task).
   `get_top_de_features` stays for now (7 call sites in summer-students-26 +
   `dev/test_parallel_timing.R`); migrate those, then drop it.
8. **Never-read S4 slots** (`training_signature`; `validation_metrics` /
   `partition_info` always `list()` on TransferablePanelResult): wiring up or
   dropping changes the class definition and can invalidate previously saved
   RDS objects — left untouched deliberately.
9. **`cv_control` / `require_sign_consistency`**: provably unused / redundant
   (`require_sign_consistency` == `sign_consistency_threshold = 0`); candidates
   for a deprecation cycle at the next breaking release.
10. **`@keywords internal` Rd files**: 117 of ~200 man pages document internal
    dot-helpers. Switching those roxygen blocks to `@noRd` would shrink
    `man/` by half; policy call.
11. **`dev/test_parallel_timing.R`** still calls `get_top_de_features()`
    (works via the shim) and is the cited source of CLAUDE.md's parallel
    benchmark numbers — update it when convenient (dev/ is untracked, so it
    was not edited from the cleanup branch).

## CLAUDE.md edits needed (file is untracked; could not be edited from the branch)

- "Pareto Solution Selection" bullet: remove `select_panel_top_sensitivity()`
  and `select_panel_by_pathway()`; keep `select_panel_inclusion_frequency()`.
- "Rule-Out Diagnostics" + "Main Entry Points": remove/replace
  `define_ruleout_objectives()` (never existed). The real recipe:
  `define_objectives(metrics = c("pauc", "num_features"), params = list(pauc =
  list(sens_floor = 0.9)))` + `constraints =
  list(min_metric_constraint("sensitivity", 0.9))`.
- Metric registry section: note that `metric_*` functions are internal; users
  reach them via registry names (`define_objectives(metrics = "auc")`).
- Architecture section: `de_features_helpers.R` gone;
  `stable_features.R`/`pareto_feature_stability.R` renames;
  `panel_selection.R`/`glm_fitting.R` split out of `fitness_cache_utils.R`;
  `utils.R` gone.
