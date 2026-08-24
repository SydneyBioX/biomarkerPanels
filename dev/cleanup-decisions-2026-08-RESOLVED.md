# Codebase cleanup 2026-08 — the 11 open decisions, resolved

Resolutions for the open items in `dev/cleanup-decisions-2026-08.md`, decided
2026-08-24. Nothing here is executed yet. Each entry records the choice, the
concrete edit list, and the consequence accepted.

Where the original doc's factual claims were wrong, the correction is noted —
several resolutions turn on those corrections.

---

## 1. CPOP subsystem — RETIRE ENTIRELY

Delete `R/cpop_select.R` (224), `R/cpop_fit.R` (176), `R/cpop_utils.R` (90),
`R/cpop_step2.R` (57), `tests/testthat/test-cpop.R` (137), the 6 `man/` pages,
2 NAMESPACE exports, 4 `Collate:` entries, and the `[select_cpop_features()]`
`@seealso` at `R/de_features.R:31`.

**Blocker:** `biomarkerPanels-analysis/5_fetalgrowth/2_FGR_ratios.Rmd:85` calls
`select_cpop_features()`. Migrate it to the standalone `/dskh/nobackup/harryr/CPOP`
repo, or pin it to an older package version, before this lands.

`fit_cpop_panel()` and `colmeans_penalty()` had zero callers outside tests.

## 2. `fitness_mode` loco / within_cohort_rotating — KEEP BOTH, ADD A LOCO TEST

Write a package test for `fitness_mode = "loco"` mirroring
`tests/testthat/test-rotating-validation.R`: runs on >=2 cohorts, errors on 1,
records `fitness_mode` in `result@control`.

**Doc correction:** the two are not equal-weight leftovers.
`within_cohort_rotating` has a test, 4 analysis call sites, and the
`n_val_splits` API; `loco` has no test and 2 call sites in
`2.1_FGR_hypervolume.Rmd`. Both are arms of the open
`dev/CV_STRATEGY_BRIEF.md` question, so neither is retired.

The ~400-line fold-plan unification stays deferred pending that brief.

## 3. `optimize_panel(scoring_fn=)` — REMOVE OUTRIGHT

Unreachable under defaults: only honoured when `fitness_cv = FALSE` AND
`regularized = FALSE`.

Delete from `R/optimize_panel.R`: signature `:115`, `@param` `:49`, validation
`:271-272`, `scoring_fn_label` / `control$scoring_function` `:262-267` + `:466`,
and collapse the custom-scoring branch `:379-390` so the non-CV path always uses
`.default_scoring_fn` / `.fit_predict_binomial_glm`. Remove the
`vignettes/getting-started.Rmd` section (~`:192`, `:220`, `:231`). Regenerate
`man/optimize_panel.Rd`.

**Do not touch `evaluate_panel(scoring_fn=)`** — it works, is tested
(`test-evaluate-single.R:32`), and is the documented no-model scoring path.

No downstream callers. API break.

## 4. Inner-CV nfolds — UNIFY TO THE 5-CAP

`R/model_fitting.R:198`: `min(10L, ...)` -> `min(5L, ...)`, matching
`cv_scoring.R:178` and `:304`. With all three sites then identical, extract
`.fit_binomial_glmnet(x, y, alpha)` wrapping the shared
`cv.glmnet(family = "binomial", type.measure = "deviance")` + nfolds rule.

**Consequence:** `fit_panel()` picks a different `lambda.min`, so final-model
numbers and `panel_metrics()` shift — re-run affected analysis results. NSGA
output is unchanged (those sites already used 5).

## 5. Cache knobs — REMOVE BOTH

Drop `cache_fitness` and `cache_max_entries` from `optimize_panel()`,
`optimize_panel_transferable()`, `.make_fitness_scaffold()`, `fitness_loco()`,
`fitness_rotating_validation()`. Delete `.validate_cache_controls()`; simplify
`.new_fitness_cache()` / `.cache_set()` to an always-on cache with no
`max_entries`, no eviction branch (`fitness_cache_utils.R:61-70`), and no
`cache$keys` vector — that vector's `c(cache$keys, key)` growth was O(n^2) over
an NSGA run.

**Consequence:** `test-fitness-cache-utils.R` exercises `cache_fitness` TRUE and
FALSE at 7 sites. The FALSE arm currently proves caching does not alter results
— replace it with a direct `.cache_get`/`.cache_set` unit test rather than
deleting the coverage.

Zero downstream callers of either knob. API break.

## 6. `define_ruleout_objectives()` — LEAVE REMOVED, FIX THE DOCS

- **CLAUDE.md** (untracked; edit in the working tree): rewrite "Rule-Out
  Diagnostics", remove the function from "Main Entry Points". Describe what
  exists: `metric_pauc`'s `sens_floor` and `specificity_at_sensitivity`'s
  `target_sensitivity`, composed via `define_objectives()`;
  `min_metric_constraint()` optional.
- **`vignettes/getting-started.Rmd`**: add the real live recipe as a worked
  example — `define_objectives(metrics = c("specificity_at_sensitivity",
  "pauc", "min_cohort_auc", "num_features"), params = list(
  specificity_at_sensitivity = list(target_sensitivity = 0.90),
  pauc = list(sens_floor = 0.90, partial_auc_correct = TRUE)))`, no constraint.

**Doc correction:** all three live call sites (`3.1_FGR_ruleout.Rmd`,
`dev_ruleout_panel_phase1.R`, `_v2.R`) use the 4-objective form with **no
constraint**, and the recipe changed between v1 and v2 (`spec_at_sens90` ->
`specificity_at_sensitivity`). Nothing stable enough to freeze into the API.

## 7. Deprecated shims — MIGRATE ALL 19 SITES, DROP BOTH

`R/deprecated.R` goes entirely, plus `man/biomarkerPanels-deprecated.Rd` and
2 NAMESPACE exports.

- `get_top_de_features` -> `select_de_features` (pure rename), **7 sites**:
  `summer-students-26/.../3_IBD/2_IBD_transferable_features.Rmd` (3),
  `5_FGR/1_FGR_features.Rmd` (3), `5_FGR/2_FGR_ratios.Rmd` (1).
  (The 8th, `dev/test_parallel_timing.R:28`, is moot — decision 11 deletes it.)
- `select_denominator_features` -> `select_batch_associated_features`,
  **12 sites** in `biomarkerPanels-analysis`:
  `3_IBD/4_IBD_feature_preselection.Rmd` (4),
  `5_fetalgrowth/15_exp2_dual_stream.Rmd` (2),
  `5_fetalgrowth/6_FGR_feature_preselection.Rmd` (2), and one each in
  `2_FGR_ratios.Rmd`, `5_FGR_final-model.Rmd`,
  `5_FGR_final-model_rotating.Rmd`, `3_IBD/11_IBD_final_model.Rmd`.
  **Not a pure rename:** also `n_denominators` -> `n_features` and every
  `$denominators` result access -> `$features`.

**Cautions:** spans two external repos with their own git state (separate
commits/PRs). `summer-students-26` is a teaching repo — confirm no student
cohort is mid-run. Re-verify the 7 analysis notebooks after the denominator
edits; they produce published results.

## 8. S4 slots — DROP `validation_metrics` ONLY

Remove the slot from `setClass TransferablePanelResult`
(`R/panel_class.R:298`) and its `@slot` doc (`:286`); drop the
`validation_metrics = list()` arguments at `R/np_fit.R:286` and
`R/threshold_calibration.R:269`. Regenerate the Rd.

**Consequence:** changes the S4 class definition, so previously saved
`TransferablePanelResult` RDS objects may fail validity checks on load. Note in
NEWS alongside the version bump.

**Doc correction — all three claims were wrong:**
- `training_signature` **is** populated with real content
  (`optimize_panel.R:493`, `optimize_panel_transferable.R:360` — n, p,
  class_balance, pool sizes, cohort labels/counts); it is simply never read by
  the package. Kept.
- `partition_info` is **not** always `list()` —
  `threshold_calibration.R:270` fills it from `panel@control$partition_info`,
  which `optimize_panel_transferable.R:351` sets. Only `np_fit.R:287`
  hardcodes empty. Kept.
- Only `validation_metrics` was genuinely always empty.

**Follow-up:** `summer-students-26/.../3_IBD/5_IBD_optimize_transferable.Rmd:215`
has a commented-out `# val_metrics <- result@validation_metrics` — delete that
line during the decision-7 migration. Lines 250-251 reading `partition_info`
can be uncommented instead.

## 9. `cv_control` / `require_sign_consistency` — DROP `cv_control`, KEEP THE OTHER

After decision 1 removes the CPOP sites, `cv_control` survives only at
`R/transferable_features.R:36` (`@param`), `:58` (signature) and `:92` (the
`modifyList` into the `cv.glmnet` arg list). Remove all three and inline the
fixed arg list. Regenerate `man/select_transferable_features.Rd`. Zero
downstream callers; gives up the `cv.glmnet` tuning escape hatch.

`require_sign_consistency` is **kept**. **Doc correction:** "provably unused" is
wrong — live callers at `5_fetalgrowth/1.3_FGR_features.Rmd:115` and
`summer-students-26/.../3_IBD/2_IBD_transferable_features.Rmd:238`, with
student-facing prose at `:263`, plus tests at `test-feature-selection.R:204`,
`:213`. The redundancy claim *is* correct (`sign_agreement >= 0` is always
TRUE, so `FALSE` == `sign_consistency_threshold = 0`) — document that at
`R/transferable_features.R:27-29` and `:224`.

## 10. Internal man pages — SWITCH ALL INTERNALS TO `@noRd`

Replace `@keywords internal` with `@noRd` in the 141 internal-only roxygen
blocks (103 dot-prefixed), run `devtools::document()`, and `git rm` the
orphaned `man/*.Rd`. The package uses `@noRd` in zero files today, so do it
wholesale in one commit. Shrinks `man/` from 220 pages / 1.2 MB to ~79 pages.

**Caveat accepted:** `@noRd` skips Rd generation, so `R CMD check` no longer
validates `@param`/`@usage` against internal signatures — drift becomes silent.
This lands right after a refactor introducing ~8 new shared helpers
(`.build_roc`, `.confusion_counts`, `.make_fitness_scaffold`, `.run_nsga`,
`.resolve_panel_solution`, `.validate_selection_threshold`,
`.min_features_required`, `.combine_partitions`), so run
`devtools::document()` + `R CMD check` **before** the switch to catch existing
mismatches while they are still checkable.

Roxygen prose stays in source either way.

## 11. `dev/test_parallel_timing.R` — RETIRE

`git rm` the script (198 lines, tracked) and soften CLAUDE.md's "Parallel
Evaluation" note to a qualitative warning: `parallel = TRUE` is usually slower
because of rmoo's per-candidate overhead, and is only worth considering when
per-candidate work is large (CV fitness + large populations). Drop the "~1.5x",
"up to 40x" and "popSize >= 128" specifics and the pointer to the script.

**Rationale:** those numbers were measured against a fitness path this work
replaces (the `.make_fitness_scaffold()` consolidation, decision 5's
unconditional cache, decision 4's nfolds change), so they are no longer
known-good. Retiring removes the stale figures and the maintenance burden
together.

**Doc correction:** `dev/` is tracked, not untracked — both `dev/` files are on
`worktree-codebase-cleanup` and editable from the PR.

---

## Cross-cutting notes

**Decisions interlock — execution order matters:**

1. Decision 1 (retire CPOP) must precede decision 9, which removes what is left
   of `cv_control` afterwards.
2. Decision 11 (retire the timing script) removes one of decision 7's call
   sites — do it first so the migration list is 7, not 8.
3. Decision 10 (`@noRd`) should run **last** among code changes, and only after
   a clean `document()` + `R CMD check`, since it removes that safety net.
4. Decision 8's slot removal pairs with the decision-7 notebook edits (the
   commented-out `validation_metrics` line).
5. Within decision 7, the **shim deletion must come after** the external-repo
   migrations. Those migrations can start at any time — both new names
   (`select_de_features`, `select_batch_associated_features`) already exist, so
   the migrated notebooks work against current code.

**Two PRs, not one.** The cleanup branch (`worktree-codebase-cleanup`, 7
commits) lands first, on its own verification. All 11 decisions execute as a
**follow-up PR** branched off `main` after that merge. Reasons: the FAIL 0 /
PASS 1142 run certifies the branch *as it stands*; the cleanup claims to be
behaviour-preserving while decisions 3/4/5/7/8/9 are deliberate breaks;
decision 10 needs the cleanup PR's clean `R CMD check` as its precondition; and
decision 1 is blocked on an external analysis-repo migration that would
otherwise hold up the whole cleanup.

**Two version bumps, then.** One in the cleanup PR (its 19 export removals
already break API), one in the follow-up with a NEWS entry covering decisions
3 (`scoring_fn`), 5 (both cache knobs), 7 (both shims), 8 (S4 slot) and
9 (`cv_control`). The S4 change is the one to call out — it can invalidate
saved `TransferablePanelResult` RDS objects.

**Three repos are touched:** `biomarkerPanels` (this one),
`biomarkerPanels-analysis` (decisions 1, 7), `summer-students-26` (decisions 7,
8). The latter two need their own commits and pushes.

**Teaching repo is clear.** Confirmed 2026-08-24: no student cohort is
currently running against `summer-students-26`, so the decision-7 edits there
are safe to make.

**CLAUDE.md is untracked** and accumulates edits from decisions 6 and 11, on top
of the ones already listed in the original doc. It cannot travel in the PR —
edit it in the working tree.
