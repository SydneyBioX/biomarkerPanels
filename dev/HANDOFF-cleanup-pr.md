# Handoff: land the codebase-cleanup PR, then execute the 11 decisions

Written 2026-08-24 for a fresh Claude Code session. Working directory:
`/dskh/nobackup/harryr/biomarkerPanels`.

Read this file, then `dev/cleanup-decisions-2026-08-RESOLVED.md` (the execution
spec for phase 2). The original survey doc is
`dev/cleanup-decisions-2026-08.md` — treat it as **background only**; several of
its factual claims were checked and found wrong, and RESOLVED.md records the
corrections.

---

## Where things stand

**Already done this session — do not redo:**

- `worktree-harden-spec-at-sens` was merged into `main` as `88cff64` (`--no-ff`).
  It hardens `metric_specificity_at_sensitivity()` so degenerate inputs return
  `0` instead of erroring or returning `NA`, matching its sibling
  `metric_sensitivity_at_specificity()`.
- Post-merge tests on `main` pass: `test-metrics-base.R` 79/79,
  `test-metrics-cohort.R` 18/18, `test-metrics-objectives.R` 35/35.
- The full suite was run on `worktree-codebase-cleanup` at its tip (`884c259`):
  **FAIL 0 | WARN 54 | SKIP 2 | PASS 1142**. The 54 warnings are benign
  small-fixture noise (`glmnet` "fewer than 8 observations", `nproc` empirical
  threshold fallback), not regressions.
- A trial merge of the cleanup branch into the new `main` is **conflict-free**
  (`git merge-tree`), and the merged `metric_specificity_at_sensitivity()` was
  inspected: it correctly combines the hardening with the branch's `.build_roc()`
  extraction.

**Current git state:**

- `main` is at `88cff64`, **2 commits ahead of `origin/main`**, working tree clean.
- Two locked worktrees, both clean, both branched off `22a975e`:
  - `.claude/worktrees/harden-spec-at-sens` — now merged, can be removed when convenient.
  - `.claude/worktrees/codebase-cleanup` — branch `worktree-codebase-cleanup`,
    7 commits, 103 files, +1759/−4559. **This is the phase-1 PR.**
- Also present and unrelated: local branches `codebase-audit-fixes`, `pr-5`.

**Two constraints that apply throughout:**

- **Never `git push` from this session.** The user's SSH key is passphrase-
  protected and ssh-agent is not unlocked here — push hangs. Ask the user to
  push from their own terminal.
- **`gh` is not installed.** PRs get opened in the browser by the user. Remote
  is `git@github.com:SydneyBioX/biomarkerPanels.git`.
- `dev/` is **gitignored** (`.gitignore:22`); the tracked `dev/` files were
  force-added. This handoff and RESOLVED.md are therefore invisible to git
  until `git add -f`'d.

---

## Phase 1 — land the cleanup PR

Do these in order, in the worktree (the branch is checked out there, so the
merge must happen inside it):

```bash
cd /dskh/nobackup/harryr/biomarkerPanels/.claude/worktrees/codebase-cleanup
git merge main          # verified conflict-free
```

1. **Merge `main` into the branch.** This picks up the harden commit and makes
   the PR diff read against current `main`.
2. **Re-run the metric tests on the merged tree.** This is the one gap in
   existing verification: the merged `tests/testthat/test-metrics-base.R`
   (the branch's −72-line rewrite combined with the harden commit's edits) has
   never been executed by anyone.
   ```bash
   R -q -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-metrics-base.R')"
   ```
   Also worth re-running `test-metrics-cohort.R` and `test-metrics-objectives.R`.
3. **`R CMD check --no-tests`** on the merged tree. Nobody has run check
   end-to-end on this branch, and it rewrites `NAMESPACE`, the `Collate:` field
   and ~50 `man/` files. Its tip commit (`884c259`) exists solely to fix a check
   complaint, so check is the gate this work was built against. Use
   `devtools::check(args = "--no-tests")` — the suite runs separately.
4. **Version bump commit** in `DESCRIPTION`. The branch removes 19 exports,
   which is breaking for two live consumers, and currently carries no bump.
5. **CLAUDE.md edits.** The file is **untracked**, so it cannot travel in the
   PR — edit it in the working tree at
   `/dskh/nobackup/harryr/biomarkerPanels/CLAUDE.md`. Three sources of edits:
   - From the original doc: remove `select_panel_top_sensitivity()` and
     `select_panel_by_pathway()` from the "Pareto Solution Selection" bullet
     (keep `select_panel_inclusion_frequency()` — it has 2 live call sites);
     note in the Metric Registry section that `metric_*` functions are now
     internal and users reach them via registry names; update the Architecture
     section for the file moves (`de_features_helpers.R` gone,
     `stable_features.R` / `pareto_feature_stability.R` renames,
     `panel_selection.R` / `glm_fitting.R` split out of `fitness_cache_utils.R`,
     `utils.R` gone).
   - From decision 6: rewrite "Rule-Out Diagnostics" and remove
     `define_ruleout_objectives()` from "Main Entry Points" — **that function
     never existed**. See RESOLVED.md §6 for what to say instead.
   - From decision 11: soften the "Parallel Evaluation" note to a qualitative
     warning, dropping the "~1.5x", "up to 40x" and "popSize >= 128" figures and
     the pointer to `dev/test_parallel_timing.R`. (Strictly this belongs with
     phase 2, which deletes that script — but it is a CLAUDE.md edit, so doing
     it here keeps all CLAUDE.md work in one pass. Either is fine; just do not
     do it twice.)
6. **Hand off for push.** Ask the user to push `main` and
   `worktree-codebase-cleanup` from their own terminal, then open the PR in the
   browser.

**PR body should cover** the 6 commit themes — dead-export deletion, un-exporting
15 registry-only metrics, feature-selection consolidation, metrics/evaluation
dedup, the shared fitness scaffold + NSGA launch path, and the file
reorganisation — plus a pointer to `dev/cleanup-decisions-2026-08.md`, which
*is* tracked on the branch. Mention the FAIL 0 / PASS 1142 run.

`/code-review ultra <PR#>` is available once the PR exists, but it is
user-triggered and billed — do not attempt to launch it.

---

## Phase 2 — the 11 decisions, as a follow-up PR

**Do not start this until the cleanup PR has merged.** The reasoning is in
RESOLVED.md's cross-cutting notes; the short version is that the cleanup's green
run certifies the branch as it stands, the cleanup claims to be
behaviour-preserving while six of the decisions are deliberate breaks, and
decision 10 needs the cleanup PR's clean `R CMD check` as its precondition.

Branch off `main` after the merge. `git add -f
dev/cleanup-decisions-2026-08-RESOLVED.md` onto that branch so the spec travels
with the work.

**RESOLVED.md is the full spec** — file paths and line numbers for every edit.
Summary of the choices:

| # | Item | Decision |
|---|---|---|
| 1 | CPOP subsystem | Retire entirely |
| 2 | `loco` / `within_cohort_rotating` | Keep both, add a `loco` test |
| 3 | `optimize_panel(scoring_fn=)` | Remove outright |
| 4 | Inner-CV nfolds | Unify to the 5-cap everywhere |
| 5 | `cache_fitness` / `cache_max_entries` | Remove both |
| 6 | `define_ruleout_objectives()` | Leave removed, fix docs |
| 7 | Both deprecated shims | Migrate 19 sites, drop both |
| 8 | S4 slots | Drop `validation_metrics` only |
| 9 | `cv_control` / `require_sign_consistency` | Drop `cv_control`, keep the other |
| 10 | 141 internal man pages | Switch all to `@noRd` |
| 11 | `dev/test_parallel_timing.R` | Retire |

**Execution order is not free — the decisions interlock:**

1. Decision 1 (retire CPOP) **before** decision 9 — CPOP removal deletes 5 of
   `cv_control`'s 7 sites, leaving only the three in `transferable_features.R`.
2. Decision 11 (retire the timing script) **before** decision 7 — it deletes one
   of the `get_top_de_features` call sites, making the migration list 19, not 20.
3. Decision 7's external-repo migrations can happen any time (both new names
   already exist), but the **shim deletion must come after** them.
4. Decision 10 (`@noRd`) **last** among code changes, and only after a clean
   `document()` + `R CMD check` — it removes exactly that safety net, so any
   pre-existing `@param`/`@usage` mismatches must be caught while they are still
   checkable.
5. Decision 8's slot removal pairs with a decision-7 notebook edit (a
   commented-out `validation_metrics` line in the teaching repo).

**Blockers and cautions:**

- **Decision 1 is blocked on an external migration.**
  `biomarkerPanels-analysis/5_fetalgrowth/2_FGR_ratios.Rmd:85` calls
  `select_cpop_features()`. Migrate it to the standalone
  `/dskh/nobackup/harryr/CPOP` repo, or pin it to an older package version,
  before the retirement lands.
- **Teaching repo is clear.** Confirmed 2026-08-24 with the user: no student
  cohort is currently running against `summer-students-26`, so decision 7's
  edits there are safe.
- **Decision 4 changes results.** `fit_panel()` will pick a different
  `lambda.min`, so final-model numbers and `panel_metrics()` shift. NSGA output
  is unchanged. Affected analysis results need re-running.
- **Decision 8 can invalidate saved RDS.** Removing the slot changes the
  `TransferablePanelResult` class definition. Call it out in NEWS.
- **Three repos are touched**, each with its own git state and its own push:
  `biomarkerPanels`, `biomarkerPanels-analysis` (decisions 1, 7),
  `summer-students-26` (decisions 7, 8).
- **Second version bump** with a NEWS entry covering decisions 3, 5, 7, 8, 9.

---

## Build and test notes

```bash
# Fast, preferred during development
R -e "devtools::load_all(); testthat::test_file('tests/testthat/test-metric-functions.R')"

# Full suite — slow, optimization tests run real NSGA loops
R -e "devtools::test()"
```

`devtools::load_all()` is sufficient for development. If the user needs the
shared library refreshed at `/enna/nobackup/biostat/Rpackages/v4`, note that an
rsession (PID 2826545 as of this session) holds `biomarkerPanels.so` — kill
blocking PIDs *before* reinstalling, never use `--no-lock` while the `.so` is
held (it wipes the install dir partway and leaves the package unloadable). The
user chose to **skip** the reinstall for the harden fix, so the installed
package is currently behind `main`.

---

## STATUS UPDATE 2026-08-24 (phase 1 executed)

Phase 1 is **complete up to the push**. On `worktree-codebase-cleanup`:

- `677f0dd` merge of `main` (conflict-free, as predicted)
- `bee1e51` fixes two roxygen Rd-link escapes found by `R CMD check`
  (`R/fitness_eval.R`, `R/utils_validation.R`)
- `b3e0da4` version bump `0.0.0.9000` -> `0.1.0` + NEWS.md entry

Verification on the merged tree: test-metrics-base 72/72, -cohort 13/13,
-objectives 34/34 (FAIL 0 everywhere); `devtools::check(--no-tests)` at
`b3e0da4`: **0 errors**, remaining warning/note are environmental only
(missing qpdf binary; timestamp verification). CLAUDE.md in the main checkout
was updated (all step-5 edits, plus removal of the deleted
`metric_conditional_score_shift` / `metric_cohort_transferability.R`
references).

Still to do: user pushes `worktree-codebase-cleanup` (`main` is already
pushed) and opens the PR in the browser. Phase 2 remains gated on that PR
merging.

---

## STATUS UPDATE 2026-08-25 (phase 2 executed)

All 11 decisions are done on branch `phase2-decisions` (14 commits, tip = the
0.2.0 bump). Verification at tip: full suite FAIL 0 | WARN 59 | SKIP 2 |
PASS 1132; `R CMD check --no-tests` 0 errors / 1 warning (qpdf, environmental)
/ 0 notes. man/ 212 -> 74 pages after the @noRd sweep. External repos:
`summer-students-26/Summer-Biomarker-2026` committed on `taine` (734237f);
`biomarkerPanels-analysis` migrations left uncommitted in the working tree
(user WIP was staged there) with recovery patch
`dev-shim-migration-2026-08-24.patch`. Shared NFS lib installed at 0.2.0.
Still to do: user pushes `phase2-decisions` + teaching-repo `taine`, opens the
PR. Optional follow-up: ~18 untagged file-overview man topics are arguably
internal (future @noRd pass).
