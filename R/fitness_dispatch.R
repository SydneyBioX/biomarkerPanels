#' Fitness Strategy Dispatch
#'
#' Maps [optimize_panel()]'s `fitness_mode` onto one of the four fitness
#' factories (`R/fitness_cv.R`, `R/fitness_validation.R`,
#' `R/fitness_rotating_validation.R`, `R/fitness_loco.R`), handing each the
#' data shape it expects. Kept out of `optimize_panel.R` so the pipeline reads
#' as one pass.
#'
#' @name fitness_dispatch
NULL

#' Dispatch to the Fitness Factory for a Given fitness_mode
#'
#' Picks among the four fitness strategies and hands each the data shape it
#' needs: the pooled train + validation matrix for `"cv"` / `"in_sample"` /
#' `"loco"` / `"within_cohort_rotating"`, and the two partitions separately for
#' `"within_cohort_val"`.
#'
#' @param fitness_mode Resolved mode string.
#' @param x_pool,truth,cohort Pooled train + validation data.
#' @param train_inputs,val_inputs Prepared partitions (`val_inputs` may be
#'   `NULL`).
#' @param train_ratio,val_ratio Partition shares, used to derive the rotating
#'   split ratio.
#' @param seed Optional seed forwarded to the rotating split generator.
#' @param ... Remaining arguments shared by every factory.
#' @return List with `wrapper` and `evaluate`, as returned by the factories.
#' @noRd
.make_optimize_panel_fitness <- function(fitness_mode,
                                         x_pool, truth, cohort,
                                         train_inputs, val_inputs,
                                         feature_pool, max_features,
                                         objectives, constraints,
                                         regularized, regularized_alpha,
                                         feature_transform,
                                         min_features_required,
                                         selection_threshold,
                                         fitness_cv_folds,
                                         n_val_splits,
                                         train_ratio, val_ratio,
                                         seed) {
  common <- list(
    feature_pool = feature_pool,
    max_features = max_features,
    objectives = objectives,
    constraints = constraints,
    regularized = regularized,
    alpha = regularized_alpha,
    feature_transform = feature_transform,
    min_features_required = min_features_required,
    selection_threshold = selection_threshold
  )

  switch(
    fitness_mode,
    cv = ,
    in_sample = do.call(.make_cv_fitness, c(
      list(x = x_pool, y = truth, cohort = cohort,
           fitness_cv = identical(fitness_mode, "cv"),
           fitness_cv_folds = fitness_cv_folds),
      common
    )),
    within_cohort_val = do.call(.make_validation_fitness, c(
      list(train_x = train_inputs$x,
           train_y = train_inputs$truth,
           train_cohort = train_inputs$cohort,
           val_x = val_inputs$x,
           val_y = val_inputs$truth,
           val_cohort = val_inputs$cohort),
      common
    )),
    within_cohort_rotating = {
      # Feature selection already happened on the full input, so pooling
      # train + val here does not leak anything the pool did not already see.
      splits <- .generate_stratified_splits(
        y = truth,
        cohort = cohort,
        n_splits = n_val_splits,
        val_ratio = val_ratio / (train_ratio + val_ratio),
        seed = seed
      )
      do.call(.make_rotating_validation_fitness, c(
        list(pool_x = x_pool, pool_y = truth, pool_cohort = cohort,
             splits = splits),
        common
      ))
    },
    loco = {
      n_cohorts <- nlevels(droplevels(cohort))
      if (n_cohorts < 2L) {
        stop(
          "fitness_mode = 'loco' requires >=2 cohorts but only ", n_cohorts,
          " was available after partitioning.",
          call. = FALSE
        )
      }
      do.call(.make_loco_fitness, c(
        list(x = x_pool, y = truth, cohort = cohort),
        common
      ))
    }
  )
}
