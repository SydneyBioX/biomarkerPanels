#' Fast GLM Design Construction and Fitting
#'
#' Internal helpers for building cohort-adjusted design matrices and fitting
#' binomial GLMs via [stats::glm.fit()] on the optimization hot path. Cohort
#' membership enters as dummy columns; `predict_cohort = "reference"` scores
#' new samples at the reference cohort (all dummies zero) while `"observed"`
#' uses each sample's own cohort.
#'
#' @name glm_fitting
#' @keywords internal
NULL

.cohort_dummy_matrix <- function(cohort, levels) {
  n <- length(cohort)
  if (length(levels) <= 1L) {
    return(NULL)
  }
  cohort_factor <- factor(cohort, levels = levels)
  dummies <- matrix(0, nrow = n, ncol = length(levels) - 1L)
  colnames(dummies) <- paste0(".cohort_", seq_len(ncol(dummies)))
  for (j in seq_along(levels[-1L])) {
    dummies[, j] <- ifelse(
      is.na(cohort_factor),
      NA_real_,
      as.numeric(cohort_factor == levels[-1L][[j]])
    )
  }
  dummies
}

.prepare_glm_design_terms <- function(cohort_train = NULL, cohort_new = NULL,
                                      predict_cohort = c("observed", "reference")) {
  predict_cohort <- match.arg(predict_cohort)
  if (is.null(cohort_train) || length(unique(cohort_train)) <= 1L) {
    return(list(train_extra = NULL, new_extra = NULL))
  }

  cohort_train <- factor(cohort_train)
  levels_train <- levels(cohort_train)
  train_extra <- .cohort_dummy_matrix(cohort_train, levels_train)

  if (is.null(cohort_new)) {
    cohort_new <- cohort_train
  }

  if (predict_cohort == "reference") {
    new_extra <- matrix(0, nrow = length(cohort_new),
                        ncol = length(levels_train) - 1L)
    colnames(new_extra) <- colnames(train_extra)
  } else {
    new_extra <- .cohort_dummy_matrix(cohort_new, levels_train)
  }

  list(
    train_extra = train_extra,
    new_extra = new_extra,
    levels = levels_train,
    predict_cohort = predict_cohort
  )
}

.prepare_cv_glm_design_terms <- function(cohort, fold_ids) {
  if (is.null(cohort)) {
    return(NULL)
  }
  unique_folds <- sort(unique(fold_ids))
  terms <- lapply(unique_folds, function(fold) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    .prepare_glm_design_terms(
      cohort_train = cohort[train_idx],
      cohort_new = cohort[test_idx],
      predict_cohort = "observed"
    )
  })
  names(terms) <- as.character(unique_folds)
  terms
}

.build_glm_design_matrix <- function(x, extra = NULL) {
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  design <- cbind("(Intercept)" = 1, x_mat)
  if (!is.null(extra) && ncol(extra) > 0L) {
    design <- cbind(design, extra)
  }
  design
}

.fit_predict_binomial_glm <- function(x_train, truth, x_new = NULL,
                                      cohort_train = NULL, cohort_new = NULL,
                                      predict_cohort = c("observed", "reference"),
                                      design_terms = NULL) {
  predict_cohort <- match.arg(predict_cohort)
  if (is.null(x_train) || ncol(x_train) == 0L) {
    stop("No features selected for GLM scoring.", call. = FALSE)
  }
  if (length(unique(truth)) < 2L) {
    stop("Response contains only one class. Cannot fit classification model.",
         call. = FALSE)
  }
  if (is.null(x_new)) {
    x_new <- x_train
  }
  if (is.null(design_terms)) {
    design_terms <- .prepare_glm_design_terms(
      cohort_train = cohort_train,
      cohort_new = if (is.null(cohort_new)) cohort_train else cohort_new,
      predict_cohort = predict_cohort
    )
  }

  x_design <- .build_glm_design_matrix(x_train, design_terms$train_extra)
  new_design <- .build_glm_design_matrix(x_new, design_terms$new_extra)
  y_vec <- as.integer(truth) - 1L

  fit <- stats::glm.fit(
    x = x_design,
    y = y_vec,
    family = stats::binomial()
  )
  coefs <- fit$coefficients
  # Rank-deficient designs give NA coefficients; zeroing them reproduces
  # stats::predict.glm's treatment of aliased terms.
  coefs[is.na(coefs)] <- 0
  preds <- stats::binomial()$linkinv(drop(new_design %*% coefs))
  as.numeric(preds)
}
