#' Fit a CPOP Panel as a BiomarkerPanelResult
#'
#' Implements CPOP step 3: refits ridge logistic regression on the selected
#' pair features for each cohort, then averages the per-cohort coefficient
#' vectors into a single linear model. At predict time this is identical to
#' averaging the per-cohort log-odds, which is how the original CPOP combines
#' its two cohort-specific models.
#'
#' The returned [BiomarkerPanelResult-class] wraps the averaged coefficients
#' inside a `cv.glmnet` carrier so that [evaluate_panel()],
#' [evaluate_panel_by_cohort()], and threshold tooling all compose with CPOP
#' panels without special-casing.
#'
#' @param cpop_result Output of [select_cpop_features()].
#' @param x_list Training cohort matrices (same shape used for selection).
#' @param y_list Aligned binary response factors (`c("No", "Yes")`).
#' @param lambda_choice `"lambda.min"` (default) or `"lambda.1se"`.
#' @param intercept Logical; passed to [glmnet::cv.glmnet()]. Default `FALSE`
#'   matches CPOP's standard configuration where pair features are already
#'   sample-centred via the ratio.
#' @param assay For `SummarizedExperiment` inputs, the assay name or index.
#' @param cv_control Optional named list forwarded to [glmnet::cv.glmnet()].
#' @return A [BiomarkerPanelResult-class] with `feature_transform =
#'   "pairwise_ratios"` so that [evaluate_panel()] applies the same ratio
#'   construction to validation data.
#' @seealso [select_cpop_features()], [evaluate_panel()]
#' @export
fit_cpop_panel <- function(cpop_result, x_list, y_list,
                           lambda_choice = c("lambda.min", "lambda.1se"),
                           intercept = FALSE,
                           assay = NULL,
                           cv_control = list()) {
  lambda_choice <- match.arg(lambda_choice)
  if (!is.list(cpop_result) ||
      is.null(cpop_result$features) ||
      is.null(cpop_result$base_features)) {
    stop("`cpop_result` must be the output of select_cpop_features().",
         call. = FALSE)
  }
  selected_pairs <- cpop_result$features
  base_features <- cpop_result$base_features
  if (!length(selected_pairs)) {
    stop("`cpop_result$features` is empty; nothing to fit.", call. = FALSE)
  }
  if (length(selected_pairs) < 2L) {
    stop("CPOP step 3 requires at least two selected pair features. ",
         "Re-run select_cpop_features() with relaxed parameters ",
         "(e.g. larger n_features, lower sign_consistency_threshold).",
         call. = FALSE)
  }

  prepared <- .cpop_prepare_inputs(x_list, y_list, assay = assay)
  pair_matrices <- prepared$pair_matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names

  missing_pairs <- setdiff(selected_pairs, colnames(pair_matrices[[1]]))
  if (length(missing_pairs)) {
    stop("Selected CPOP pair features missing from training data: ",
         paste(utils::head(missing_pairs, 5L), collapse = ", "),
         if (length(missing_pairs) > 5L) ", ..." else "",
         call. = FALSE)
  }

  ridge_models <- vector("list", length(pair_matrices))
  for (k in seq_along(pair_matrices)) {
    cv_args <- utils::modifyList(list(
      x = pair_matrices[[k]][, selected_pairs, drop = FALSE],
      y = as.integer(responses[[k]]) - 1L,
      alpha = 0,
      family = "binomial",
      intercept = intercept
    ), cv_control)
    ridge_models[[k]] <- do.call(glmnet::cv.glmnet, cv_args)
  }

  coef_matrix <- vapply(ridge_models, function(m) {
    cm <- as.matrix(stats::coef(m, s = lambda_choice))
    cm[-1L, 1L][selected_pairs]
  }, numeric(length(selected_pairs)))
  if (is.null(dim(coef_matrix))) {
    coef_matrix <- matrix(coef_matrix, nrow = length(selected_pairs))
  }
  rownames(coef_matrix) <- selected_pairs
  colnames(coef_matrix) <- cohort_names

  intercepts <- vapply(ridge_models, function(m) {
    as.matrix(stats::coef(m, s = lambda_choice))[1L, 1L]
  }, numeric(1))

  avg_beta <- rowMeans(coef_matrix)
  avg_intercept <- mean(intercepts)

  shell <- .cpop_build_shell_model(pair_matrices, responses, selected_pairs,
                                   avg_beta, avg_intercept, intercept)

  pooled_n <- sum(vapply(responses, length, integer(1)))
  pooled_y_chr <- unlist(lapply(responses, as.character))

  new(
    "BiomarkerPanelResult",
    base_features = base_features,
    features = selected_pairs,
    metrics = numeric(0),
    objectives = data.frame(
      solution_id = 1L,
      features = I(list(selected_pairs)),
      stringsAsFactors = FALSE
    ),
    control = list(
      feature_transform = "pairwise_ratios",
      regularized = TRUE,
      regularized_alpha = 0,
      positive_class = "Yes",
      cpop = list(
        per_cohort_coefficients = coef_matrix,
        per_cohort_intercepts = stats::setNames(intercepts, cohort_names),
        averaged_coefficients = avg_beta,
        averaged_intercept = avg_intercept,
        lambda_choice = lambda_choice,
        intercept = intercept
      )
    ),
    training_data = list(
      n = pooled_n,
      p = length(base_features),
      class_balance = table(factor(pooled_y_chr, levels = c("No", "Yes"))),
      feature_pool_size = length(selected_pairs),
      num_cohorts = length(pair_matrices)
    ),
    model = shell
  )
}

#' Build a cv.glmnet shell carrying averaged CPOP coefficients
#'
#' Fits a placeholder `cv.glmnet` on pooled cohort data, then overwrites every
#' column of `glmnet.fit$beta` and every entry of `glmnet.fit$a0` with the
#' averaged CPOP coefficients. With identical coefficients at every lambda,
#' downstream `predict(model, newx, s = "lambda.min")` calls used by
#' [evaluate_panel()] produce CPOP's averaged log-odds regardless of which
#' lambda index is selected.
#'
#' @keywords internal
.cpop_build_shell_model <- function(pair_matrices, responses, selected_pairs,
                                    avg_beta, avg_intercept, intercept) {
  pooled_x <- do.call(rbind, lapply(pair_matrices, function(m) {
    m[, selected_pairs, drop = FALSE]
  }))
  pooled_y <- as.integer(unlist(lapply(responses, function(r) {
    as.integer(r) - 1L
  })))

  shell <- glmnet::cv.glmnet(
    x = pooled_x,
    y = pooled_y,
    alpha = 0,
    family = "binomial",
    intercept = intercept
  )

  n_lambda <- ncol(shell$glmnet.fit$beta)
  new_beta <- matrix(avg_beta, nrow = length(selected_pairs), ncol = n_lambda,
                     dimnames = list(selected_pairs, NULL))
  shell$glmnet.fit$beta <- methods::as(new_beta, "CsparseMatrix")
  shell$glmnet.fit$a0 <- stats::setNames(
    rep(avg_intercept, n_lambda),
    paste0("s", seq_len(n_lambda) - 1L)
  )
  shell$biomarkerPanels_meta <- list(
    feature_names = selected_pairs,
    cohort_info = NULL,
    alpha = 0
  )
  shell
}
