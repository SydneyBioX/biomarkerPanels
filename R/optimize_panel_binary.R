#' Optimize Biomarker Panels with Binary Genetic Algorithm
#'
#' Alternative to [optimize_panel()] that uses binary encoding for feature
#' selection. Each chromosome is a binary vector where 1 indicates feature
#' inclusion and 0 indicates exclusion. This eliminates the continuous-to-
#' discrete mapping problem of the standard NSGA-II approach.
#'
#' @param x Matrix-like object, `SummarizedExperiment`, or list of matrices /
#'   experiments representing one or more cohorts.
#' @param y Binary response (`factor`, `character`, or `logical`) aligned with
#'   `x`. When `x` is a list, `y` must be a list of the same length.
#' @param loss Primary loss function to optimize (default: "auc").
#' @param max_features Maximum number of biomarkers permitted in a panel.
#' @param min_features Minimum number of biomarkers required. Default 2.
#' @param feature_pool Optional subset of feature identifiers to consider.
#' @param cohort_aggregator Transformation applied to cohort feature matrices.
#'   Defaults to `"pairwise_ratios"`.
#' @param feature_alignment Strategy for aligning features across cohorts.
#' @param ga_control Named list of arguments passed to [GA::ga()]. Defaults
#'   to `list(popSize = 100, maxiter = 100, pmutation = 0.1, pcrossover = 0.8)`.
#' @param assay For `SummarizedExperiment` inputs, assay name or index to use.
#' @param seed Optional integer seed for reproducibility.
#' @param fitness_cv Logical; if `TRUE` (default), use cross-validation when
#'   evaluating candidate solutions.
#' @param fitness_cv_folds Number of cross-validation folds for fitness
#'   evaluation when `fitness_cv = TRUE`. Default is 5.
#' @param regularized_alpha Elastic net mixing parameter. Default is 0.5.
#' @param parallel Whether to run fitness evaluations in parallel. Default FALSE.
#' @return A `BiomarkerPanelResult` with the best solution found.
#' @export
optimize_panel_binary <- function(x, y,
                                   loss = "auc",
                                   max_features = 5L,
                                   min_features = 2L,
                                   feature_pool = NULL,
                                   cohort_aggregator = "pairwise_ratios",
                                   feature_alignment = c("intersection", "majority", "impute_median"),
                                   ga_control = list(),
                                   assay = NULL,
                                   seed = NULL,
                                   fitness_cv = TRUE,
                                   fitness_cv_folds = 5L,
                                   regularized_alpha = 0.5,
                                   parallel = FALSE) {
  feature_alignment <- match.arg(feature_alignment)
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop("The 'GA' package is required. Install it via install.packages('GA').",
         call. = FALSE)
  }

  # Set seed if provided

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --------------------------------------------------------------------------
  # Preprocess inputs using shared helper (same as optimize_panel)
  # --------------------------------------------------------------------------
  inputs <- .prepare_cohort_inputs(
    x, y,
    assay = assay,
    aggregator = cohort_aggregator,
    feature_subset = feature_pool,
    feature_alignment = feature_alignment
  )

  x_combined <- inputs$x
  truth <- inputs$truth
  cohort_vec <- inputs$cohort
  feature_names_pool <- colnames(x_combined)

  n_features <- length(feature_names_pool)
  if (n_features == 0L) {
    stop("No features available after filtering.", call. = FALSE)
  }
  if (n_features < min_features) {
    stop("Number of features (", n_features, ") is less than min_features (",
         min_features, ").", call. = FALSE)
  }

  # Limit max_features to available features
  max_features <- min(max_features, n_features)

  # Create CV folds for fitness evaluation
  if (fitness_cv) {
    fold_ids <- .create_stratified_folds(truth, fitness_cv_folds)
  } else {
    fold_ids <- NULL
  }

  # --------------------------------------------------------------------------
  # Define fitness function for binary GA
  # --------------------------------------------------------------------------
  # GA maximizes fitness, so return the loss metric directly
  loss_entry <- .loss_registry[[loss]]
  if (is.null(loss_entry)) {
    stop("Unknown loss function: ", loss, ". Use loss_registry() to see available options.",
         call. = FALSE)
  }
  loss_fn <- loss_entry$fun
  loss_direction <- loss_entry$direction

  fitness_fn <- function(chromosome) {
    selected_idx <- which(chromosome == 1)
    n_selected <- length(selected_idx)

    # Penalty for invalid panel sizes
    if (n_selected < min_features || n_selected > max_features) {
      return(-1e6)
    }

    selected_features <- feature_names_pool[selected_idx]
    x_selected <- x_combined[, selected_idx, drop = FALSE]

    # Compute scores using regularized regression
    if (fitness_cv && !is.null(fold_ids)) {
      scores <- .compute_cv_scores(
        x_selected = x_selected,
        truth = truth,
        fold_ids = fold_ids,
        cohort = cohort_vec,
        regularized = TRUE,
        alpha = regularized_alpha
      )
    } else {
      scores <- .default_scoring_fn(
        x_selected = x_selected,
        selected_features = selected_features,
        truth = truth,
        cohort = cohort_vec,
        regularized = TRUE,
        alpha = regularized_alpha
      )
    }

    # Compute loss value
    loss_value <- tryCatch({
      loss_fn(
        truth = truth,
        scores = scores,
        selected = selected_features
      )
    }, error = function(e) {
      # Track failures in verbose mode - penalty will guide search away
      # from problematic solutions
      if (getOption("biomarkerPanels.verbose", FALSE)) {
        message("Fitness evaluation failed for candidate: ", conditionMessage(e))
      }
      if (loss_direction == "maximize") -1e6 else 1e6
    })

    # Return fitness (GA maximizes, so negate if loss should be minimized)
    if (loss_direction == "minimize") {
      return(-loss_value)
    } else {
      return(loss_value)
    }
  }

  # --------------------------------------------------------------------------
  # Run binary GA
  # --------------------------------------------------------------------------
  ga_defaults <- list(
    popSize = 100L,
    maxiter = 100L,
    pmutation = 0.1,
    pcrossover = 0.8,
    elitism = base::max(1, round(100 * 0.05)),
    monitor = FALSE
  )
  ga_params <- modifyList(ga_defaults, ga_control)

  ga_result <- GA::ga(
    type = "binary",
    fitness = fitness_fn,
    nBits = n_features,
    popSize = ga_params$popSize,
    maxiter = ga_params$maxiter,
    pmutation = ga_params$pmutation,
    pcrossover = ga_params$pcrossover,
    elitism = ga_params$elitism,
    parallel = parallel,
    monitor = ga_params$monitor,
    seed = seed
  )

  # Extract best solution
  best_chromosome <- as.integer(ga_result@solution[1, ])
  selected_idx <- which(best_chromosome == 1)
  selected_features <- feature_names_pool[selected_idx]
  x_selected <- x_combined[, selected_idx, drop = FALSE]

  # --------------------------------------------------------------------------
  # Fit final model and compute metrics
  # --------------------------------------------------------------------------
  final_model <- .fit_final_model_regularized(
    x_selected, truth, cohort_vec, alpha = regularized_alpha
  )

  # Compute final scores for metrics
  if (fitness_cv && !is.null(fold_ids)) {
    final_scores <- .compute_cv_scores(
      x_selected = x_selected,
      truth = truth,
      fold_ids = fold_ids,
      cohort = cohort_vec,
      regularized = TRUE,
      alpha = regularized_alpha
    )
  } else {
    final_scores <- .default_scoring_fn(
      x_selected = x_selected,
      selected_features = selected_features,
      truth = truth,
      cohort = cohort_vec,
      regularized = TRUE,
      alpha = regularized_alpha
    )
  }

  # Compute all standard metrics
  metrics <- list(
    num_features = length(selected_features),
    auc = loss_auc(truth, final_scores, selected_features),
    sensitivity = loss_sensitivity(truth, final_scores, selected_features),
    specificity = loss_specificity(truth, final_scores, selected_features),
    best_fitness = ga_result@fitnessValue
  )

  # Map aggregated features (e.g., "A--B") back to original feature names
  original_features <- unique(unlist(strsplit(selected_features, "--", fixed = TRUE)))

  # --------------------------------------------------------------------------
  # Build result object
  # --------------------------------------------------------------------------
  objectives_df <- data.frame(
    solution_id = 1L,
    auc = metrics$auc,
    sensitivity = metrics$sensitivity,
    specificity = metrics$specificity,
    num_features = metrics$num_features,
    ga_fitness = metrics$best_fitness,
    stringsAsFactors = FALSE
  )

  control <- list(
    method = "binary_ga",
    loss = loss,
    max_features = max_features,
    min_features = min_features,
    cohort_aggregator = cohort_aggregator,
    feature_alignment = feature_alignment,
    fitness_cv = fitness_cv,
    fitness_cv_folds = fitness_cv_folds,
    regularized_alpha = regularized_alpha,
    ga_control = ga_params,
    n_generations = ga_result@iter,
    seed = seed,
    # Store aggregated feature names for reference
    aggregated_features = selected_features
  )

  # Convert metrics list to named numeric vector
  metrics_vec <- c(
    num_features = as.numeric(metrics$num_features),
    auc = as.numeric(metrics$auc),
    sensitivity = as.numeric(metrics$sensitivity),
    specificity = as.numeric(metrics$specificity),
    best_fitness = as.numeric(metrics$best_fitness)
  )

  methods::new(
    "BiomarkerPanelResult",
    features = original_features,
    metrics = metrics_vec,
    objectives = objectives_df,
    control = control,
    training_data = list(),
    model = final_model
  )
}
