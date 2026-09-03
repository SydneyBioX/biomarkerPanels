#' Annotation-free co-expression diversity for Pareto solutions
#'
#' Curated pathway-diversity scores (e.g. counting MSigDB Hallmark gene sets a
#' panel touches) collapse to `NA` whenever the feature identifiers are not gene
#' symbols -- a routine situation for proteomics panels (SomaScan/aptamer target
#' names, UniProt accessions). [compute_diversity()] is an annotation-free
#' replacement: it clusters the training feature pool into co-expression modules
#' and scores each panel by how many distinct modules it spans.

#' Score Pareto solutions by co-expression diversity
#'
#' Appends a co-expression diversity score to a [summarize_solutions()] data
#' frame. The feature pool stored in the `OptimizationResult` is clustered once
#' into `n_modules` co-expression modules (correlation distance `1 - |r|`,
#' hierarchical clustering), and each solution is scored as the number of
#' distinct modules its base features span divided by the panel size. The result
#' lies in `(0, 1]`: values near `1` indicate the panel draws from many
#' independent co-expression programs (spread / non-redundant), while low values
#' indicate the features co-cluster (co-regulated / redundant).
#'
#' Because every pool feature belongs to exactly one module (there is no
#' "unmapped" bucket), the score is never `NA` for panels whose features are in
#' the pool -- unlike annotation-based pathway diversity on non-symbol features.
#'
#' @param solutions A data frame from [summarize_solutions()] (must contain a
#'   `solution_id` column). As a convenience, an `OptimizationResult` may be
#'   passed directly, in which case its summary is computed internally.
#' @param optimization_result The `OptimizationResult` the `solutions` summary
#'   came from. Supplies both the training feature matrix (`@aggregated_x`, whose
#'   columns are the base-feature pool) and the per-solution base features.
#'   Optional when `solutions` is itself an `OptimizationResult`.
#' @param n_modules Number of co-expression modules to cut the pool into
#'   (default `10`). Capped at the number of usable pool features (a warning is
#'   issued if reduced). Note the score can only reach `1.0` when
#'   `n_modules >= panel size`, and degenerates to all-`1.0` once
#'   `n_modules >= pool size`.
#' @param method Clustering method. Currently only `"hierarchical"` is
#'   supported. (See the TODO in Details -- as an open-source project we should
#'   add alternative methods such as k-means.)
#' @param linkage Agglomeration method passed to [stats::hclust()]
#'   (default `"ward.D2"`).
#' @param cor_method Correlation type for the feature distance, one of
#'   `"pearson"` (default) or `"spearman"`.
#' @param column Name of the appended diversity column
#'   (default `"coexpr_diversity"`).
#'
#' @details
#' **Multi-cohort handling.** When the optimization spans multiple cohorts,
#' `@aggregated_x` stacks all cohorts (and, for the transferable pipeline, the
#' train + validation partitions). Correlating that pooled matrix directly would
#' conflate genuine within-cohort co-expression with between-cohort batch
#' mean-shifts (the effect \code{metric_cohort_leakage} is designed to detect),
#' biasing features to look more redundant. To avoid this, each feature is
#' **centered within its cohort** before the correlation matrix is computed, so
#' modules reflect co-expression rather than batch structure. For a single
#' cohort this is a no-op (correlation is invariant to centering).
#'
#' **TODO (open-source extension point).** Only hierarchical clustering is
#' implemented. Other clustering strategies -- k-means, DBSCAN, graph/community
#' detection on the correlation network -- would be natural additions, ideally
#' exposed via a small clustering-method registry mirroring the feature-transform
#' registry in `R/feature_transforms.R` (`new.env()` + `register_*()` + getter,
#' with built-ins registered inline). Until then `method` accepts a single value.
#'
#' @return The `solutions` data frame with one additional numeric column
#'   (named by `column`) giving each solution's co-expression diversity in
#'   `(0, 1]`. Row order and all existing columns are preserved. Solutions whose
#'   features cannot be mapped to the pool receive `NA`.
#'
#' @seealso [optimize_panel()], [summarize_solutions()],
#'   [analyze_feature_stability()]
#'
#' @examples
#' \dontrun{
#' opt <- optimize_panel(
#'   x = list(Cohort1 = x1, Cohort2 = x2),
#'   y = list(Cohort1 = y1, Cohort2 = y2),
#'   objectives = define_objectives(metrics = c("auc", "cohort_auc_gap")),
#'   feature_pool = ruv_features, max_features = 5,
#'   feature_transform = "pairwise_ratios", algorithm = "NSGA-III",
#'   fitness_mode = "within_cohort_val", train_ratio = 0.7, val_ratio = 0.2
#' )
#' solutions <- summarize_solutions(opt)
#' solutions <- compute_diversity(solutions, opt)
#' solutions[order(-solutions$coexpr_diversity), ]
#' }
#'
#' @importFrom stats cor hclust cutree as.dist var
#' @export
compute_diversity <- function(solutions,
                              optimization_result,
                              n_modules = 10L,
                              method = c("hierarchical"),
                              linkage = "ward.D2",
                              cor_method = c("pearson", "spearman"),
                              column = "coexpr_diversity") {
  method <- match.arg(method)
  cor_method <- match.arg(cor_method)

  # Convenience front-door: allow an OptimizationResult as the first argument.
  if (inherits(solutions, "OptimizationResult")) {
    if (missing(optimization_result) || is.null(optimization_result)) {
      optimization_result <- solutions
    }
    solutions <- summarize_solutions(optimization_result)
  }

  if (!is.data.frame(solutions)) {
    stop("`solutions` must be a data frame from summarize_solutions() ",
         "(or an OptimizationResult).", call. = FALSE)
  }
  if (!"solution_id" %in% names(solutions)) {
    stop("`solutions` must contain a `solution_id` column.", call. = FALSE)
  }
  if (missing(optimization_result) ||
        !inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult from ",
         "optimize_panel().", call. = FALSE)
  }
  if (!is.character(column) || length(column) != 1L || !nzchar(column)) {
    stop("`column` must be a single non-empty string.", call. = FALSE)
  }
  if (column %in% names(solutions)) {
    stop("Column '", column, "' already exists in `solutions`; choose a ",
         "different `column` name.", call. = FALSE)
  }

  mat <- optimization_result@aggregated_x
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    warning("`optimization_result` has no training matrix (@aggregated_x); ",
            "returning NA diversity.", call. = FALSE)
    solutions[[column]] <- NA_real_
    return(solutions)
  }

  modules <- .coexpr_modules(
    mat = mat,
    cohort = optimization_result@aggregated_cohort,
    n_modules = n_modules,
    linkage = linkage,
    cor_method = cor_method
  )

  if (is.null(modules)) {
    solutions[[column]] <- NA_real_
    return(solutions)
  }

  # Recover each solution's base features, matched back to the summary by id.
  sol_df <- optimization_result@solutions
  idx <- match(solutions$solution_id, sol_df$solution_id)
  bf_list <- sol_df$base_features[idx]

  solutions[[column]] <- vapply(
    bf_list,
    function(bf) .coexpr_diversity(bf, modules),
    numeric(1)
  )
  solutions
}

#' Cluster a feature pool into co-expression modules
#'
#' Centers each feature within its cohort (to strip batch mean-shifts), drops
#' zero-variance columns, then cuts a `1 - |cor|` hierarchical clustering into
#' modules. Returns a named integer vector (feature -> module id) or `NULL` if
#' fewer than two usable features remain.
#'
#' @noRd
.coexpr_modules <- function(mat, cohort, n_modules, linkage, cor_method) {
  centered <- .center_within_cohort(mat, cohort)

  # Constant columns yield NA correlations and break hclust -- drop them.
  vars <- apply(centered, 2L, stats::var)
  keep <- is.finite(vars) & vars > 0
  if (any(!keep)) {
    warning(sum(!keep), " constant feature(s) dropped before clustering; ",
            "they map to NA modules and are ignored in scoring.", call. = FALSE)
    centered <- centered[, keep, drop = FALSE]
  }

  p <- ncol(centered)
  if (p < 2L) {
    warning("Fewer than two usable features in the pool; cannot build ",
            "co-expression modules. Returning NA diversity.", call. = FALSE)
    return(NULL)
  }

  K <- as.integer(n_modules)
  if (K < 1L) {
    stop("`n_modules` must be a positive integer.", call. = FALSE)
  }
  if (K > p) {
    warning("`n_modules` (", n_modules, ") exceeds the number of usable ",
            "features (", p, "); capping at ", p, ".", call. = FALSE)
    K <- p
  }

  cmat <- stats::cor(centered, method = cor_method)
  hc <- stats::hclust(stats::as.dist(1 - abs(cmat)), method = linkage)
  stats::cutree(hc, k = K)
}

#' Subtract per-cohort feature means
#'
#' @param mat Numeric matrix (samples x features).
#' @param cohort Cohort membership for the rows of `mat`, or `NULL`.
#' @return `mat` with each feature centered within its cohort. A no-op for a
#'   single cohort (correlation is invariant to centering).
#' @noRd
.center_within_cohort <- function(mat, cohort) {
  if (is.null(cohort)) {
    return(mat)
  }
  cohort <- as.factor(cohort)
  if (nlevels(droplevels(cohort)) < 2L || length(cohort) != nrow(mat)) {
    return(mat)
  }
  out <- mat
  for (lv in levels(cohort)) {
    rows <- which(cohort == lv)
    if (length(rows) > 0L) {
      cm <- colMeans(mat[rows, , drop = FALSE])
      out[rows, ] <- sweep(mat[rows, , drop = FALSE], 2L, cm, `-`)
    }
  }
  out
}

#' Co-expression diversity of a single panel
#'
#' @param base_features Character vector of base feature names in the panel.
#' @param modules Named integer vector mapping feature -> module id.
#' @return Number of distinct modules spanned divided by panel size, in
#'   `(0, 1]`; `NA_real_` if no feature maps to a module or the panel is empty.
#' @noRd
.coexpr_diversity <- function(base_features, modules) {
  if (is.null(base_features) || length(base_features) == 0L) {
    return(NA_real_)
  }
  m <- modules[base_features]
  m <- m[!is.na(m)]
  if (length(m) == 0L) {
    return(NA_real_)
  }
  length(unique(m)) / length(base_features)
}
