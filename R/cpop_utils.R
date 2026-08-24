#' CPOP Shared Utilities
#'
#' Internal helpers shared across the CPOP feature-selection workflow
#' ([select_cpop_features()], [fit_cpop_panel()]): input preparation, penalty
#' normalization, and MST-based pruning of redundant ratio features.
#'
#' @name cpop_utils
NULL

#' MST pruning of CPOP ratio features
#'
#' Removes redundant pair features by retaining a spanning tree of the gene
#' graph implied by the `"A--B"` edges. Equivalent to CPOP's `mst_lratio` but
#' implemented with a base-R union-find Kruskal pass — no `igraph` dependency.
#' All edges share equal weight; ties break by input order.
#'
#' @keywords internal
.mst_lratio_pruning <- function(pair_features) {
  if (!length(pair_features)) return(pair_features)

  parts <- strsplit(pair_features, "--", fixed = TRUE)
  if (any(lengths(parts) != 2L)) {
    stop("CPOP pair feature names must use the 'A--B' separator.",
         call. = FALSE)
  }
  edges <- do.call(rbind, parts)
  nodes <- unique(c(edges))
  node_idx <- stats::setNames(seq_along(nodes), nodes)
  parent <- seq_along(nodes)

  find_root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }

  keep <- logical(nrow(edges))
  for (i in seq_len(nrow(edges))) {
    a <- find_root(node_idx[[edges[i, 1L]]])
    b <- find_root(node_idx[[edges[i, 2L]]])
    if (a != b) {
      parent[a] <- b
      keep[i] <- TRUE
    }
  }

  ekept <- edges[keep, , drop = FALSE]
  sorted_pairs <- t(apply(ekept, 1L, sort))
  paste0(sorted_pairs[, 1L], "--", sorted_pairs[, 2L])
}

.p_transform <- function(penalty) {
  p <- length(penalty)
  pen <- pmax(penalty, 0)
  if (sum(pen) > 0) {
    return(pen * p / sum(pen))
  }
  rep(1, p)
}

.cpop_prepare_inputs <- function(x_list, y_list, assay = NULL) {
  if (length(.as_cohort_list(x_list)) < 2L) {
    stop("CPOP requires at least two cohorts.", call. = FALSE)
  }

  prepared <- .prepare_selection_inputs(
    x_list, y_list,
    assay = assay,
    response = "factor",
    feature_order = "sorted"
  )
  raw_matrices <- prepared$matrices
  responses <- prepared$responses
  cohort_names <- prepared$cohort_names
  ordered_features <- prepared$feature_names

  if (length(ordered_features) < 2L) {
    stop("CPOP requires at least two shared base features to construct ",
         "pairwise ratios.", call. = FALSE)
  }

  pair_matrices <- lapply(raw_matrices, pairwise_col_diff)
  pair_features <- colnames(pair_matrices[[1]])

  list(
    pair_matrices = pair_matrices,
    responses = responses,
    cohort_names = cohort_names,
    pair_features = pair_features
  )
}
