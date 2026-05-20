#' CPOP step 2 sign-agreement pruning
#'
#' Refits ridge logistic regression on the CPOP step-1 candidate set for each
#' cohort, then prunes features whose coefficient signs disagree across
#' cohorts. Iterated until the surviving set stabilises or all features are
#' eliminated, matching the original CPOP loop.
#'
#' @keywords internal
.cpop_step2 <- function(pair_matrices, responses, cpop1_features, n_iter,
                        lambda_choice, sign_threshold, family, cv_control) {
  remaining <- cpop1_features
  coef_matrix <- NULL

  if (length(remaining) < 2L) {
    return(list(
      features = remaining,
      coefficients = matrix(numeric(0), nrow = 0L,
                            ncol = length(pair_matrices))
    ))
  }

  for (j in seq_len(n_iter)) {
    coef_list <- vector("list", length(pair_matrices))
    for (k in seq_along(pair_matrices)) {
      cv_args <- utils::modifyList(list(
        x = pair_matrices[[k]][, remaining, drop = FALSE],
        y = if (family == "binomial") {
          as.integer(responses[[k]]) - 1L
        } else responses[[k]],
        alpha = 0,
        family = family,
        intercept = FALSE
      ), cv_control)
      cv_fit <- do.call(glmnet::cv.glmnet, cv_args)
      coef_mat <- as.matrix(stats::coef(cv_fit, s = lambda_choice))
      cv <- coef_mat[-1L, 1L]
      coef_list[[k]] <- cv[remaining]
    }
    coef_matrix <- do.call(cbind, coef_list)
    rownames(coef_matrix) <- remaining

    sign_mat <- sign(coef_matrix)
    agreement <- apply(sign_mat, 1L, function(s) {
      nz <- s[s != 0]
      if (!length(nz)) return(1.0)
      mean(nz == sign(stats::median(nz)))
    })
    keep <- agreement >= sign_threshold

    if (all(keep)) break

    remaining <- remaining[keep]
    if (!length(remaining)) break
  }

  list(features = remaining, coefficients = coef_matrix)
}

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
