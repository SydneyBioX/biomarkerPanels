#' Batch-Associated Feature Selection
#'
#' Functions for identifying genes that capture batch effects but not biological
#' signal, using explicitly supplied cohort/batch labels.
#'
#' @name batch_associated_features
NULL

#' Select Batch-Associated Features
#'
#' Identifies genes that are associated with batch/cohort effects but NOT with
#' the biological outcome, using PCA plus per-component ANOVA against supplied
#' cohort labels. The returned genes act as proxies for the batch axis: they
#' are useful as a confounding diagnostic in their own right, and as log-ratio
#' denominators, where `log(Signal) - log(BatchProxy)` cancels batch-specific
#' scaling factors (see [construct_dual_stream_ratios()]).
#'
#' @param x Expression matrix (samples x genes). Should contain log-transformed
#'   expression values.
#' @param cohort Factor of cohort/batch labels, one per row of `x`.
#' @param y Binary outcome labels for biology filtering. Set to `NULL` to skip
#'   biology filtering (features selected purely on batch association).
#' @param n_features Number of batch-associated genes to select (default 20).
#' @param n_pcs Number of principal components to analyze (default 50).
#' @param batch_pvalue_threshold P-value threshold for batch association. PCs
#'   with ANOVA p < this threshold are considered batch-associated (default
#'   0.10).
#' @param biology_pvalue_threshold P-value threshold for biology non-association.
#'   PCs with biology p > this threshold are considered biology-free (default
#'   0.20).
#' @param biology_penalty_quantile Quantile of the per-gene biology scores used
#'   to stabilise the batch/biology ratio (default `0.5`, the median). This acts
#'   as a data-driven estimate of the biology-score noise floor; see Details.
#'   Larger values weaken the biology penalty, with `1.0` reducing the ranking
#'   to `batch_score` alone.
#'
#' @return A list with components:
#' \describe{
#'   \item{features}{Character vector of selected batch-associated gene names.}
#'   \item{biology_epsilon}{The stabilising constant actually added to the
#'     biology scores when computing `ratio`.}
#'   \item{diagnostics}{Data frame with per-PC statistics including variance
#'     explained, batch p-value/R^2, and biology p-value/R^2.
#'   }
#'   \item{problematic_cohort}{The cohort most associated with batch effects.}
#'   \item{one_vs_rest}{Data frame with per-cohort batch separation metrics.}
#'   \item{pure_batch_pcs}{Integer vector of PC indices identified as pure batch.}
#'   \item{gene_scores}{Data frame with per-gene batch and biology scores.}
#'   \item{pca}{The `prcomp` object from PCA on input data.}
#' }
#'
#' @details
#' The algorithm proceeds as follows:
#' 1. PCA on combined expression matrix
#' 2. Test each PC for batch association (ANOVA) and biology association (if
#'    `y` provided)
#' 3. Identify "pure batch" PCs: batch-significant AND biology-nonsignificant
#' 4. If no pure batch PCs found, fall back to LDA-based cohort discrimination
#' 5. Extract genes with highest loadings on pure batch PCs
#' 6. Score genes by ratio of batch to biology contribution
#' 7. Return top `n_features` genes
#'
#' Step 6 ranks by `batch_score / (biology_score + epsilon)`. Both scores are
#' variance-weighted means of absolute PCA loadings, so they share a scale that
#' depends on the number of genes and the data at hand. `epsilon` is therefore
#' set to the `biology_penalty_quantile` quantile of the observed biology
#' scores rather than a fixed constant: a fixed tiny epsilon lets genes whose
#' biology score is pure numerical noise (batch and biology both near zero)
#' generate arbitrarily large ratios and outrank genuinely batch-driven genes.
#' The median default assumes most genes are biologically null, making it a
#' robust estimate of the biology-score noise floor.
#'
#' Unlike [select_ruv_features()], which estimates a *latent* unwanted-variation
#' subspace from negative control genes, this function requires explicit
#' cohort/batch labels and targets the batch axis directly.
#'
#' @seealso [construct_dual_stream_ratios()] for using the selected genes as
#'   ratio denominators, [select_ruv_features()] for the label-free alternative.
#'
#' @examples
#' \dontrun{
#' # Assume x is samples x genes, cohort and y are factors
#' result <- select_batch_associated_features(
#'   x = expression_matrix,
#'   cohort = cohort_labels,
#'   y = disease_labels,
#'   n_features = 20
#' )
#'
#' # Use the batch-proxy genes as denominators
#' ratios <- construct_dual_stream_ratios(
#'   x = expression_matrix,
#'   numerators = signal_genes,
#'   denominators = result$features
#' )
#' }
#'
#' @export
#' @importFrom stats prcomp aov t.test var
#' @importFrom MASS lda
select_batch_associated_features <- function(
    x,
    cohort,
    y = NULL,
    n_features = 20L,
    n_pcs = 50L,
    batch_pvalue_threshold = 0.10,
    biology_pvalue_threshold = 0.20,
    biology_penalty_quantile = 0.5) {
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be a matrix or data.frame.", call. = FALSE)
  }
  x <- as.matrix(x)

  if (is.null(colnames(x))) {
    stop("`x` must have column names (gene identifiers).", call. = FALSE)
  }

  n_features <- .validate_positive_integer(n_features, "n_features")
  n_pcs <- .validate_positive_integer(n_pcs, "n_pcs")
  .validate_probability(biology_penalty_quantile, "biology_penalty_quantile",
    bounds = "closed"
  )

  cohort <- as.factor(cohort)
  if (length(cohort) != nrow(x)) {
    stop("`cohort` must have the same length as the number of rows in `x`.",
      call. = FALSE
    )
  }

  if (!is.null(y)) {
    y <- as.factor(y)
    if (length(y) != nrow(x)) {
      stop("`y` must have the same length as the number of rows in `x`.",
        call. = FALSE
      )
    }
  }

  # Step 1: PCA on combined expression matrix
  pca_result <- stats::prcomp(x, center = TRUE, scale. = FALSE)
  n_pcs <- min(n_pcs, ncol(pca_result$x))

  # Step 2: Per-PC diagnostics
  pc_diagnostics <- data.frame(
    PC = paste0("PC", seq_len(n_pcs)),
    PC_idx = seq_len(n_pcs),
    var_explained = (pca_result$sdev[seq_len(n_pcs)]^2) / sum(pca_result$sdev^2),
    batch_pvalue = NA_real_,
    batch_r2 = NA_real_,
    biology_pvalue = NA_real_,
    biology_r2 = NA_real_
  )

  for (i in seq_len(n_pcs)) {
    pc_scores <- pca_result$x[, i]

    # Test batch association (ANOVA)
    batch_anova <- summary(stats::aov(pc_scores ~ cohort))
    pc_diagnostics$batch_pvalue[i] <- batch_anova[[1]][["Pr(>F)"]][1]
    ss_batch <- batch_anova[[1]][["Sum Sq"]][1]
    ss_total <- sum(batch_anova[[1]][["Sum Sq"]])
    pc_diagnostics$batch_r2[i] <- ss_batch / ss_total

    # Test biology association (if outcome provided)
    if (!is.null(y)) {
      bio_anova <- summary(stats::aov(pc_scores ~ y))
      pc_diagnostics$biology_pvalue[i] <- bio_anova[[1]][["Pr(>F)"]][1]
      ss_bio <- bio_anova[[1]][["Sum Sq"]][1]
      ss_total_bio <- sum(bio_anova[[1]][["Sum Sq"]])
      pc_diagnostics$biology_r2[i] <- ss_bio / ss_total_bio
    }
  }

  # Step 3: Identify "pure batch" PCs
  if (!is.null(y)) {
    pure_batch_pcs <- which(
      pc_diagnostics$batch_pvalue < batch_pvalue_threshold &
        pc_diagnostics$biology_pvalue > biology_pvalue_threshold
    )
  } else {
    pure_batch_pcs <- which(pc_diagnostics$batch_pvalue < batch_pvalue_threshold)
  }

  # Fallback: if no pure batch PCs found, use ratio criterion
  if (length(pure_batch_pcs) == 0L && !is.null(y)) {
    message(
      "No strict pure-batch PCs found; using ratio criterion ",
      "(batch_r2 > 2*biology_r2)"
    )
    pure_batch_pcs <- which(
      pc_diagnostics$batch_r2 > (2 * pc_diagnostics$biology_r2) &
        pc_diagnostics$batch_pvalue < 0.20
    )
  }

  # If still empty, fall back to LDA-based approach
  if (length(pure_batch_pcs) == 0L) {
    message("Using LDA-based fallback to identify batch-associated PCs")
    cohort_df <- data.frame(cohort = cohort, pca_result$x[, seq_len(n_pcs)])
    var_criteria <- apply(cohort_df[, -1, drop = FALSE], 2, stats::var)
    cohort_df <- cohort_df[, c(TRUE, var_criteria > 1e-5), drop = FALSE]

    if (ncol(cohort_df) > 2L && nlevels(cohort) > 1L) {
      lda_fit <- MASS::lda(cohort ~ ., data = cohort_df)
      lda_importance <- rowSums(abs(lda_fit$scaling))
      top_lda_pcs <- names(sort(lda_importance, decreasing = TRUE))[
        seq_len(min(5L, length(lda_importance)))
      ]
      pure_batch_pcs <- as.integer(gsub("PC", "", top_lda_pcs))
    } else {
      # Last resort: use top variance PCs with any batch association
      pure_batch_pcs <- head(which(pc_diagnostics$batch_pvalue < 0.50), 3L)
    }
  }

  pc_diagnostics$is_pure_batch <- pc_diagnostics$PC_idx %in% pure_batch_pcs

  # Step 4: One-vs-rest analysis to identify problematic cohort
  cohort_levels <- levels(cohort)
  one_vs_rest <- data.frame(
    cohort = cohort_levels,
    mean_separation = NA_real_,
    n_associated_pcs = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(cohort_levels)) {
    is_target <- cohort == cohort_levels[i]
    n_sig <- 0L
    total_sep <- 0

    for (pc_idx in pure_batch_pcs) {
      pc_scores <- pca_result$x[, pc_idx]
      t_result <- stats::t.test(pc_scores[is_target], pc_scores[!is_target])
      if (t_result$p.value < 0.05) n_sig <- n_sig + 1L
      # Cohen's d for effect size
      d <- (mean(pc_scores[is_target]) - mean(pc_scores[!is_target])) /
        stats::sd(pc_scores)
      total_sep <- total_sep + abs(d)
    }

    one_vs_rest$mean_separation[i] <- total_sep / max(length(pure_batch_pcs), 1L)
    one_vs_rest$n_associated_pcs[i] <- n_sig
  }

  problematic_cohort <- cohort_levels[which.max(one_vs_rest$mean_separation)]

  # Step 5: Extract batch-associated genes from PCA loadings
  if (length(pure_batch_pcs) > 0L) {
    batch_loadings <- abs(pca_result$rotation[, pure_batch_pcs, drop = FALSE])

    # Weight by variance explained
    weights <- pc_diagnostics$var_explained[pure_batch_pcs]
    weights <- weights / sum(weights)

    # Weighted sum of absolute loadings
    gene_batch_scores <- as.vector(batch_loadings %*% weights)
    names(gene_batch_scores) <- colnames(x)
  } else {
    warning(
      "No batch-associated principal components were found after all fallbacks; ",
      "batch scores are all zero, so the returned features carry no batch ",
      "information and any ratios built from them will not cancel batch effects. ",
      "Inspect the data or relax the batch-PC criteria.",
      call. = FALSE
    )
    gene_batch_scores <- rep(0, ncol(x))
    names(gene_batch_scores) <- colnames(x)
  }

  # Step 6: Compute biology scores (to filter out)
  gene_biology_scores <- rep(0, ncol(x))
  names(gene_biology_scores) <- colnames(x)

  if (!is.null(y)) {
    bio_pcs <- which(pc_diagnostics$biology_pvalue < 0.10)

    if (length(bio_pcs) > 0L) {
      bio_loadings <- abs(pca_result$rotation[, bio_pcs, drop = FALSE])
      bio_weights <- pc_diagnostics$var_explained[bio_pcs]
      bio_weights <- bio_weights / sum(bio_weights)
      gene_biology_scores <- as.vector(bio_loadings %*% bio_weights)
      names(gene_biology_scores) <- colnames(x)
    }
  }

  # Stabilise the batch/biology ratio. A fixed tiny epsilon lets genes whose
  # biology score is pure numerical noise (both scores ~0) produce enormous
  # ratios and outrank genuinely batch-driven genes. Anchoring epsilon to a
  # quantile of the observed biology scores makes it scale-free: genes below
  # the biology noise floor are ordered by batch_score among themselves, while
  # genes above it are still penalised.
  biology_epsilon <- unname(stats::quantile(
    gene_biology_scores,
    probs = biology_penalty_quantile,
    na.rm = TRUE
  ))
  if (!is.finite(biology_epsilon) || biology_epsilon <= 0) {
    # No biology signal detected (y = NULL, no biology-associated PCs, or a
    # degenerate quantile). Fall back to a tiny constant, which makes `ratio`
    # monotone in batch_score - i.e. rank on batch association alone.
    biology_epsilon <- .Machine$double.eps
  }

  # Create gene score data frame
  gene_scores <- data.frame(
    gene = colnames(x),
    batch_score = gene_batch_scores,
    biology_score = gene_biology_scores,
    ratio = gene_batch_scores / (gene_biology_scores + biology_epsilon),
    stringsAsFactors = FALSE
  )

  # Select features: high batch score, low biology score
  gene_scores <- gene_scores[order(-gene_scores$ratio, -gene_scores$batch_score), ]

  features <- gene_scores$gene[seq_len(min(n_features, nrow(gene_scores)))]

  # Step 7: Validation - compute simple DE p-values for the selected features
  if (!is.null(y)) {
    levels_y <- levels(y)
    selected_de_pvalues <- vapply(features, function(g) {
      stats::t.test(
        x[y == levels_y[1], g],
        x[y == levels_y[2], g]
      )$p.value
    }, FUN.VALUE = numeric(1))

    gene_scores$de_pvalue <- NA_real_
    gene_scores$de_pvalue[match(features, gene_scores$gene)] <- selected_de_pvalues
  }

  list(
    features = features,
    biology_epsilon = biology_epsilon,
    diagnostics = pc_diagnostics,
    problematic_cohort = problematic_cohort,
    one_vs_rest = one_vs_rest,
    pure_batch_pcs = pure_batch_pcs,
    gene_scores = gene_scores,
    pca = pca_result
  )
}
