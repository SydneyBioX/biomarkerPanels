#' Pareto solution utilities and selection strategies
#'
#' These helpers support downstream triage of Pareto-optimal solutions produced
#' by [`optimize_panel()`]. They implement the three post-hoc decision rules
#' used in the fetal growth restriction (FGR) analysis: maximising sensitivity
#' at a fixed specificity, prioritising biomarkers with high selection
#' frequency, and incorporating pathway-level biological priors.

#' Extract solution-level feature sets from a `BiomarkerPanelResult`.
#'
#' @param panel_result A [`BiomarkerPanelResult`].
#' @return A list where each element corresponds to a Pareto solution and
#'   contains the unique biomarker identifiers selected in that solution. The
#'   list is named by `solution_id`.
#' @keywords internal
.extract_solution_features <- function(panel_result) {
  stopifnot(inherits(panel_result, "BiomarkerPanelResult"))
  objectives <- panel_result@objectives
  if (!nrow(objectives)) {
    return(list())
  }

  split_idx <- split(seq_len(nrow(objectives)), objectives$solution_id)
  lapply(split_idx, function(idx) {
    feats <- objectives$features[idx]
    unique(unlist(feats, use.names = FALSE))
  })
}

#' Select the top-sensitivity panel at 90% specificity.
#'
#' Given solution-level performance summaries this function identifies the panel
#' that maximises external sensitivity while constraining the loss relative to
#' internal cross-validation performance.
#'
#' @param performance Data frame where each row summarises a solution. Must
#'   contain at minimum columns referenced by `solution_col`, `external_col`,
#'   and `internal_col`.
#' @param solution_col Column name holding the `solution_id`.
#' @param external_col Column name holding sensitivity at 90% specificity for
#'   the external validation cohort.
#' @param internal_col Column name holding the corresponding internal
#'   cross-validation sensitivity.
#' @param max_loss Maximum allowable drop in sensitivity between internal and
#'   external evaluations (`internal - external`). Defaults to `Inf`.
#' @param panel_size_col Optional column containing the number of biomarkers in
#'   each solution. Used purely as a tie-breaker (smaller panels preferred).
#' @return A single-row data frame describing the selected solution. An extra
#'   column named `performance_loss` (internal minus external sensitivity) is
#'   appended for convenience. Returns `NULL` if no solution satisfies
#'   `max_loss` or if required columns are missing.
#' @export
select_panel_top_sensitivity <- function(performance,
                                         solution_col = "solution_id",
                                         external_col = "external_sens_at_spec90",
                                         internal_col = "internal_sens_at_spec90",
                                         max_loss = Inf,
                                         panel_size_col = NULL) {
  required_cols <- c(solution_col, external_col, internal_col)
  missing_cols <- setdiff(required_cols, colnames(performance))
  if (length(missing_cols)) {
    warning("performance data is missing required column(s): ",
            paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  perf <- performance
  perf$performance_loss <- perf[[internal_col]] - perf[[external_col]]
  perf <- perf[is.finite(perf$performance_loss), , drop = FALSE]
  perf <- perf[perf$performance_loss <= max_loss, , drop = FALSE]
  perf <- perf[!is.na(perf[[external_col]]) & !is.na(perf[[internal_col]]), , drop = FALSE]

  if (!nrow(perf)) {
    return(NULL)
  }

  ordering <- order(
    -perf[[external_col]],
    abs(perf$performance_loss),
    -perf[[internal_col]],
    if (!is.null(panel_size_col) && panel_size_col %in% colnames(perf)) {
      perf[[panel_size_col]]
    } else {
      seq_len(nrow(perf))
    }
  )
  perf[ordering[1L], , drop = FALSE]
}

#' Compute biomarker inclusion frequencies across optimisation runs.
#'
#' @param panels List of [`BiomarkerPanelResult`] objects, character vectors of
#'   biomarkers, or a mixture of both.
#' @return Data frame with columns `feature`, `count`, and `proportion`.
#' @export
compute_inclusion_frequencies <- function(panels) {
  if (is.null(panels) || !length(panels)) {
    return(data.frame(
      feature = character(),
      count = integer(),
      proportion = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  collected <- lapply(panels, function(p) {
    if (inherits(p, "BiomarkerPanelResult")) {
      sol <- .extract_solution_features(p)
      sol <- lapply(sol, unique)
      list(
        features = unlist(sol, use.names = FALSE),
        solutions = length(sol)
      )
    } else if (is.character(p)) {
      list(
        features = unique(p),
        solutions = 1L
      )
    } else {
      stop("Unsupported panel entry: must be character vector or BiomarkerPanelResult.",
           call. = FALSE)
    }
  })

  flattened <- unlist(lapply(collected, `[[`, "features"), use.names = FALSE)

  if (!length(flattened)) {
    return(data.frame(
      feature = character(),
      count = integer(),
      proportion = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  freq_table <- sort(table(flattened), decreasing = TRUE)
  total_solutions <- sum(vapply(collected, `[[`, integer(1), "solutions"))
  df <- data.frame(
    feature = names(freq_table),
    count = as.integer(freq_table),
    stringsAsFactors = FALSE
  )
  df$proportion <- df$count / total_solutions
  df
}

#' Assemble a panel by inclusion frequency.
#'
#' @param frequency_df Data frame returned by
#'   [compute_inclusion_frequencies()].
#' @param panel_size Target number of biomarkers.
#' @return Character vector of selected biomarker identifiers (length up to
#'   `panel_size`).
#' @export
select_panel_inclusion_frequency <- function(frequency_df, panel_size = 4L) {
  if (!nrow(frequency_df)) {
    return(character())
  }
  panel_size <- .validate_positive_integer(panel_size, "panel_size")
  head(frequency_df$feature, n = min(panel_size, nrow(frequency_df)))
}

#' Select biomarkers guided by pathway information.
#'
#' @param high_sensitivity_panels List of character vectors (or
#'   `BiomarkerPanelResult` objects) representing high-performing panels.
#' @param inclusion_frequency Data frame from [compute_inclusion_frequencies()].
#' @param pathway_mapping Data frame with at least columns `feature` and
#'   `pathway`. Optional columns `pathway_name` and `source` provide additional
#'   context.
#' @param relevant_pathways Character vector of pathway identifiers (matched
#'   against `pathway_mapping$pathway`) that are known to be implicated in the
#'   disease under study.
#' @param feature_performance Data frame containing per-feature performance
#'   summaries. Must include columns `feature` and the column referenced by
#'   `sensitivity_col`.
#' @param panel_size Desired final panel size (default `4`).
#' @param min_sensitivity Minimum acceptable sensitivity at 90% specificity for
#'   biomarkers selected via pathway enrichment.
#' @param sensitivity_col Column name in `feature_performance` holding the
#'   sensitivity estimates (defaults to `"external_sens_at_spec90"`).
#' @param top_frequency Optional number of highly frequent biomarkers to add to
#'   the candidate pool (defaults to `panel_size`).
#' @return A list containing `features` (selected biomarker identifiers) and
#'   `pathway_assignments` describing the mapping between selected biomarkers
#'   and pathways. The returned vector may be shorter than `panel_size` if the
#'   constraints cannot be satisfied.
#' @export
select_panel_by_pathway <- function(high_sensitivity_panels,
                                    inclusion_frequency,
                                    pathway_mapping,
                                    relevant_pathways,
                                    feature_performance,
                                    panel_size = 4L,
                                    min_sensitivity = 0.5,
                                    sensitivity_col = "external_sens_at_spec90",
                                    top_frequency = panel_size) {
  panel_size <- .validate_positive_integer(panel_size, "panel_size")
  if (!is.data.frame(pathway_mapping) || !nrow(pathway_mapping)) {
    stop("`pathway_mapping` must be a non-empty data frame.", call. = FALSE)
  }
  required_pm <- c("feature", "pathway")
  missing_pm <- setdiff(required_pm, colnames(pathway_mapping))
  if (length(missing_pm)) {
    stop("`pathway_mapping` is missing column(s): ",
         paste(missing_pm, collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(feature_performance) || !nrow(feature_performance)) {
    stop("`feature_performance` must be a non-empty data frame.", call. = FALSE)
  }
  if (!sensitivity_col %in% colnames(feature_performance)) {
    stop("`feature_performance` is missing the column ", sQuote(sensitivity_col),
         call. = FALSE)
  }

  hs_features <- lapply(high_sensitivity_panels, function(p) {
    if (inherits(p, "BiomarkerPanelResult")) {
      unique(unlist(.extract_solution_features(p), use.names = FALSE))
    } else if (is.character(p)) {
      p
    } else {
      stop("Unsupported entry in `high_sensitivity_panels`. Use character vectors or BiomarkerPanelResult objects.",
           call. = FALSE)
    }
  })
  hs_features <- unique(unlist(hs_features, use.names = FALSE))

  freq_candidates <- character()
  if (!missing(inclusion_frequency) && nrow(inclusion_frequency)) {
    if (is.null(top_frequency) || !is.finite(top_frequency)) {
      top_frequency <- panel_size
    }
    top_frequency <- max(panel_size, as.integer(top_frequency))
    freq_candidates <- head(inclusion_frequency$feature, n = min(top_frequency, nrow(inclusion_frequency)))
  }

  candidate_features <- unique(c(hs_features, freq_candidates))
  mapping_subset <- pathway_mapping[pathway_mapping$feature %in% candidate_features,
                                    , drop = FALSE]
  if (!length(relevant_pathways)) {
    relevant_pathways <- unique(mapping_subset$pathway)
  }
  mapping_subset <- mapping_subset[mapping_subset$pathway %in% relevant_pathways,
                                   , drop = FALSE]
  if (!nrow(mapping_subset)) {
    return(list(
      features = character(),
      pathway_assignments = data.frame()
    ))
  }

  perf_subset <- feature_performance[feature_performance$feature %in% candidate_features,
                                     c("feature", sensitivity_col),
                                     drop = FALSE]
  perf_subset <- perf_subset[!is.na(perf_subset[[sensitivity_col]]), , drop = FALSE]

  if (!nrow(perf_subset)) {
    return(list(
      features = character(),
      pathway_assignments = data.frame()
    ))
  }

  selection <- character()
  pathway_assignments <- list()
  for (path in unique(mapping_subset$pathway)) {
    feats <- mapping_subset$feature[mapping_subset$pathway == path]
    feats <- intersect(feats, perf_subset$feature)
    if (!length(feats)) {
      next
    }
    perf_rows <- perf_subset[match(feats, perf_subset$feature), , drop = FALSE]
    perf_rows <- perf_rows[perf_rows[[sensitivity_col]] >= min_sensitivity, , drop = FALSE]
    if (!nrow(perf_rows)) {
      next
    }
    best_idx <- which.max(perf_rows[[sensitivity_col]])
    chosen <- perf_rows$feature[best_idx]
    if (!chosen %in% selection) {
      selection <- c(selection, chosen)
      pathway_assignments[[chosen]] <- data.frame(
        feature = chosen,
        pathway = path,
        stringsAsFactors = FALSE
      )
    }
    if (length(selection) >= panel_size) {
      break
    }
  }

  if (length(selection) < panel_size) {
    remaining <- setdiff(perf_subset$feature[perf_subset[[sensitivity_col]] >= min_sensitivity],
                         selection)
    if (length(remaining)) {
      ordering <- order(perf_subset[[sensitivity_col]][match(remaining, perf_subset$feature)],
                        decreasing = TRUE)
      fill <- remaining[ordering]
      selection <- unique(c(selection, head(fill, panel_size - length(selection))))
    }
  }

  assignments_df <- if (length(pathway_assignments)) {
    do.call(rbind, pathway_assignments)
  } else {
    data.frame()
  }

  list(
    features = selection,
    pathway_assignments = assignments_df
  )
}
