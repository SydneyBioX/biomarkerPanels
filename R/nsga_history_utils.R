#' NSGA Per-Generation History Utilities
#'
#' Internal helpers for capturing and materializing the per-generation
#' population/fitness/front history during an NSGA-II/III run. The capture
#' machinery is shared verbatim by [optimize_panel()] and
#' [optimize_panel_transferable()]; these factories remove that duplication.
#'
#' @name nsga_history_utils
#' @keywords internal
NULL

#' Create an Empty History Buffer
#'
#' Allocates a parent-less environment used to accumulate per-generation
#' snapshots without copying on each push.
#'
#' @return An environment with an empty `gens` list.
#' @keywords internal
.make_history_buffer <- function() {
  buffer <- new.env(parent = emptyenv())
  buffer$gens <- list()
  buffer
}

#' Create a Per-Generation History Monitor
#'
#' Builds the `monitor` closure passed to `rmoo::nsga2()`/`nsga3()`. Each
#' generation it pushes the population, sign-corrected fitness, and front rank
#' into `buffer`. Decision-weight vectors are stored raw; base features are
#' reconstructed later in [.materialize_history()] to keep per-iteration
#' overhead low.
#'
#' @param buffer A history buffer from [.make_history_buffer()].
#' @param objective_directions Character vector of per-objective directions
#'   (`"maximize"`/`"minimize"`); maximized columns are sign-flipped back to
#'   their natural scale.
#' @param objective_names Column names to assign to the fitness matrix.
#' @return A function with signature `function(object, ...)`.
#' @keywords internal
.make_history_monitor <- function(buffer, objective_directions, objective_names) {
  function(object, ...) {
    iter <- object@iter
    fit <- object@fitness
    pop <- object@population
    if (is.null(dim(fit))) fit <- matrix(fit, nrow = 1L)
    if (is.null(dim(pop))) pop <- matrix(pop, nrow = 1L)
    fr <- as.integer(object@front)
    for (j in seq_along(objective_directions)) {
      if (objective_directions[j] == "maximize") fit[, j] <- -fit[, j]
    }
    colnames(fit) <- objective_names
    buffer$gens[[length(buffer$gens) + 1L]] <- list(
      iter = as.integer(iter),
      fit = fit,
      front = fr,
      pop = pop
    )
    invisible(NULL)
  }
}

#' Materialize Captured History into a Data Frame
#'
#' Converts the raw per-generation snapshots in `buffer` into a long data frame,
#' running `select_base_features()` on each individual so the recorded panels
#' exactly match what NSGA was scoring.
#'
#' @param buffer A history buffer populated by the monitor.
#' @param select_base_features Function mapping a decision vector to the
#'   selected base features (the same selection logic the optimizer used).
#' @return A data frame with one row per individual per generation, or an empty
#'   list when no history was captured.
#' @keywords internal
.materialize_history <- function(buffer, select_base_features) {
  if (!length(buffer$gens)) {
    return(list())
  }
  gen_dfs <- lapply(buffer$gens, function(g) {
    n_ind <- nrow(g$pop)
    bf <- lapply(seq_len(n_ind), function(i) select_base_features(g$pop[i, ]))
    df <- data.frame(
      generation = g$iter,
      individual = seq_len(n_ind),
      rank = g$front,
      is_pareto = g$front == 1L,
      g$fit,
      n_features = lengths(bf),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    df$base_features <- bf
    df
  })
  do.call(rbind, gen_dfs)
}
