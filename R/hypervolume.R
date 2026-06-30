#' @include panel_class.R
NULL

#' Hypervolume of an Optimization Result
#'
#' Compute the dominated hypervolume of the Pareto front contained in an
#' `OptimizationResult`. Hypervolume is a standard quality indicator for
#' multi-objective optimization: larger values indicate a front that
#' simultaneously improves on more objectives.
#'
#' Objective values are converted to a common minimization form using the
#' directions stored in `control$objective_directions` (i.e. `"maximize"`
#' objectives are negated). The reference point must therefore strictly dominate
#' every point in minimization form. By default the nadir of the observed front
#' is used and shifted outward by `ref_buffer * range` per objective so all
#' points contribute non-zero volume.
#'
#' Requires the suggested `emoa` package.
#'
#' @param object An `OptimizationResult`.
#' @param ref Optional numeric reference point in the **user-facing** direction
#'   (same scale as the columns in `solutions(object)`). If `NULL`, defaults to
#'   the front's nadir shifted by `ref_buffer`.
#' @param ref_buffer Fractional outward shift applied per objective when `ref`
#'   is auto-computed. Default `0.05`.
#' @return A single numeric value: the dominated hypervolume in minimization
#'   form. The reference point used is attached as `attr(., "ref")` (still in
#'   user-facing direction).
#' @seealso [hypervolume_history()] for HV trajectory across NSGA generations.
#' @export
setGeneric("hypervolume", function(object, ref = NULL, ref_buffer = 0.05)
  standardGeneric("hypervolume"))

#' @describeIn hypervolume Hypervolume of the final Pareto front.
#' @export
setMethod(
  "hypervolume",
  signature = "OptimizationResult",
  definition = function(object, ref = NULL, ref_buffer = 0.05) {
    directions <- object@control$objective_directions
    if (is.null(directions) || !length(directions)) {
      stop("OptimizationResult has no objective_directions in @control; ",
           "cannot compute hypervolume.", call. = FALSE)
    }

    sols <- object@solutions
    if (!nrow(sols)) {
      stop("OptimizationResult contains no solutions.", call. = FALSE)
    }
    if (!all(names(directions) %in% names(sols))) {
      stop("Objective columns missing from solutions data frame.", call. = FALSE)
    }

    points_user <- as.matrix(sols[, names(directions), drop = FALSE])
    .compute_hypervolume(points_user, directions, ref, ref_buffer)
  }
)

#' Hypervolume Trajectory Across NSGA Generations
#'
#' Compute the dominated hypervolume of the rank-1 (Pareto) population at each
#' recorded generation. Useful for diagnosing convergence — a flattening curve
#' suggests further iterations are unlikely to expand the front.
#'
#' Requires the optimization to have been run with `record_history = TRUE` and
#' the suggested `emoa` package to be installed. The reference point is held
#' fixed across generations so HV values are comparable: if not supplied, the
#' nadir of the *final* generation's Pareto set is used (with `ref_buffer`).
#'
#' @param object An `OptimizationResult` with non-empty `nsga_history()`.
#' @param ref Optional numeric reference point in user-facing direction. If
#'   `NULL`, derived from the final generation's Pareto set.
#' @param ref_buffer Fractional outward shift applied per objective when `ref`
#'   is auto-computed. Default `0.05`.
#' @param include_dominated Logical; if `TRUE`, compute HV over the full
#'   generation population rather than only the rank-1 individuals. Default
#'   `FALSE`.
#' @return A data frame with columns `generation` and `hypervolume`. The
#'   reference point used is attached as `attr(., "ref")`.
#' @export
setGeneric("hypervolume_history",
  function(object, ref = NULL, ref_buffer = 0.05, include_dominated = FALSE)
    standardGeneric("hypervolume_history"))

#' @describeIn hypervolume_history HV trajectory by generation.
#' @export
setMethod(
  "hypervolume_history",
  signature = "OptimizationResult",
  definition = function(object, ref = NULL, ref_buffer = 0.05,
                        include_dominated = FALSE) {
    hist <- nsga_history(object)
    if (!is.data.frame(hist) || !nrow(hist)) {
      stop("No NSGA history recorded. Re-run optimize_panel() with ",
           "record_history = TRUE.", call. = FALSE)
    }

    directions <- object@control$objective_directions
    obj_names <- names(directions)
    if (!all(obj_names %in% names(hist))) {
      stop("Objective columns missing from nsga_history().", call. = FALSE)
    }

    if (is.null(ref)) {
      final_gen <- max(hist$generation)
      final_pareto <- hist[hist$generation == final_gen & hist$is_pareto, obj_names,
                           drop = FALSE]
      if (!nrow(final_pareto)) {
        final_pareto <- hist[hist$generation == final_gen, obj_names, drop = FALSE]
      }
      ref <- .auto_reference(as.matrix(final_pareto), directions, ref_buffer)
    }

    gens <- sort(unique(hist$generation))
    hv <- vapply(gens, function(g) {
      rows <- if (isTRUE(include_dominated)) {
        hist$generation == g
      } else {
        hist$generation == g & hist$is_pareto
      }
      pts <- as.matrix(hist[rows, obj_names, drop = FALSE])
      if (!nrow(pts)) return(NA_real_)
      .compute_hypervolume(pts, directions, ref = ref, ref_buffer = ref_buffer)
    }, numeric(1))

    out <- data.frame(generation = gens, hypervolume = hv)
    attr(out, "ref") <- stats::setNames(ref, obj_names)
    out
  }
)

# Internal: convert user-facing points to minimization form, validate the
# reference point, and call emoa. `points_user` is a matrix with rows = points
# and cols = objectives in user-facing direction. `ref` (if supplied) is in the
# same direction; the auto-derived reference is also in that direction. The
# returned scalar is in minimization form (i.e. directly the volume that emoa
# computes — sign-flipping objectives does not change the hypervolume value).
.compute_hypervolume <- function(points_user, directions, ref, ref_buffer) {
  if (!requireNamespace("emoa", quietly = TRUE)) {
    stop("Package 'emoa' is required for hypervolume(). Install it with ",
         "install.packages('emoa').", call. = FALSE)
  }

  if (is.null(ref)) {
    ref <- .auto_reference(points_user, directions, ref_buffer)
  } else if (length(ref) != ncol(points_user)) {
    stop("`ref` must have length ", ncol(points_user), " (one per objective).",
         call. = FALSE)
  }

  pts_min <- .to_minimization(points_user, directions)
  ref_min <- .to_minimization(matrix(ref, nrow = 1L), directions)[1L, ]

  worse <- pts_min > matrix(ref_min, nrow = nrow(pts_min), ncol = length(ref_min),
                            byrow = TRUE)
  keep <- !apply(worse, 1L, any)
  if (!any(keep)) {
    out <- 0
    attr(out, "ref") <- stats::setNames(ref, names(directions))
    return(out)
  }
  pts_min <- pts_min[keep, , drop = FALSE]

  # emoa expects objectives as rows, points as columns.
  hv <- emoa::dominated_hypervolume(t(pts_min), ref_min)
  attr(hv, "ref") <- stats::setNames(ref, names(directions))
  hv
}

# Convert a points matrix from user-facing direction to minimization form by
# negating columns whose direction is "maximize".
.to_minimization <- function(points_user, directions) {
  out <- points_user
  for (j in seq_along(directions)) {
    if (directions[j] == "maximize") out[, j] <- -out[, j]
  }
  out
}

# Auto reference point in user-facing direction: take the nadir (worst observed
# value per objective) and shift outward by ref_buffer * (range of column).
# A constant column gets a fallback shift of `ref_buffer` (so the buffer is at
# least non-zero, e.g. for AUC plateaued at 1.0).
.auto_reference <- function(points_user, directions, ref_buffer) {
  ref <- numeric(ncol(points_user))
  for (j in seq_along(directions)) {
    col <- points_user[, j]
    rng <- diff(range(col))
    shift <- if (rng > 0) ref_buffer * rng else ref_buffer
    ref[j] <- if (directions[j] == "maximize") min(col) - shift else max(col) + shift
  }
  ref
}
