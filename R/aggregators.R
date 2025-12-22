#' Cohort Aggregator Registry
#'
#' Functions for registering and retrieving cohort aggregation strategies
#' used to transform feature matrices before optimization.
#'
#' @name aggregator-registry
NULL

# -----------------------------------------------------------------------------
# Registry Environment
# -----------------------------------------------------------------------------

.aggregator_registry <- new.env(parent = emptyenv())

.register_default_aggregator <- function(name, fun, description) {

  assign(name, list(fun = fun, description = description),
         envir = .aggregator_registry)
}

#' Register a cohort aggregator function
#'
#' Add a custom aggregation strategy to the registry for use with
#' [optimize_panel()] and [evaluate_panel()].
#'
#' @param name Unique identifier for the aggregator.
#' @param fun Function implementing the aggregator. Must accept a matrix `x`
#'   with column names and return a transformed matrix.
#' @param description Human-readable description of the aggregator.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#'
#' @examples
#' # Register a custom centering aggregator
#' register_aggregator(
#'   "center",
#'   function(x) scale(x, center = TRUE, scale = FALSE),
#'   "Center each feature (subtract column means)"
#' )
#'
#' # List all registered aggregators
#' aggregator_registry()
#'
#' @seealso [aggregator_registry()], [optimize_panel()]
#' @export
register_aggregator <- function(name, fun, description = name, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  if (!overwrite && exists(name, envir = .aggregator_registry, inherits = FALSE)) {
    stop(sprintf("Aggregator '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(fun = fun, description = description),
         envir = .aggregator_registry)
  invisible(name)
}

#' List registered aggregator functions
#'
#' @return Named list of aggregator registrations (`fun`, `description`).
#'
#' @examples
#' aggregator_registry()
#'
#' @seealso [register_aggregator()]
#' @export
aggregator_registry <- function() {
  if (length(ls(.aggregator_registry)) == 0L) {
    return(list())
  }
  mget(ls(.aggregator_registry), envir = .aggregator_registry)
}

#' @keywords internal
.get_aggregator <- function(name) {
  if (!exists(name, envir = .aggregator_registry, inherits = FALSE)) {
    stop(
      "Aggregator '", name, "' is not registered. ",
      "Available: ", paste(ls(.aggregator_registry), collapse = ", "), ". ",
      "Use register_aggregator() to add custom aggregators.",
      call. = FALSE
    )
  }
  get(name, envir = .aggregator_registry, inherits = FALSE)
}

# -----------------------------------------------------------------------------
# Aggregator Implementations
# -----------------------------------------------------------------------------

#' Identity aggregator (no transformation)
#'
#' Returns the input matrix unchanged.
#'
#' @param x A numeric matrix with column names.
#' @return The input matrix unchanged.
#' @keywords internal
aggregator_none <- function(x) {
  x
}

#' Pairwise difference aggregator
#'
#' Computes pairwise differences between all feature pairs. This is the default
#' aggregator that helps dampen batch effects across cohorts by generating
#' within-sample contrasts.
#'
#' @param x A numeric matrix with column names.
#' @return Matrix of pairwise differences with (p choose 2) columns.
#'   Column names follow the pattern "FeatureA--FeatureB".
#' @keywords internal
aggregator_pairwise_ratios <- function(x) {
  if (ncol(x) < 2L) {
    warning(
      "Aggregator 'pairwise_ratios' requires at least two features; ",
      "returning original matrix.",
      call. = FALSE
    )
    return(x)
  }
  pairwise_col_diff(x)
}

#' Pairwise log-ratio aggregator
#'
#' Computes pairwise log-ratios between all feature pairs. Useful for
#' compositional data or when multiplicative relationships are expected.
#'
#' @param x A numeric matrix with column names. Values should be positive.
#' @return Matrix of log-ratios with (p choose 2) columns.
#'   Column names follow the pattern "FeatureA/FeatureB".
#' @keywords internal
aggregator_pairwise_log_ratios <- function(x) {
  if (ncol(x) < 2L) {
    warning(
      "Aggregator 'pairwise_log_ratios' requires at least two features; ",
      "returning original matrix.",
      call. = FALSE
    )
    return(x)
  }
  if (is.null(colnames(x))) {
    stop("x must have column names", call. = FALSE)
  }

  # Sort columns for consistency

col_order <- order(colnames(x))
  x <- x[, col_order, drop = FALSE]
  col_names <- colnames(x)

  n <- nrow(x)
  p <- ncol(x)
  n_pairs <- (p * (p - 1L)) / 2L

  result <- matrix(NA_real_, nrow = n, ncol = n_pairs)
  result_names <- character(n_pairs)

  idx <- 1L
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      ratio <- x[, i] / x[, j]
      result[, idx] <- log(ratio)
      result_names[idx] <- paste0(col_names[i], "/", col_names[j])
      idx <- idx + 1L
    }
  }

  colnames(result) <- result_names
  result
}

#' Reference normalization aggregator
#'
#' Normalizes all features relative to a specified reference feature (e.g., a
#' housekeeping gene). The reference feature must be specified via the
#' `"reference_feature"` attribute on the input matrix.
#'
#' @param x A numeric matrix with column names. Must have a `"reference_feature"`
#'   attribute specifying which column to use as the reference.
#' @return Matrix with (p - 1) columns, one for each non-reference feature.
#'   Column names follow the pattern "Feature--Reference".
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' colnames(x) <- c("GeneA", "GeneB", "Housekeeping")
#' attr(x, "reference_feature") <- "Housekeeping"
#' result <- aggregator_reference_norm(x)
#' }
#'
#' @keywords internal
aggregator_reference_norm <- function(x) {
  ref_feature <- attr(x, "reference_feature")
  if (is.null(ref_feature)) {
    stop(
      "Aggregator 'reference_norm' requires 'reference_feature' attribute. ",
      "Set via: attr(x, 'reference_feature') <- 'FeatureName'",
      call. = FALSE
    )
  }

  if (!ref_feature %in% colnames(x)) {
    stop(
      "Reference feature '", ref_feature, "' not found in data. ",
      "Available: ", paste(colnames(x), collapse = ", "),
      call. = FALSE
    )
  }

  ref_col <- x[, ref_feature]
  col_names <- colnames(x)
  other_cols <- setdiff(col_names, ref_feature)

  if (length(other_cols) == 0L) {
    warning(
      "Aggregator 'reference_norm' requires at least one non-reference feature; ",
      "returning original matrix.",
      call. = FALSE
    )
    return(x)
  }

  n <- nrow(x)
  p <- length(other_cols)

  result <- matrix(NA_real_, nrow = n, ncol = p)
  result_names <- character(p)

  for (i in seq_along(other_cols)) {
    feat_name <- other_cols[i]
    result[, i] <- x[, feat_name] - ref_col
    result_names[i] <- paste0(feat_name, "--", ref_feature)
  }

  colnames(result) <- result_names
  result
}

# -----------------------------------------------------------------------------
# Default Registrations
# -----------------------------------------------------------------------------

.register_default_aggregator(
  "none",
  aggregator_none,
"No transformation (identity)"
)

.register_default_aggregator(
  "pairwise_ratios",
  aggregator_pairwise_ratios,
  "Pairwise differences (A - B) for all feature pairs"
)

.register_default_aggregator(
  "pairwise_log_ratios",
  aggregator_pairwise_log_ratios,
  "Pairwise log-ratios log(A / B) for all feature pairs"
)

.register_default_aggregator(
  "reference_norm",
  aggregator_reference_norm,
  "Reference normalization (A - ref) for each feature relative to reference"
)
