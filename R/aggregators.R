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
#' within-sample normalized differences.
#'
#' @details
#' For log-normalized expression data (which is standard for microarray and
#' RNA-seq), pairwise differences are mathematically equivalent to log-ratios.
#'
#' @param x A numeric matrix with column names. For best results, data should be
#'   log-normalized (e.g., log2 CPM, log2 RPKM, or microarray RMA values).
#' @return Matrix of pairwise differences with (p choose 2) columns.
#'   Column names follow the pattern "FeatureA--FeatureB".
#' @keywords internal
aggregator_pairwise_ratios <- function(x) {
  if (ncol(x) < 2L) {
    stop(
      "Aggregator 'pairwise_ratios' requires at least two features. ",
      "Input matrix has ", ncol(x), " column(s). ",
      "Use aggregator = 'none' if pairwise ratios are not needed.",
      call. = FALSE
    )
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
#'   Column names follow the pattern "FeatureA--FeatureB".
#' @keywords internal
aggregator_pairwise_log_ratios <- function(x) {
  if (ncol(x) < 2L) {
    stop(
      "Aggregator 'pairwise_log_ratios' requires at least two features. ",
      "Input matrix has ", ncol(x), " column(s). ",
      "Use aggregator = 'none' if pairwise log-ratios are not needed.",
      call. = FALSE
    )
  }
  if (is.null(colnames(x))) {
    stop("x must have column names", call. = FALSE)
  }

  # Check for non-positive values that would cause log issues
  if (any(x <= 0, na.rm = TRUE)) {
    n_nonpositive <- sum(x <= 0, na.rm = TRUE)
    stop(
      "Input matrix contains ", n_nonpositive, " non-positive value(s). ",
      "Log-ratios require strictly positive values. ",
      "Either add a pseudocount (x + 1), use 'pairwise_ratios' for log-transformed data, ",
      "or remove/impute non-positive values.",
      call. = FALSE
    )
  }

  # Apply log transformation first, then use pairwise_col_diff
  # Since log(A/B) = log(A) - log(B), this is mathematically equivalent
  pairwise_col_diff(log(x))
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

  ref_col <- x[, ref_feature]
  other_matrix <- x[, other_cols, drop = FALSE]

  .reference_diff_cpp(other_matrix, ref_col, ref_feature, other_cols)
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
