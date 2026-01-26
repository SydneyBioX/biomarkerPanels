#' Feature Transform Registry
#'
#' Functions for registering and retrieving feature transformation strategies
#' used to transform feature matrices during optimization.
#'
#' @name transform-registry
NULL

# -----------------------------------------------------------------------------
# Registry Environment
# -----------------------------------------------------------------------------

.transform_registry <- new.env(parent = emptyenv())

.register_default_transform <- function(name, fun, description) {

  assign(name, list(fun = fun, description = description),
         envir = .transform_registry)
}

#' Register a feature transform function
#'
#' Add a custom transformation strategy to the registry for use with
#' [optimize_panel()] and [evaluate_panel()].
#'
#' @param name Unique identifier for the transform.
#' @param fun Function implementing the transform. Must accept a matrix `x`
#'   with column names and return a transformed matrix.
#' @param description Human-readable description of the transform.
#' @param overwrite Logical; set to `TRUE` to replace an existing registration.
#' @return Invisibly, the registered name.
#'
#' @examples
#' # Register a custom centering transform
#' register_feature_transform(
#'   "center",
#'   function(x) scale(x, center = TRUE, scale = FALSE),
#'   "Center each feature (subtract column means)"
#' )
#'
#' # List all registered transforms
#' feature_transform_registry()
#'
#' @seealso [feature_transform_registry()], [optimize_panel()]
#' @export
register_feature_transform <- function(name, fun, description = name, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.function(fun))
  if (!overwrite && exists(name, envir = .transform_registry, inherits = FALSE)) {
    stop(sprintf("Feature transform '%s' is already registered.", name), call. = FALSE)
  }
  assign(name, list(fun = fun, description = description),
         envir = .transform_registry)
  invisible(name)
}

#' List registered feature transform functions
#'
#' @return Named list of transform registrations (`fun`, `description`).
#'
#' @examples
#' feature_transform_registry()
#'
#' @seealso [register_feature_transform()]
#' @export
feature_transform_registry <- function() {
  if (length(ls(.transform_registry)) == 0L) {
    return(list())
  }
  mget(ls(.transform_registry), envir = .transform_registry)
}

#' @keywords internal
.get_feature_transform <- function(name) {
  if (!exists(name, envir = .transform_registry, inherits = FALSE)) {
    stop(
      "Feature transform '", name, "' is not registered. ",
      "Available: ", paste(ls(.transform_registry), collapse = ", "), ". ",
      "Use register_feature_transform() to add custom transforms.",
      call. = FALSE
    )
  }
  get(name, envir = .transform_registry, inherits = FALSE)
}

# -----------------------------------------------------------------------------
# Transform Implementations
# -----------------------------------------------------------------------------

#' Identity transform (no transformation)
#'
#' Returns the input matrix unchanged.
#'
#' @param x A numeric matrix with column names.
#' @return The input matrix unchanged.
#' @keywords internal
transform_none <- function(x) {
  x
}

#' Pairwise difference transform
#'
#' Computes pairwise differences between all feature pairs. This is the default
#' transform that helps dampen batch effects across cohorts by generating
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
transform_pairwise_ratios <- function(x) {
  if (ncol(x) < 2L) {
    stop(
      "Transform 'pairwise_ratios' requires at least two features. ",
      "Input matrix has ", ncol(x), " column(s). ",
      "Use feature_transform = 'none' if pairwise ratios are not needed.",
      call. = FALSE
    )
  }
  pairwise_col_diff(x)
}

#' Pairwise log-ratio transform
#'
#' Computes pairwise log-ratios between all feature pairs. Useful for
#' compositional data or when multiplicative relationships are expected.
#'
#' @param x A numeric matrix with column names. Values should be positive.
#' @return Matrix of log-ratios with (p choose 2) columns.
#'   Column names follow the pattern "FeatureA--FeatureB".
#' @keywords internal
transform_pairwise_log_ratios <- function(x) {
  if (ncol(x) < 2L) {
    stop(
      "Transform 'pairwise_log_ratios' requires at least two features. ",
      "Input matrix has ", ncol(x), " column(s). ",
      "Use feature_transform = 'none' if pairwise log-ratios are not needed.",
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

#' Reference normalization transform
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
#' result <- transform_reference_norm(x)
#' }
#'
#' @keywords internal
transform_reference_norm <- function(x) {
  ref_feature <- attr(x, "reference_feature")
  if (is.null(ref_feature)) {
    stop(
      "Transform 'reference_norm' requires 'reference_feature' attribute. ",
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
      "Transform 'reference_norm' requires at least one non-reference feature; ",
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

.register_default_transform(
  "none",
  transform_none,
"No transformation (identity)"
)

.register_default_transform(
  "pairwise_ratios",
  transform_pairwise_ratios,
  "Pairwise differences (A - B) for all feature pairs"
)

.register_default_transform(
  "pairwise_log_ratios",
  transform_pairwise_log_ratios,
  "Pairwise log-ratios log(A / B) for all feature pairs"
)

.register_default_transform(
  "reference_norm",
  transform_reference_norm,
  "Reference normalization (A - ref) for each feature relative to reference"
)
