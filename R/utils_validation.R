#' Validate Binary Response
#'
#' Ensures that the response vector conforms to the required binary labels.
#' 
#' @param y Vector containing two outcome levels. Accepted types include factor,
#'   character, logical, or numeric. Level names are coerced to the canonical
#'   `c("No", "Yes")` order.
#' @param positive Optional character vector of value(s) to treat as the
#'   positive (`"Yes"`) class. Comparisons are case-insensitive.
#' @param negative Optional character vector of value(s) to treat as the
#'   negative (`"No"`) class. Comparisons are case-insensitive.
#' @return Factor with levels `c("No", "Yes")`.
#' @export
#'
#' @examples
#' ensure_binary_response(c("control", "fgr", "control"))
#' ensure_binary_response(factor(c("case", "control")), positive = "case")
ensure_binary_response <- function(y, positive = NULL, negative = NULL) {
  if (is.null(y)) {
    stop("`y` cannot be NULL.", call. = FALSE)
  }

  # Short-circuit: if y is already a factor with levels c("No", "Yes"),
  # return it unchanged. This handles subsetted factors from cohort-aware
  # metric functions where only one class may be present in the data.
  if (is.factor(y) && identical(levels(y), c("No", "Yes"))) {
    return(y)
  }

  # Convert supported input types to a character vector for processing.
  if (is.factor(y)) {
    y_chr <- as.character(y)
  } else if (is.logical(y)) {
    y_chr <- ifelse(is.na(y), NA_character_, ifelse(y, "Yes", "No"))
  } else if (is.numeric(y)) {
    y_chr <- ifelse(is.na(y), NA_character_, as.character(y))
  } else {
    y_chr <- as.character(y)
  }

  tokens <- tolower(trimws(y_chr))
  unique_tokens <- unique(tokens[!is.na(tokens)])

  if (length(unique_tokens) > 2) {
    stop("`y` must contain exactly two outcome values.", call. = FALSE)
  }

  yes_tokens <- c("yes", "y", "1", "true", "positive", "pos", "case",
                  "fgr", "disease", "affected", "present")
  no_tokens <- c("no", "n", "0", "false", "negative", "neg", "control",
                 "healthy", "unaffected", "absent")

  if (!is.null(positive)) {
    positive <- tolower(trimws(as.character(positive)))
  } else {
    positive <- unique_tokens[unique_tokens %in% yes_tokens]
  }

  if (!is.null(negative)) {
    negative <- tolower(trimws(as.character(negative)))
  } else {
    negative <- unique_tokens[unique_tokens %in% no_tokens]
  }

  # Fall back to assigning remaining levels when heuristics are inconclusive.
  unresolved <- setdiff(unique_tokens, c(positive, negative))
  used_fallback <- FALSE
  if (length(unresolved)) {
    if (length(unresolved) > 2) {
      stop("Unable to determine binary class mapping for values: ",
           paste(sort(unique(trimws(y_chr[ tokens %in% unresolved ]))), collapse = ", "),
           call. = FALSE)
    }
    if (!length(negative)) {
      negative <- unresolved[1]
      unresolved <- setdiff(unresolved, negative)
      used_fallback <- TRUE
    }
    if (length(unresolved) && !length(positive)) {
      positive <- unresolved[1]
      unresolved <- setdiff(unresolved, positive)
      used_fallback <- TRUE
    }
  }

  if (!length(positive) || !length(negative)) {
    stop(
      "Unable to infer which level represents the positive versus negative class. ",
      "Provide explicit `positive` and `negative` mappings.",
      call. = FALSE
    )
  }

  mapped <- rep(NA_character_, length(tokens))
  mapped[tokens %in% positive] <- "Yes"
  mapped[tokens %in% negative] <- "No"

  if (used_fallback) {
    warning(
      "Unable to infer positive/negative classes from values. ",
      "Coercing ", sQuote(negative[1]), " to 'No' and ",
      sQuote(positive[1]), " to 'Yes'. ",
      "Supply `positive`/`negative` to override.",
      call. = FALSE
    )
  }

  unmapped <- which(!is.na(tokens) & is.na(mapped))
  if (length(unmapped)) {
    stop("Unable to map the following response values: ",
         paste(sort(unique(y_chr[unmapped])), collapse = ", "),
         ". Supply `positive`/`negative` arguments to resolve this.",
         call. = FALSE)
  }

  factor(mapped, levels = c("No", "Yes"))
}

#' Ensure Input is a List of Cohorts
#'
#' Wraps single objects in a list for consistent multi-cohort handling.
#'
#' @param x An object or list of objects.
#' @return A list.
#' @keywords internal
.as_cohort_list <- function(x) { # nocov start
  if (is.null(x)) {
    stop("Input cannot be NULL.", call. = FALSE)
  }

  if (is.list(x)) {
    return(x)
  }
  list(x)
} # nocov end

#' Validate Positive Integer Parameter
#'
#' Checks that a parameter is a single integer no smaller than `min` and coerces
#' it. With the default `min = 1L` this enforces a positive integer; pass a
#' larger `min` (e.g. `2L` for a fold or split count) to require a higher floor.
#'
#' @param x The value to validate.
#' @param name Parameter name for error messages.
#' @param min Minimum allowed value (default `1L`).
#' @return An integer value.
#' @keywords internal
.validate_positive_integer <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min) {
    msg <- if (min <= 1L) {
      paste0("`", name, "` must be a positive integer.")
    } else {
      paste0("`", name, "` must be an integer >= ", min, ".")
    }
    stop(msg, call. = FALSE)
  }
  as.integer(x)
}

#' Validate a Numeric Scalar in the Unit Interval
#'
#' Checks that a parameter is a single non-missing numeric value inside the unit
#' interval. Use `bounds = "open"` for a probability-like value in (0, 1) and
#' `bounds = "closed"` for a value in [0, 1] where the endpoints are allowed
#' (e.g. an elastic-net mixing parameter or a proportion).
#'
#' @param x The value to validate.
#' @param name Parameter name for error messages.
#' @param bounds Either `"open"` (0 and 1 excluded) or `"closed"` (0 and 1
#'   allowed).
#' @return The value coerced with [as.numeric()].
#' @keywords internal
.validate_probability <- function(x, name, bounds = c("open", "closed")) {
  bounds <- match.arg(bounds)
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x)
  if (ok) {
    ok <- if (bounds == "open") x > 0 && x < 1 else x >= 0 && x <= 1
  }
  if (!ok) {
    stop("`", name, "` must be a numeric scalar between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

#' Validate a Numeric Scalar
#'
#' Checks that a parameter is a single numeric value.
#'
#' @param x The value to validate.
#' @param name Parameter name for error messages.
#' @return The value coerced with [as.numeric()].
#' @keywords internal
.validate_numeric_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop("`", name, "` must be a numeric scalar.", call. = FALSE)
  }
  as.numeric(x)
}
