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
                 "cohort", "healthy", "unaffected", "absent")

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
