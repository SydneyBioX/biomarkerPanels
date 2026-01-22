#' Data Partitioning for Cross-Validation
#'
#' Internal functions for partitioning multi-cohort data into train, validation,
#' and held-out sets while maintaining stratified class balance.
#'
#' @name data_partitioning
#' @keywords internal
NULL

#' Validate Partition Ratios for Train/Validation/Held-out Split
#'
#' Ensures that train, validation, and held-out ratios satisfy minimum
#' requirements for meaningful model fitting and evaluation.
#'
#' @param train_ratio Proportion of data for training (must be >= 0.5).
#' @param val_ratio Proportion of data for validation (must be >= 0.1).
#' @return Invisibly returns TRUE if valid; otherwise throws an error.
#' @keywords internal
.validate_partition_ratios <- function(train_ratio, val_ratio) {
  if (!is.numeric(train_ratio) || length(train_ratio) != 1L ||
      is.na(train_ratio) || train_ratio < 0 || train_ratio > 1) {
    stop("`train_ratio` must be a single numeric value between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(val_ratio) || length(val_ratio) != 1L ||
      is.na(val_ratio) || val_ratio < 0 || val_ratio > 1) {
    stop("`val_ratio` must be a single numeric value between 0 and 1.",
         call. = FALSE)
  }

  heldout_ratio <- 1 - train_ratio - val_ratio

  if (train_ratio < 0.5) {
    stop("`train_ratio` must be at least 0.5 to ensure sufficient training data.",
         call. = FALSE)
  }
  if (val_ratio < 0.1) {
    stop("`val_ratio` must be at least 0.1 to ensure meaningful validation.",
         call. = FALSE)
  }
  if (heldout_ratio < 0.05) {
    stop("Held-out ratio (1 - train_ratio - val_ratio) must be at least 0.05. ",
         "Current: ", round(heldout_ratio, 3), call. = FALSE)
  }

  invisible(TRUE)
}

#' Stratified Partitioning of Multi-Cohort Data
#'
#' Partitions each cohort's data into train, validation, and held-out sets
#' while maintaining class balance within each partition.
#'
#' @param x_list List of feature matrices (one per cohort).
#' @param y_list List of binary response factors (one per cohort).
#' @param train_ratio Proportion of data for training.
#' @param val_ratio Proportion of data for validation.
#' @return List containing:
#'   \describe{
#'     \item{train_x}{List of training matrices (one per cohort)}
#'     \item{train_y}{List of training response factors}
#'     \item{val_x}{List of validation matrices}
#'     \item{val_y}{List of validation response factors}
#'     \item{heldout_x}{List of held-out matrices}
#'     \item{heldout_y}{List of held-out response factors}
#'     \item{cohort_names}{Character vector of cohort names}
#'     \item{partition_info}{Data frame with partition sizes per cohort}
#'   }
#' @keywords internal
.stratified_partition_cohorts <- function(x_list, y_list, train_ratio, val_ratio) {
  if (!is.list(x_list) || !is.list(y_list)) {
    stop("`x_list` and `y_list` must be lists.", call. = FALSE)
  }
  if (length(x_list) != length(y_list)) {
    stop("`x_list` and `y_list` must have the same length.", call. = FALSE)
  }

  n_cohorts <- length(x_list)
  cohort_names <- names(x_list)
  if (is.null(cohort_names) || any(cohort_names == "")) {
    cohort_names <- sprintf("cohort_%02d", seq_len(n_cohorts))
  }

  train_x <- vector("list", n_cohorts)
  train_y <- vector("list", n_cohorts)
  val_x <- vector("list", n_cohorts)
  val_y <- vector("list", n_cohorts)
  heldout_x <- vector("list", n_cohorts)
  heldout_y <- vector("list", n_cohorts)

  partition_info <- data.frame(
    cohort = cohort_names,
    n_total = integer(n_cohorts),
    n_train = integer(n_cohorts),
    n_val = integer(n_cohorts),
    n_heldout = integer(n_cohorts),
    train_yes = integer(n_cohorts),
    train_no = integer(n_cohorts),
    val_yes = integer(n_cohorts),
    val_no = integer(n_cohorts),
    heldout_yes = integer(n_cohorts),
    heldout_no = integer(n_cohorts),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_cohorts)) {
    x_i <- x_list[[i]]
    y_i <- ensure_binary_response(y_list[[i]])

    n <- nrow(x_i)
    if (n != length(y_i)) {
      stop("Cohort ", i, ": `x` and `y` have different sample sizes.", call. = FALSE)
    }

    # Stratified sampling within each class
    yes_idx <- which(y_i == "Yes")
    no_idx <- which(y_i == "No")

    # Partition each class
    partition_class <- function(idx) {
      n_class <- length(idx)
      if (n_class == 0L) {
        return(list(train = integer(0), val = integer(0), heldout = integer(0)))
      }

      shuffled <- sample(idx)
      n_train <- max(1L, round(n_class * train_ratio))
      n_val <- max(1L, round(n_class * val_ratio))
      n_heldout <- n_class - n_train - n_val

      # Ensure at least 1 sample in held-out if possible

      if (n_heldout < 1L && n_class >= 3L) {
        n_heldout <- 1L
        if (n_train > n_val) {
          n_train <- n_train - 1L
        } else {
          n_val <- n_val - 1L
        }
      }

      list(
        train = shuffled[seq_len(n_train)],
        val = if (n_val > 0L) shuffled[(n_train + 1L):min(n_train + n_val, n_class)] else integer(0),
        heldout = if (n_heldout > 0L) shuffled[(n_train + n_val + 1L):n_class] else integer(0)
      )
    }

    yes_parts <- partition_class(yes_idx)
    no_parts <- partition_class(no_idx)

    train_idx <- c(yes_parts$train, no_parts$train)
    val_idx <- c(yes_parts$val, no_parts$val)
    heldout_idx <- c(yes_parts$heldout, no_parts$heldout)

    # Check for small partitions
    min_per_class <- 2L
    warn_threshold <- 5L

    check_partition_size <- function(part_name, yes_n, no_n) {
      if (yes_n < min_per_class || no_n < min_per_class) {
        stop(
          "Cohort '", cohort_names[i], "' has insufficient samples in ", part_name, " set. ",
          "Yes: ", yes_n, ", No: ", no_n, ". Minimum required: ", min_per_class, " per class.",
          call. = FALSE
        )
      }
      if (yes_n < warn_threshold || no_n < warn_threshold) {
        warning(
          "Cohort '", cohort_names[i], "' has few samples in ", part_name, " set. ",
          "Yes: ", yes_n, ", No: ", no_n, ". Results may be unreliable.",
          call. = FALSE
        )
      }
    }

    check_partition_size("training", length(yes_parts$train), length(no_parts$train))
    check_partition_size("validation", length(yes_parts$val), length(no_parts$val))
    check_partition_size("held-out", length(yes_parts$heldout), length(no_parts$heldout))

    # Store partitions
    train_x[[i]] <- x_i[train_idx, , drop = FALSE]
    train_y[[i]] <- y_i[train_idx]
    val_x[[i]] <- x_i[val_idx, , drop = FALSE]
    val_y[[i]] <- y_i[val_idx]
    heldout_x[[i]] <- x_i[heldout_idx, , drop = FALSE]
    heldout_y[[i]] <- y_i[heldout_idx]

    # Record partition info
    partition_info$n_total[i] <- n
    partition_info$n_train[i] <- length(train_idx)
    partition_info$n_val[i] <- length(val_idx)
    partition_info$n_heldout[i] <- length(heldout_idx)
    partition_info$train_yes[i] <- length(yes_parts$train)
    partition_info$train_no[i] <- length(no_parts$train)
    partition_info$val_yes[i] <- length(yes_parts$val)
    partition_info$val_no[i] <- length(no_parts$val)
    partition_info$heldout_yes[i] <- length(yes_parts$heldout)
    partition_info$heldout_no[i] <- length(no_parts$heldout)
  }

  names(train_x) <- cohort_names
  names(train_y) <- cohort_names
  names(val_x) <- cohort_names
  names(val_y) <- cohort_names
  names(heldout_x) <- cohort_names
  names(heldout_y) <- cohort_names

  list(
    train_x = train_x,
    train_y = train_y,
    val_x = val_x,
    val_y = val_y,
    heldout_x = heldout_x,
    heldout_y = heldout_y,
    cohort_names = cohort_names,
    partition_info = partition_info
  )
}
