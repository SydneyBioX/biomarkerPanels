#' Simulate Log-scale Gene Expression Datasets with Distributional Shifts
#'
#' Generate multiple expression matrices with dataset-specific mean shifts and
#' paired binary responses. Useful for benchmarking panel optimization routines
#' under `n << p` conditions.
#'
#' @param p Number of features (genes).
#' @param n Number of samples per dataset.
#' @param k Number of datasets to simulate.
#' @param seed Optional random seed for reproducibility.
#' @param shift_scale Magnitude of dataset-specific feature shifts.
#' @return A list with `x_list`, `y_list`, and `metadata`.
#' @export
simulate_expression_data <- function(p = 5000L, n = 100L, k = 4L, seed = NULL,
                                     shift_scale = 0.6) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  dataset_ids <- seq_len(k)
  gene_names <- sprintf("gene_%04d", seq_len(p))

  base_mean <- stats::rnorm(p, mean = 5, sd = 0.75)
  feature_sd <- stats::runif(p, min = 0.4, max = 0.9)

  make_shift_vector <- function() {
    shift <- numeric(p)
    idx <- sample(p, size = floor(0.12 * p))
    direction <- sample(c(-1, 1), length(idx), replace = TRUE)
    shift[idx] <- stats::rnorm(length(idx), mean = direction * shift_scale, sd = 0.1)
    shift
  }

  global_shift <- seq(-0.4, 0.4, length.out = k)
  shift_matrix <- vapply(dataset_ids, function(i) make_shift_vector(), numeric(p))

  generate_dataset <- function(mean_vec) {
    noise <- matrix(stats::rnorm(n * p), nrow = n)
    scaled <- sweep(noise, 2, feature_sd, `*`)
    sweep(scaled, 2, mean_vec, `+`)
  }

  informative_idx <- sample(p, 60L)
  beta <- stats::rnorm(length(informative_idx), mean = 0, sd = 0.25)

  generate_labels <- function(x_matrix) {
    eta <- as.numeric(x_matrix[, informative_idx, drop = FALSE] %*% beta)
    eta <- eta - mean(eta)
    eta <- eta + stats::rnorm(length(eta), sd = 0.3)
    prob <- stats::plogis(eta)
    labels <- ifelse(stats::runif(length(prob)) < prob, "Yes", "No")
    factor(labels, levels = c("No", "Yes"))
  }

  x_list <- vector("list", k)
  y_list <- vector("list", k)

  for (i in dataset_ids) {
    mean_vec <- base_mean + shift_matrix[, i] + global_shift[i]
    x_mat <- generate_dataset(mean_vec)
    colnames(x_mat) <- gene_names
    rownames(x_mat) <- sprintf("sample_%03d", seq_len(n))
    x_list[[i]] <- x_mat
    y_list[[i]] <- generate_labels(x_mat)
  }

  names(x_list) <- sprintf("x%d", dataset_ids)
  names(y_list) <- sprintf("y%d", dataset_ids)

  list(
    x_list = x_list,
    y_list = y_list,
    metadata = list(
      seed = seed,
      n_samples = n,
      n_features = p,
      informative_genes = gene_names[informative_idx],
      global_shift = global_shift
    )
  )
}
