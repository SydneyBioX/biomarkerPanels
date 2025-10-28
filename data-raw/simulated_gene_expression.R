## code to prepare `simulated_gene_expression` dataset goes here

#!/usr/bin/env Rscript

# Simulate four log-scale gene expression datasets with dataset-specific shifts
# and save them alongside binary response labels.

set.seed(20240220)

n_features <- 500L
n_samples <- 100L
n_datasets <- 4L
dataset_ids <- seq_len(n_datasets)
gene_names <- sprintf("gene_%04d", seq_len(n_features))

base_mean <- rnorm(n_features, mean = 5, sd = 0.75)
feature_sd <- runif(n_features, min = 0.4, max = 0.9)

make_shift_vector <- function(num_features, prop_shift = 0.12, shift_scale = 0.6) {
  shift <- numeric(num_features)
  idx <- sample(num_features, size = floor(prop_shift * num_features))
  direction <- sample(c(-1, 1), length(idx), replace = TRUE)
  shift[idx] <- rnorm(length(idx), mean = direction * shift_scale, sd = 0.1)
  shift
}

global_shift <- seq(-0.4, 0.4, length.out = n_datasets)
shift_matrix <- vapply(dataset_ids, function(i) make_shift_vector(n_features), numeric(n_features))

generate_dataset <- function(mean_vector) {
  noise <- matrix(rnorm(n_samples * n_features), nrow = n_samples)
  scaled <- sweep(noise, 2, feature_sd, `*`)
  sweep(scaled, 2, mean_vector, `+`)
}

signal_gene_count <- 12L
informative_idx <- seq_len(signal_gene_count)
informative_genes <- gene_names[informative_idx]
beta <- rep(c(1.5, -1.5), length.out = length(informative_idx))
signal_shift <- rep(c(2.0, -2.0), length.out = length(informative_idx))

generate_labels <- function(x_matrix) {
  signal_block <- scale(x_matrix[, informative_idx, drop = FALSE], center = TRUE, scale = TRUE)
  eta <- as.numeric(signal_block %*% beta)
  eta <- eta - mean(eta)
  eta_sd <- sd(eta)
  if (eta_sd > 1e-08) {
    eta <- eta / eta_sd
  }
  eta <- eta + rnorm(length(eta), sd = 0.25)
  prob <- plogis(eta)
  labels <- ifelse(runif(length(prob)) < prob, "Yes", "No")
  factor(labels, levels = c("No", "Yes"))
}

x_list <- vector("list", n_datasets)
y_list <- vector("list", n_datasets)

for (i in dataset_ids) {
  mean_vec <- base_mean + shift_matrix[, i] + global_shift[i]
  x_mat <- generate_dataset(mean_vec)
  colnames(x_mat) <- gene_names
  rownames(x_mat) <- sprintf("sample_%03d", seq_len(n_samples))
  y_vec <- generate_labels(x_mat)
  yes_idx <- y_vec == "Yes"
  no_idx <- y_vec == "No"
  if (any(yes_idx)) {
    x_mat[yes_idx, informative_idx] <- sweep(
      x_mat[yes_idx, informative_idx, drop = FALSE],
      2, signal_shift, `+`
    )
  }
  if (any(no_idx)) {
    x_mat[no_idx, informative_idx] <- sweep(
      x_mat[no_idx, informative_idx, drop = FALSE],
      2, -signal_shift, `+`
    )
  }
  x_list[[i]] <- x_mat
  y_list[[i]] <- y_vec
}

names(x_list) <- sprintf("x%d", dataset_ids)
names(y_list) <- sprintf("y%d", dataset_ids)

simulated_gene_expression <- list(
  metadata = list(
    seed = 20240220,
    n_samples = n_samples,
    n_features = n_features,
    informative_genes = informative_genes,
    signal_coefficients = stats::setNames(beta, informative_genes),
    signal_shift = stats::setNames(signal_shift, informative_genes),
    global_shift = global_shift
  ),
  x_list = x_list,
  y_list = y_list
)

cat("Simulation complete. Saved to simulated_gene_expression.Rds\n")
usethis::use_data(simulated_gene_expression, compress = "xz")
