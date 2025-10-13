#!/usr/bin/env Rscript

# Simulate four log-scale gene expression datasets with dataset-specific shifts
# and save them alongside binary response labels.

set.seed(20240220)

n_features <- 5000L
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

informative_idx <- sample(n_features, 60L)
beta <- rnorm(length(informative_idx), mean = 0, sd = 0.25)

generate_labels <- function(x_matrix) {
  eta <- as.numeric(x_matrix[, informative_idx, drop = FALSE] %*% beta)
  eta <- eta - mean(eta)
  eta <- eta + rnorm(length(eta), sd = 0.3)
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
  x_list[[i]] <- x_mat
  y_list[[i]] <- generate_labels(x_mat)
}

names(x_list) <- sprintf("x%d", dataset_ids)
names(y_list) <- sprintf("y%d", dataset_ids)

output <- list(
  metadata = list(
    seed = 20240220,
    n_samples = n_samples,
    n_features = n_features,
    informative_genes = gene_names[informative_idx],
    global_shift = global_shift
  ),
  x_list = x_list,
  y_list = y_list
)

saveRDS(output, file = "simulated_gene_expression.Rds")

cat("Simulation complete. Saved to simulated_gene_expression.Rds\n")
