## code to prepare `fake_gene_expression` dataset goes here

pkgload::load_all()

# Generate a small but realistic multi-cohort dataset
# We use p=250 to keep the file size well under Bioconductor limits while
# still being useful for vignettes and tests.
fake_gene_expression <- simulate_expression_data(
  p = 250L,
  n = 100L,
  k = 4L,
  seed = 20240220
)

usethis::use_data(fake_gene_expression, overwrite = TRUE)
