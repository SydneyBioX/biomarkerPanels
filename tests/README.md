# Test Data Assets

The `data/` directory stores reproducible fixtures used by future tests.
`fake_gene_expression.Rds` mirrors the structure produced by `simulate_gene_expression.R`:
- `x_list`: four `100 x 5000` log-expression matrices with dataset-specific shifts.
- `y_list`: matching binary response factors (`"No"`/`"Yes"`).
- `metadata`: generation parameters for reproducibility.

Regenerate the fixture by running `Rscript simulate_gene_expression.R` and copying the
resulting `simulated_gene_expression.Rds` into this folder with the same filename.
