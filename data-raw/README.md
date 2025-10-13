# Data Preparation Scripts

This directory holds scripts that (re)generate datasets included in the
package. Scripts should read from raw sources, perform necessary processing,
and write summarized objects into `inst/extdata/` or serialize `.rda` objects
into `data/`.

Use `usethis::use_data_raw()` to scaffold additional data pipelines as the
package grows.
