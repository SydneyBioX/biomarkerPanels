---
name: r-check
description: Runs roxygen2::roxygenise() and devtools::check() via check.sh
---

# R Package Check Workflow

Use this skill when you need to document the package and run a clean R package check using `devtools::check()`. 
This workflow ensures that the environment PATH avoids interference from custom `~/.zshrc` compiler toolchains on macOS.

## Instructions

1.  **Generate Documentation**
    First, ensure the documentation is updated by executing:
    ```bash
    Rscript -e "roxygen2::roxygenise()"
    ```

2.  **Execute the Check Script**
    Next, execute the custom check script:
    ```bash
    chmod +x check.sh
    ./check.sh
    ```
    This script (`check.sh`) sets a clean `PATH` and calls `devtools::check(document = FALSE)`.

3.  **Review the Output**
    Wait for the command to finish. If there are any `WARNING`s or `ERROR`s reported during the check, review the output logs or the `biomarkerPanels.Rcheck` directory to diagnose the problems.

## Notes
* Both `check.sh` and `check_output.txt` are included in `.gitignore` and `.Rbuildignore` so they do not accidentally get included in the R package build.
