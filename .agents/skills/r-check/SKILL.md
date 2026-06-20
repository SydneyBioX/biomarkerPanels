---
name: r-check
description: Runs roxygen2::roxygenise() and devtools::check() via check.sh
---

# R Package Check Workflow

Use this skill when you need to document the package and run a clean R package check using `devtools::check()`. 
This workflow ensures that the environment PATH avoids interference from custom `~/.zshrc` compiler toolchains on macOS.

## Instructions

1.  **Format Code (Recommended)**
    Before checking, ensure code formatting complies with the repository standards by running:
    ```bash
    air format .
    ```
    *(Note: Ensure the `air` command is installed on your system).*

2.  **Generate Documentation**
    First, ensure the documentation is updated by executing:
    ```bash
    Rscript -e "roxygen2::roxygenise()"
    ```

3.  **Execute the Check Script**
    Next, execute the custom check script:
    ```bash
    chmod +x check.sh
    ./check.sh
    ```
    This script (`check.sh`) sets a clean `PATH` and calls `devtools::check(document = FALSE)`.

4.  **Check Pkgdown Documentation (Optional but Recommended)**
    To ensure the website documentation renders without errors, run:
    ```bash
    Rscript -e "pkgdown::check_pkgdown()"
    ```

5.  **Review the Output**
    Wait for the command to finish. If there are any `WARNING`s or `ERROR`s reported during the check, review the output logs or the `biomarkerPanels.Rcheck` directory to diagnose the problems.

## Notes
* Both `check.sh` and `check_output.txt` are included in `.gitignore` and `.Rbuildignore` so they do not accidentally get included in the R package build.
