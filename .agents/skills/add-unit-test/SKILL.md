---
name: add-unit-test
description: Streamline test-driven development for R package unit tests
---

# Add Unit Test Workflow

Use this skill when you need to write new unit tests for the R package to ensure high test coverage and correctness.

## Instructions

1.  **Locate Test File**
    Determine the base name of the `R/{name}.R` file being tested. Its corresponding test file must be located at `tests/testthat/test-{name}.R`. If the file does not exist, create it.

2.  **Placement & Style**
    - Place the new test block (`test_that(...)`) next to similar existing tests in the file.
    - Write the test assertion to clearly document the expected behavior.
    - Strive to keep tests minimal with few comments. The test description (first argument to `test_that`) should serve as the primary documentation.

3.  **Execute**
    Run the specific test file active file to verify that the test passes:
    ```bash
    Rscript -e "devtools::test_active_file('tests/testthat/test-{name}.R')"
    ```

4.  **Format**
    After making changes or successfully adding the test, ensure the file is properly formatted:
    ```bash
    air format .
    ```
    *(Note: Ensure `air` is installed on your system).*

5.  **Troubleshooting**
    If the test fails, use the `fix-unit-tests` skill to iterate on a solution.
