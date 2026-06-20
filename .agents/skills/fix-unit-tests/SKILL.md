---
name: fix-unit-tests
description: Fix unit tests
---

# Fix Unit Tests Workflow

When following this workflow to finish fixing the unit tests, perform the following steps for each unit test file (`.R` file):

1. **Document Expected Behavior**: For each of the `.R` file and each of the unit test, make sure that each of the tests is well-documented in terms of what the expected behavior is. Focus on clarity and ensuring the purpose of the test is easily understandable. Strive to keep tests minimal with few comments, only adding what's necessary for clarity.

2. **Execute and Verify**: Execute the specific test file using `devtools::test_active_file()` instead of the full test suite for faster iteration.
   ```bash
   Rscript -e "devtools::test_active_file('tests/testthat/test-{name}.R')"
   ```
   Check that the test can be executed without an error.

3. **Handle Failures**: If the test returns an error or failure:
   - Document a potential solution to fixing that particular unit test.
   - Alternatively, propose removing that test altogether if it is no longer relevant.
   - Put that test into an implementation plan if it requires substantial changes.

4. **Format Code**: After making any changes to the test code, you must run `air format .` to ensure the changes adhere to the project's formatting rules. Ensure the `air` command is installed on the system before executing it.
