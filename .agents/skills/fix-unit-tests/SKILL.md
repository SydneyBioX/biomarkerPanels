---
name: fix-unit-tests
description: Fix unit tests
---

When following this workflow to finish fixing the unit tests, perform the following steps for each unit test file (`.R` file):

1. **Document Expected Behavior**: For each of the `.R` file and each of the unit test, make sure that each of the tests is well-documented in terms of what the expected behavior is. Focus on clarity and ensuring the purpose of the test is easily understandable.

2. **Execute and Verify**: Check that the test actually can be executed in the unit test suite without an error.

3. **Handle Failures**: If the test returns an error or failure:
   - Put that test into an implementation plan.
   - Document a potential solution to fixing that particular unit test.
   - Alternatively, propose removing that test altogether if it is no longer relevant.
