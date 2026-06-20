---
name: format-code
description: Enforce coding standards automatically using air format . and rules checks
---

# Format Code Workflow

Use this skill to automatically enforce coding standards for the project. 

## Instructions

1.  **Run Formatter**
    Run the project's default formatter on all files in the current directory:
    ```bash
    air format .
    ```
    *Note: Ensure `air` is installed or accessible in your current environment. If `zsh: command not found: air` occurs, wait for the user to install it or skip this step if instructed.*

2.  **Verify Pipe Operators**
    Scan any newly created or modified files to ensure the base pipe operator (`|>`) is used instead of the magrittr pipe (`%>%`). This is a strict project requirement.

3.  **Verify Anonymous Functions**
    Scan for single-line anonymous functions (e.g. `\(x) x+1` or `function(x) x+1`) and ensure they strictly use the `function() {...}` syntax.

4.  **Report Violations**
    Fix any non-compliant code using IDE editing tools, then report the changes made to the user.
