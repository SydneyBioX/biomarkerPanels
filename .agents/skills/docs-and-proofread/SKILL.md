---
name: docs-and-proofread
description: Automate R code documentation and text proofreading
---

# Documentation & Proofreading Workflow

Use this skill when documenting R code or proofreading extensive text files (such as vignettes, READMEs, or other Markdown files) to ensure they meet the project's quality standards.

## 1. Code Documentation Workflow

When creating or modifying R functions, perform the following:
*   **Export**: Ensure every user-facing function is exported. Internal functions should *not* be exported or have roxygen documentation.
*   **Roxygen Formatting**: Add or update `roxygen2` comments for the function. Ensure these comments are wrapped strictly at 80 characters.
*   **Pkgdown Metadata**: Whenever you add a new user-facing (non-internal) documentation topic, you must add the topic to `_pkgdown.yml`.
*   **Compile**: Always re-document the package after changing a roxygen2 comment:
    ```bash
    Rscript -e "devtools::document()"
    ```
    *(Alternatively, trigger the `r-check` skill if a full check is preferred).*

## 2. Text Proofreading Workflow

When asked to proofread a file, act as an expert proofreader and editor with a deep understanding of clear, engaging, and well-structured writing.
*   **Sentence Case & Dialect**: Use sentence case for headings and British English for spelling.
*   **Iterative Review**: Work paragraph by paragraph. Always start by making a TODO list that includes individual items for each top-level heading.
*   **Autofix Minor Issues**: Fix spelling, grammar, and other minor problems directly without asking the user.
*   **Flag Ambiguity**: Label any unclear, confusing, or ambiguous sentences with a `FIXME` comment.
*   **Report**: Only report what you have changed.
