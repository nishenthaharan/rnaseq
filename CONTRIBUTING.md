# Contributing to RNASeq Explorer

Thank you for improving RNASeq Explorer. Contributions should strengthen the
scientific correctness, reproducibility, usability, or maintainability of the
application.

## Before starting

1. Search the issue tracker for existing work.
2. Open or comment on an issue describing the problem and acceptance criteria.
3. Keep each pull request focused on one reviewable change.

Never commit patient information, unpublished research data, access tokens, or
other confidential material. Use simulated or fully anonymised test fixtures.

## Development setup

RNASeq Explorer requires R 4.3 or newer. From the repository root, install the
dependencies and run the application:

```bash
Rscript install_dependencies.R
R -q -e 'shiny::runApp()'
```

The Docker workflow is an alternative when a local R environment is not
available:

```bash
docker build -t rnaseq-explorer-r .
docker run --rm -p 3838:3838 rnaseq-explorer-r
```

## Testing

Add regression tests for bug fixes and focused tests for new behaviour. Run:

```bash
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
Rscript -e 'invisible(lapply(c("app.R", list.files("R", pattern = "[.]R$", full.names = TRUE)), parse))'
```

For changes affecting the container, confirm that `docker build .` succeeds.

## Commit and pull-request quality

- Use an imperative, scoped message such as `fix(io): reject overflowing counts`.
- Separate tests, implementation, and documentation when they are independently
  reviewable; do not split trivial changes merely to create more commits.
- Explain scientific assumptions and user-visible behaviour in the pull request.
- Link the issue with `Closes #123` when the work fully resolves it.
- Keep generated files, credentials, and local R state out of the repository.

## Scientific changes

Changes to filtering, statistical models, transformations, classifications, or
plots must describe their rationale and expected interpretation. Where
possible, cite the relevant package documentation or peer-reviewed method and
include a deterministic test fixture.
