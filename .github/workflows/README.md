# GitHub Actions

Adapted from the TreeTools workflow suite to match Consensus's situation:
a pure-R package (`R >= 4.1.0`) with no compiled code *yet*.

## Active workflows

| Workflow | Purpose |
| --- | --- |
| `R-CMD-check.yml` | `R CMD check` on Windows, macOS (Intel + ARM), Linux (ARM release + devel), plus a `--as-cran` release check. Runs on push/PR and a weekly Thursday cron. Coverage uploaded to Codecov from the Windows job. |
| `pkgdown.yml` | Builds and deploys the documentation site; swaps `dev/` docs for redirects on release. |
| `codemeta.yml` | Regenerates `codemeta.json` when `DESCRIPTION` changes. |
| `r-universe-text.yml` | Tests the r-universe build. |
| `copilot-setup-steps.yml` | Pre-installs the dev toolchain for GitHub Copilot coding agent. |
| `rhub.yaml` | Manual (`workflow_dispatch`) R-hub multi-platform checks. |

## Required repository secrets

- `CODECOV_TOKEN` — for the coverage upload in `R-CMD-check.yml`.
- `RHUB_TOKEN` — for `rhub.yaml` (`rhub::rhub_setup()` sets this up).

## Deferred until C++ lands

The DESCRIPTION anticipates incorporating C++ from FACT/FACT2/FDCT, but there is
no `src/` directory yet. The following TreeTools workflows are intentionally
**not** included, because their path filters also trigger on
`tests/testthat/**.R` — with no compiled code they would misfire on every test
edit (and `RcppDeepState` sets `fail_ci_if_error: true`):

- `memcheck.yml` (valgrind)
- `ASan.yml` (AddressSanitizer)
- `rchk.yml`
- `RcppDeepState.yml`
- `benchmark.yml` (also needs a `benchmark/` harness)

Add these in the same change that introduces `src/`.

Also skipped: `update-csl.yml` (downloads a TreeTools-specific CSL file Consensus
does not use).
