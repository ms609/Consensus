# Red-team review: package infrastructure & wrappers
**Date:** 2026-06-03
**Verdict:** CONDITIONAL SHIP — one CRAN-blocking NOTE (stale NAMESPACE import),
one stale `@importFrom` that will produce a second NOTE, no correctness bugs.
Wrappers, RcppExports, DESCRIPTION, and Rd files are all sound.

---

## F1 — Stale `importFrom(TreeTools,CompatibleSplits)` in NAMESPACE (MAJOR)

**Location:** `NAMESPACE:23`

`CompatibleSplits` is listed as `importFrom(TreeTools,CompatibleSplits)` in
`NAMESPACE` but there is no `@importFrom TreeTools CompatibleSplits` tag
anywhere in `R/`.  The function is used only in `tests/testthat/test-loose.R`
(via `TreeTools::` qualified calls) and in `dev/` scripts — never in package
source code.

`R CMD check --as-cran` will fire a NOTE:

```
NOTE: package imports 'CompatibleSplits' from 'TreeTools' but does not use it
```

This is a CRAN-blocking NOTE for a first submission.  NAMESPACE is out of sync
with the `@importFrom` tags: running `roxygen2::roxygenise()` would remove this
line.

**Fix:** run `roxygen2::roxygenise()`.  Alternatively, delete the line manually:
```
importFrom(TreeTools,CompatibleSplits)
```
(Do not add a redundant `@importFrom` in R/ — the function is not used in
package code.)

---

## F2 — Stale `@importFrom ape write.tree` in `R/local.R` (MINOR)

**Location:** `R/local.R:53`

```r
#' @importFrom ape read.tree write.tree
```

`write.tree` is declared but never called in `R/local.R` (confirmed by
searching all of `R/`).  It is used only in `dev/oracle/local/oracle_local.R`
and other `dev/` scripts, which are not shipped.

Running `roxygen2::roxygenise()` would silently retain this line because
roxygen trusts the explicit `@importFrom` tag and does not verify usage.  It
will appear in `NAMESPACE` as `importFrom(ape,write.tree)` and generate a
second NOTE on CRAN check:

```
NOTE: package imports 'write.tree' from 'ape' but does not use it
```

**Fix:** Remove `write.tree` from that `@importFrom` tag:
```r
#' @importFrom ape read.tree
```

---

## F3 — `MR` alias not tested (INFO)

**Location:** `tests/testthat/test-wrappers.R`

`test-wrappers.R` checks `MajorityRule(trees) == Majority(trees)` but never
exercises `MR`.  Since `MR <- Majority` these are the same object and the risk
is low.  However, a future accidental redefinition of `MR` (e.g. via a
namespace collision with another package that exports `MR`) would not be caught
by the suite.

**Fix:** add one line to the existing test:
```r
expect_equal(MR(trees), Majority(trees))
```

---

## F4 — No test for bare `phylo` passed to `Strict()` or `Majority()` (INFO)

**Location:** `tests/testthat/test-wrappers.R`

The single-tree test wraps the tree in `list(tree)`.  `TreeTools::Consensus`
accepts a bare `phylo` and returns it unchanged (confirmed by tracing the
function body).  The wrappers delegate directly, so a bare `phylo` is handled
correctly — but this is not exercised by any test.

**Suggested addition:**
```r
test_that("Strict() and Majority() accept a bare phylo", {
  tree <- ape::as.phylo(0, 8)
  expect_equal(Strict(tree), tree)
  expect_equal(Majority(tree), tree)
})
```

---

## Verified fine

### Wrapper correctness
- `Strict(trees)` -> `Consensus(trees, p = 1)`: correct.
- `Majority(trees, p = 0.5)` -> `Consensus(trees, p = p)`: correct.  The `p =
  0.5` default means "strictly more than 50%" (verified empirically: two
  conflicting trees with each split in exactly 50% return a star).  The
  docstring "more than half" is accurate.
- `MajorityRule <- Majority` and `MR <- Majority`: correct aliases, both
  exported in NAMESPACE, both aliased in `Majority.Rd`.
- Delegation to `TreeTools::Consensus` means tip-label validation, `p` range
  validation (`p < 0.5 || p > 1` error), and single-tree handling are all
  inherited from TreeTools.

### NAMESPACE vs. `@importFrom` (excluding F1 and F2)
Every other function in NAMESPACE traces back to a `@importFrom` tag in `R/`:
- `TreeTools::Consensus` (wrappers.R, Quartet.R)
- `TreeTools::DescendantEdges`, `NTip`, `Preorder`, `RootTree` (selection.R, adams.R)
- `TreeTools::KeepTip` (adams.R)
- `TreeTools::NSplits` (Quartet.R `@examples`)
- `TreeTools::RenumberTips`, `TipLabels`, `StarTree` (multiple files)
- `TreeTools::as.Splits` (BHV.R, Quartet.R)
- `ape::as.phylo` (BHV.R, Quartet.R — generic needed for S3 dispatch to TreeTools method)
- `ape::read.tree` (selection.R, adams.R, rstar.R, local.R)
- `ape::collapse.singles`, `is.rooted`, `unroot` (BHV.R)
- `ape::bionj`, `cophenetic.phylo`, `fastme.bal`, `fastme.ols`, `nj`, `root`, `unroot` (Average.R)
- `Rcpp::sourceCpp` (ConsTree-package.R, BHV.R)
- `Rdpack::reprompt` (ConsTree-package.R)
- `stats::as.dist` (BHV.R, Average.R)
- `utils::modifyList` (Average.R)

### RcppExports consistency
`R/RcppExports.R` and `src/RcppExports.cpp` are fully consistent: 12 functions
each, matching names, argument names, and types.  The generator token
(`10BE3573-1514-4C36-9D1C-5A225CD40393`) is identical in both files.
`Rcpp::compileAttributes()` would not change either file.

### DESCRIPTION correctness
- All `Imports:` packages (`ape >= 5.6`, `Rcpp >= 1.0.0`, `Rdpack >= 2.6`,
  `TreeTools >= 2.1.0`) are used in package code.
- `LinkingTo: Rcpp` is sufficient: all C++ includes are standard library +
  `<Rcpp.h>` only (no BH, RcppArmadillo, or other header packages needed).
- `TreeSearch` is correctly in `Suggests` (used only when `method = "ls"` in
  `Average()`, guarded by `requireNamespace()`).
- `Quartet` is correctly in `Suggests` (used via `Quartet::` in tests, guarded
  by `skip_if_not_installed("Quartet")`).
- `knitr`, `rmarkdown`, `spelling`, `testthat` are all correctly in `Suggests`.
- The `Additional_repositories` field pointing to `https://ms609.github.io/packages/`
  is present for TreeTools pre-CRAN; this is correct practice.

### Documentation
- All exported functions have Rd pages generated by roxygen.
- `@examples` for `Strict()`, `Majority()`, `Loose()`, `Greedy()`,
  `MajorityPlus()`, `Frequency()` all use `ape::as.phylo(0:5, 8)` — runnable
  without qualification issues.
- `cpp_quartet_consensus` is `\keyword{internal}` and not exported — correct.
- `Average()` `@examples` wraps the `method = "ls"` path in `\donttest{}` with
  a `requireNamespace("TreeSearch")` guard — correct for a Suggested dependency.
- `Majority.Rd` aliases `MajorityRule` and `MR` and shows their full `\usage`
  — correct.
- No `\dontrun{}` examples were found; all examples should run during `R CMD check`.

### Edge-case tracing for `Strict()` and `Majority()`
- **Single bare `phylo`**: `TreeTools::Consensus(phylo, p = 1)` returns the
  tree unchanged (verified).
- **Single-element list**: handled by the `length(trees) == 1L` guard in
  `TreeTools::Consensus` — returns `trees[[1]]`.
- **NULL entries**: `TreeTools::Consensus` coerces via `lapply(c(trees), ...)`,
  where `c(list(tree, NULL))` drops the NULL (R list concatenation drops
  NULLs), so `Strict(list(tree, NULL))` silently drops the NULL and returns
  `tree`. This is consistent with `.PrepareTrees()` behaviour in the FACT-family
  methods.
- **Mismatched tip labels**: `TreeTools::Consensus` calls `RenumberTips` and
  silently drops tips not in the smallest tree (warning emitted). The R
  wrappers inherit this behaviour without extra validation — acceptable for
  methods that delegate entirely to TreeTools.
- **`p` out of range**: `TreeTools::Consensus` validates `p >= 0.5 && p <= 1`
  and throws an error — so `Majority(trees, p = 0.3)` is properly rejected.

---

## Coverage gaps

1. **`MR` alias** — no test calls `MR()` directly (see F3).
2. **Bare `phylo` input** to `Strict()` / `Majority()` — delegated correctly
   but untested (see F4).
3. **`p` boundary** — no test exercises `Majority(trees, p = 1)` returning the
   same result as `Strict(trees)` for a non-trivial conflict case (the existing
   test uses this but only for equality, not to confirm splits are dropped).
4. **Tip-label mismatch** in `Strict()` / `Majority()` — the TreeTools warning
   path ("removing leaves not in smallest") is never exercised by the test suite.

---

## Notes for next reviewer of this area

- The stale NAMESPACE entry (F1) is the only issue that would block CRAN
  submission; both F1 and F2 are fixed by running `roxygen2::roxygenise()` and
  removing `write.tree` from `local.R:53`.
- The wrappers are extremely thin (2 lines each) so there is very little that
  can go wrong in them; the real complexity is in `selection.R` and the C++
  ports reviewed in earlier slots.
- `Rdpack` glue in `ConsTree-package.R` is the standard pattern; the
  `inst/REFERENCES.bib` was not reviewed in this pass — a future round could
  check BibTeX key consistency against all `\insertCite{}` usages.
- The `ape::as.phylo` import in `BHV.R` and `Quartet.R` registers the generic
  so S3 dispatch reaches `TreeTools::as.phylo.Splits` — this is correct and
  not a duplication issue.
