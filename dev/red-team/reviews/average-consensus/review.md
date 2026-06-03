# Red-team review: Average consensus
**Date:** 2026-06-03
**Reviewed:** `R/Average.R`, `tests/testthat/test-Average.R` (no C++ backing)
**Verdict:** CONDITIONAL SHIP — two MAJOR input-validation holes and two MINOR silent-misbehaviours confirmed empirically; no correctness bug in the algorithm path itself.

All findings were reproduced by `probe-edge.R` and `_probe-*.R` in this directory.

---

## F1 — n<3 tips crashes with an opaque ape error (MAJOR)

`Average()` has no guard on leaf count. When `length(trees) > 1` and the trees
have fewer than 3 tips, execution reaches `ape::nj` / `ape::bionj` /
`ape::fastme.bal` / `ape::fastme.ols`, which each throw their own unhelpful
errors ("cannot build an NJ tree with less than 3 observations").

Separately, the **single-tree early-return** path calls `.RootResult()`, which
calls `ape::unroot()`. A 2-tip or 1-tip tree passed as `Average(list(tr))` or
`Average(tr)` crashes with "cannot unroot a tree with less than three edges"
— again from ape, with no `Average()`-level context.

Reproduced:
```
Average(list(read.tree("(A,B);"), read.tree("(A,B);")), method = "nj")
# Error: cannot build an NJ tree with less than 3 observations

Average(list(read.tree("(A,B);")))
# Error: cannot unroot a tree with less than three edges.
```

The other split-based methods all route through `.PrepareTrees()` which has an
explicit `< 4` leaf guard returning a trivial shortcut; `Average()` has no
equivalent guard.

**Fix:** add a leaf-count check in `.AverageTrees()` (or `Average()`) that
stops with a clear message before reaching the distance fitters.

---

## F2 — `check.labels = FALSE` with a superset silently drops taxa (MAJOR)

When `check.labels = FALSE` is passed, `RenumberTips()` is skipped but the
label-set check is also skipped. The average is computed over the labels of
`trees[[1]]` only (via `cophenetic.phylo(tr)[labs, labs]`). If a subsequent
tree has *more* tips than `trees[[1]]`, the extra tips are silently ignored:
`ape::cophenetic.phylo` returns a named matrix, and `[labs, labs]` indexes
only the names in `labs` — no warning, no error.

Reproduced:
```
t4 <- ape::read.tree(text = "((A,B),(C,D));")       # 4 tips
t5 <- ape::read.tree(text = "(((A,B),(C,D)),E);")   # 5 tips — E is dropped
Average(list(t4, t5), check.labels = FALSE)
# returns a 4-tip tree, E silently absent
```

With a *different* (non-superset) mismatch, e.g. `{A,B,C,D}` vs `{A,B,C,E}`,
`cophenetic.phylo(t2)["D",]` throws "subscript out of bounds" — so the damage
is inconsistent: superset → silent wrong result; disjoint → cryptic error.

**Fix:** Even with `check.labels = FALSE`, assert that every tree's tip labels
are a superset of `labs1` and report any discrepancy. Or at minimum document
the restriction clearly.

---

## F3 — Non-logical `edgeLengths` silently falls through to edge-counting (MINOR)

`edgeLengths` is documented as a logical. The implementation uses:

```r
useLengths <- if (is.na(edgeLengths)) all(haveLengths) else isTRUE(edgeLengths)
```

`isTRUE()` requires a length-1 `TRUE`; anything else returns `FALSE`. So
`edgeLengths = "yes"` and `edgeLengths = 1L` both silently fall through to
`useLengths = FALSE` (count-edges mode) with no warning or error.

Reproduced:
```
# edgeLengths = "yes" silently counts edges (same as FALSE)
Average(rmtree(3, 6), edgeLengths = "yes")  # no error
# edgeLengths = 1L  silently counts edges
Average(rmtree(3, 6), edgeLengths = 1L)     # no error
```

`edgeLengths = c(TRUE, FALSE)` does error ("the condition has length > 1"),
but only because the `if(is.na(...))` test triggers R's base vector-condition
warning that becomes an error in recent R — so protection from length-2 input
is incidental, not intentional.

**Fix:** validate `edgeLengths` explicitly: `stopifnot(is.logical(edgeLengths) && length(edgeLengths) == 1L)`.

---

## F4 — NA branch lengths give a confusing ape error (INFO)

A tree with one `NA` edge length passes the `edgeLengths = NA` (default) path
and `cophenetic.phylo` propagates NAs into the distance matrix, which the
fitters then reject with "missing values are not allowed in the distance
matrix / Consider using njs()". The error is correct in that it does error, but
the message mentions `njs()` and gives no hint that the problem is a partial NA
in a branch-length vector.

No fix required unless a cleaner message is wanted; `edgeLengths = FALSE` is
the documented workaround for trees without reliable branch lengths.

---

## Verified fine

- **Algorithm (core path).** The two-step Lapointe & Cucumel (1997) procedure is
  faithfully implemented: `cophenetic.phylo` computes patristic distances
  (path-length distances with branch lengths, edge counts without); the weighted
  average matrix is computed as a convex combination (`weights / sum(weights)`);
  `fastme.bal` / `nj` / `bionj` / `fastme.ols` / `LeastSquaresTree` are each a
  reasonable choice for step 2.

- **Scale = "max" guard.** `if (maxVal > 0)` correctly prevents a division-by-zero
  when all branch lengths are zero (the all-zero matrix is left as-is, which is
  correct — NJ on a zero matrix gives an unresolved star). Verified by P3.

- **Permutation invariance.** Reordering tips in every input tree gives identical
  output distances (to 1e-10). The `cophenetic.phylo(tr)[labs, labs]` indexing is
  safe because `cophenetic` returns a named matrix. Verified by P5 and P16.

- **Weights edge cases.** `weights = c(1, 0)` produces a distance matrix equal to
  tree 1's (RSS = 0); `weights = c(0, 0)` errors. Verified by P9 and P13.

- **Identical-tree idempotence.** `Average(list(t, t, t), method = "ls")` recovers
  the generating tree with RSS ≈ 7e-30 (numerical zero). Verified by P10.

- **Single-tree path (n ≥ 3 tips).** Returns unrooted, preserves distances exactly.
  Verified by P11.

- **`method = "ls"` guard.** The `requireNamespace` / `exists("LeastSquaresTree")`
  double-check is correct and unreachable in normal installs (nocov). Validated.

- **`lsControl` validation and forwarding.** A non-list errors; a non-empty list is
  correctly forwarded via `modifyList`. Verified by test suite.

---

## Coverage gaps

1. **n < 3 tips** — no test at all (now F1 above); the test suite starts at 5 or 6
   leaves.

2. **`check.labels = FALSE` with mismatched / superset tip sets** — no test. F2
   above shows this can silently drop taxa.

3. **`edgeLengths` with non-logical input** ("yes", 1L) — no test; F3 above shows
   these silently count edges.

4. **Polytomous (non-binary) input trees.** The distance methods (NJ etc.) accept
   any `cophenetic` output regardless of tree resolution. Star trees with real branch
   lengths work (P3 uses zero-length star). A test with a genuine polytomy would close
   this gap.

5. **`outgroup` that is not a tip in the result.** `ape::root()` raises an error;
   no test checks this.

6. **Trees with a single-tip polytomy or multifurcating root.** Not tested.

7. **Weighted average with `scale = "max"`** — no test combining both.
