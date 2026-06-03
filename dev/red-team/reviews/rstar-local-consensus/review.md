# Red-Team Review: R* Consensus and Local Consensus

**Slot:** 6 of 6 (final)
**Date:** 2026-06-03
**Reviewer:** automated adversarial pass (slot 6 of 6)
**Verdict:** CONDITIONAL SHIP — one major silent-wrong-result bug (F1), one
medium silent-wrong-result path (F2), one minor dead code (F3), two coverage
gaps (F4, F5).

## Scope

Files reviewed:
- `src/rstar.cpp` — R* triplet-tally + strong-cluster assembly (~330 lines)
- `R/rstar.R` — `RStar()` R wrapper
- `tests/testthat/test-rstar.R` — shipped test suite
- `src/local_consensus.cpp` — MinRLC/MinILC exact-exponential port (~386 lines)
- `R/local.R` — `Local()` R wrapper
- `tests/testthat/test-local.R` — shipped test suite

Related context consulted: `R/selection.R` (`.PrepareTrees()` pattern),
previous slot reviews (slots 1–5).

---

## F1 — Missing cross-tree taxon-label validation — MAJOR

**Where:** `R/rstar.R:107`, `R/local.R:83` (identical pattern in both wrappers).

**What:**
Both `RStar()` and `Local()` derive the canonical label vector from the **first**
tree only:
```r
labels <- TipLabels(trees[[1L]])
edgeList <- lapply(trees, function(tr) {
  tr <- Preorder(RenumberTips(tr, labels))  # ← no check that tr has the same labels
  tr[["edge"]]
})
```
`TreeTools::RenumberTips(tr, labels)` reorders `tr`'s tips to match `labels`.
If `tr` carries a taxon absent from `labels` (or vice versa), the mismatch is
handled silently: TreeTools maps each of `tr`'s labels to the first position in
`labels` that matches, leaving unmatched taxa at an arbitrary or NA slot.  The
resulting integer edge matrix has tip indices that do not correspond to the
intended taxa.

Consequence for `RStar()`: the O(k n³) triplet tally accumulates counts for
the wrong leaf identities.  `R_maj` is corrupted, the strong-cluster
computation sees spurious or missing relationships, and the returned tree is
wrong — no warning or error is emitted.

Consequence for `Local()`: the common-triplet table is built from wrong tip
orderings, the Aho-graph edges connect wrong pairs, and the DP produces the
wrong (minimum local) consensus, silently.

**Repro sketch:**
```r
library(ConsTree)
t1 <- ape::read.tree(text = "(((a,b),c),d);")
t2 <- ape::read.tree(text = "(((a,b),c),X);")   # X not in t1
# Currently: no error, produces a tree with wrong topology
RStar(list(t1, t2))
Local(list(t1, t2))
```

**Fix:** add a per-tree taxon check after removing NULLs:
```r
if (any(vapply(trees[-1], function(tr)
       !setequal(TipLabels(tr), labels), logical(1)))) {
  stop("`trees` must all share the same tip labels.")
}
```
(Pattern: reuse `.PrepareTrees()` approach once that function gains the same
guard — it currently lacks it too, so the same fix is needed there for the
FACT-family methods.)

---

## F2 — `buildParentDepth` preorder assumption undocumented and unguarded in rstar.cpp — MEDIUM

**Where:** `src/rstar.cpp:63-74` (`buildParentDepth`) and the loop at
`rStarConsensus:270-305` that calls it.

**What:**
`buildParentDepth` computes depths in a single forward pass:
```cpp
depth[c] = depth[p] + 1;
```
This is correct only when the edge matrix is in Preorder (parent appears before
child).  The sister file `local_consensus.cpp` has an explicit comment at its
identical function (line 27): *"Since ape::Preorder() is guaranteed, a single
forward pass suffices."*  `rstar.cpp`'s copy has no such comment.

More critically, `rStarConsensus` is reachable as `ConsTree:::rStarConsensus`
and can be called with an arbitrary (non-Preorder) edge list.  If parent `p`
appears after child `c` in the edge list, `depth[c] = 0 + 1 = 1` (wrong).
All LCA depths computed via the precomputed `D` matrix are then wrong.  A
triplet whose correct resolution is `ab|c` may be tallied as `ac|b`, corrupting
`R_maj` and silently returning a wrong consensus.

The R wrapper (`R/rstar.R:108`) **does** call `Preorder()`, so normal API use
is safe.  The bug surfaces only via the internal C++ entry point.

**Example path to wrong result:**
```r
# Edge matrix with parent AFTER child (non-Preorder)
e <- matrix(c(4,1, 4,2, 5,3, 5,4), ncol=2, byrow=TRUE)  # root=5, non-Preorder
ConsTree:::rStarConsensus(list(e), 3L)   # produces wrong depths → wrong triplet
```

**Fix (defence-in-depth):** add a comment matching `local_consensus.cpp`'s
wording:
```cpp
// Preorder is guaranteed by the R wrapper (Preorder() is called before
// extracting edge matrices); a single forward pass suffices.
```
Optionally add a debug-mode assertion that `depth[p]` > 0 whenever p is not the
root (i.e., `parent[p] != -1`).

---

## F3 — Dead `seen` vector in `assembleRStar` — MINOR

**Where:** `src/rstar.cpp:231`

**What:**
```cpp
std::vector<char> seen(nNode, 0);   // ← allocated, never read or written
while (!stack.empty()) {
  int u = stack.back(); stack.pop_back();
  order.push_back(u);
  for (int c : children[u]) stack.push_back(c);
}
```
`seen` is declared alongside the post-order traversal loop but is never
referenced after initialization.  The iterative DFS does not need a visited
guard (the tree is a DAG, every node is reachable from the root exactly once).
The allocation is harmless (bounded by `nNode` ≤ 2n+1 ≤ 401 bytes at n=200),
but the dead variable is misleading and would confuse a future reader trying to
understand why it exists.

**Fix:** delete the line.

---

## F4 — No n=3 test for `RStar()` — COVERAGE

**Where:** `tests/testthat/test-rstar.R`

**What:**
All shipped RStar tests use n=4..30 leaves.  For n=3 the algorithm exercises
unique code paths:
- Only one 3-taxon subset exists, so the triplet loop runs exactly once.
- The threshold loop in `assembleRStar` runs exactly **one** iteration (theta
  goes 1..n-2 = 1).
- Two distinct cases should be tested:
  1. **Resolved:** both trees agree on the single triplet → one strong cluster of
     size 2 → a fully resolved tree `((a,b),c)`.
  2. **Tied:** trees split 1-1 → no plurality winner → no strong cluster → star
     `(a,b,c)`.

Neither case is currently exercised.  A bug in the threshold-loop boundary
(`theta <= n-2`) would not be caught by the existing suite.

**Suggested additions:**
```r
test_that("RStar() resolves a 3-leaf set when both trees agree", {
  both_ab <- c(ape::read.tree("((a,b),c);"), ape::read.tree("((a,b),c);"))
  expect_setequal(cladeSet(RStar(both_ab)), "a,b")
})
test_that("RStar() returns a star for a tied 3-leaf triple", {
  tied <- c(ape::read.tree("((a,b),c);"), ape::read.tree("((a,c),b);"))
  expect_length(cladeSet(RStar(tied)), 0L)
})
```

---

## F5 — No test for mismatched-taxon-set input — COVERAGE

**Where:** `tests/testthat/test-rstar.R`, `tests/testthat/test-local.R`

**What:**
Neither test file has a case where the input list contains trees with different
tip labels.  Given that F1 identifies this as a silent wrong-result path (before
the fix), a regression test is needed to:
1. Confirm the function errors after F1 is fixed, and
2. Prevent future refactors from silently reverting to the broken behaviour.

**Suggested addition (after F1 fix is applied):**
```r
test_that("RStar() errors when trees have different tip labels", {
  t1 <- ape::read.tree(text = "(((a,b),c),d);")
  t2 <- ape::read.tree(text = "(((a,b),c),X);")
  expect_error(RStar(list(t1, t2)), "tip labels")
})
test_that("Local() errors when trees have different tip labels", {
  t1 <- ape::read.tree(text = "((a,b),(c,d));")
  t2 <- ape::read.tree(text = "((a,b),(c,X));")
  expect_error(Local(list(t1, t2)), "tip labels")
})
```

---

## Algorithm correctness: passed

The following aspects were checked adversarially and found to be correct:

- **Triplet identification** (`rstar.cpp:292-304`): the three if/else-if branches
  correctly cover all four cases (fan, `ij|k`, `ik|j`, `jk|i`).  Any two leaves
  in a rooted tree always have equal or ordered LCA depths with respect to a
  third; the fan condition (all three equal) is caught first, and the remaining
  three are mutually exclusive and exhaustive.  Verified for binary and
  multifurcating inputs.

- **`triIdx` convention**: the encoding `triIdx(a,b,c,n) = (a*n+b)*n+c` with
  close pair `{a,b}` (a<b) and outgroup c is used consistently in the tally
  loop, plurality selection loop, `inRmaj`, and the similarity computation.  The
  ordering constraint a<b is guaranteed at every call site.

- **Plurality selection** (`rstar.cpp:308-326`): ties (nMax > 1) and all-fan
  cases (mx == 0) both correctly zero out all three cells, leaving those triples
  unresolved.  Strict plurality (nMax == 1) is correctly enforced.

- **Laminarity of strong clusters** (`rstar.cpp:170-188`): the defensive
  laminarity check is correct and unreachable for valid input by Jansson et al.
  (2016) Lemma 1.1.

- **Tree assembly** (`rstar.cpp:192-253`): the parent-finding loop correctly
  identifies the smallest strictly-containing node for each non-root node.  The
  rootIdx (all-leaves node) is guaranteed to contain every other node, so `best`
  is always found.  The iterative post-order traversal is correct (DFS pre-order
  reversed).

- **`lcaDepth` safety**: the climbing loop terminates because depth is
  non-negative and every climbing step strictly decreases depth; it never
  accesses `parent[root] = -1` because `depth[root] = 0` prevents the
  condition `depth[u] > depth[v]` from being true at the root.

- **Strong-cluster assembly edge cases**: n=3 (single theta iteration), all-fan
  input (s(a,b)=0 → no candidates → star), and single-resolution input (one
  strong cluster → resolved tree) all produce correct output by trace-through.

- **Local consensus DP termination**: the backtracking sign convention in
  `buildNewick` terminates correctly; the case where `dpBT` returns 0 (the
  initialised value) is never reached in `buildNewick` because the DP always
  sets `dpBT[Lbitmask][Dbitmask]` for all Dbitmask with popcount ≥ 2 (proved by
  induction from j=2), and the single-bit case is handled as a tail after the
  loop exits.

- **`nextBitPerm` (Gosper's hack)**: correct for all inputs in [1, upperBm).
  The loop bound `subX < upperBm` prevents wraparound.

- **MinRLC vs MinILC cost term** (`local_consensus.cpp:308-310, 347-352`):
  the `comb2(cc) * (i - cc)` induced-local cost term is applied only for
  `minrs == false`, correctly distinguishing the two variants.

---

## Summary

| ID | Severity | Description |
|----|----------|-------------|
| F1 | MAJOR    | No cross-tree taxon validation; mismatched taxa silently corrupt triplet tally |
| F2 | MEDIUM   | `buildParentDepth` preorder assumption unguarded/undocumented in rstar.cpp |
| F3 | MINOR    | Dead `seen` vector in `assembleRStar` (rstar.cpp:231) |
| F4 | COVERAGE | No n=3 test for RStar() |
| F5 | COVERAGE | No test for mismatched-taxon-set input in RStar() or Local() |

Core algorithm (triplet tally, plurality selection, strong-cluster assembly,
Local DP, Newick generation): **correct**.  Ship after addressing F1 and F2;
F3 and the coverage gaps can follow in a hardening pass.
