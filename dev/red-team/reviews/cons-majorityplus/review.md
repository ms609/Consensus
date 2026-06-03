# Red-team review: cons_majorityplus.cpp
**Date:** 2026-06-03
**Verdict:** SHIP

## Scope

`src/cons_majorityplus.cpp` and the `MajorityPlus()` R wrapper in `R/selection.R`.
Compared line-by-line against the FACT reference source in
`dev/oracle/fact-src/majorityplus.cpp` (plus `majority.cpp:169` for `majContract`).
All claims below are backed by both static analysis and a live oracle run.

## Oracle / empirical result

Build identity confirmed by `assertConsTreeBuild()`.
`Rscript.exe dev/oracle/check-oracle.R` on this build:

- All cross-validation rows (n9/n10/conflict-n8, rooted=0 and rooted=1): **MATCH**
- `majorityPlus at n > 60` hard-assertion block (n=80 and n=137): idempotent **TRUE**, FACT-exact **MATCH**, non-trivial split count (51 and 95 respectively) — `stopifnot(all(mpPass))` passed.

`Rscript.exe dev/red-team/reviews/majorityplus-fact/repro-01-strict-boundary.R`:
**PASS** — strict `>` boundary is FACT-exact.

Additional edge-case checks run inline:
- Two trees sharing a split but each contradicting the other: star (NSplits = 0). PASS.
- Three identical trees: idempotent. PASS.
- Majority split (4 of 5 trees display, 1 contradicts): split survives (NSplits >= 1). PASS.

## Hunt items

### 1. Display vs. contradiction counts — three-way distinction

**Correct.**

In `updateCounter`, for each cluster of A queried against B, the outcome is:

| Condition | Effect on `goodLabel` | Semantics |
|---|---|---|
| POS[b] - POS[a] + 1 != A.size (range not contiguous) | `--` | B contradicts cluster |
| Path query: deepest anchor mismatch | `--` | B contradicts cluster |
| `a == b` (same node; cluster is a subset of B's subtree) | `++` | B displays cluster |
| `a != b` but `B.parent[a] == B.parent[b]` | no change | B compatible but not displaying |
| `a != b` and `B.parent[a] != B.parent[b]` | `--` | B contradicts cluster |

The no-change case (compatible, not displayed) contributes 0 — exactly the required
three-way split. The port uses a `decided` flag in place of FACT's `goto end`; the
reachable states are identical and the comparison table above matches the FACT source
comment "Cannot be inserted, so we minus" at the contradiction branches.

### 2. Incremental merge — double-counting

**Correct.**

`majorityPlusConsensus` runs two passes:

1. **Merge pass** (lines 419-426): builds the union of all candidate clusters via
   `majorityPlusMerge`. The running counter from this pass is **discarded** at line 428
   (`ret.goodLabel[i] = 0` for all i).
2. **Recount pass** (line 429): calls `updateCounter(ret, T[i])` for every tree over
   the final merged structure, starting from 0. Each tree contributes exactly once to
   each cluster's tally.

There is no double-counting. The reset between passes is explicit in both FACT and the
port at the same logical location.

### 3. Boundary condition: strict > vs >=

**Correct — strict `>`.**

`majContract` keeps node `i` when `X.goodLabel[i] > 0` (lines 53 and 57). After the
full recount, `goodLabel[i] = display_count - contradict_count`. So:

- `goodLabel > 0` iff `display > contradict` (strict inequality). CORRECT.
- `goodLabel == 0` (display == contradict) → node dropped. A `>=` would keep it.

Verified at the source, and independently confirmed by repro-01 (2 display, 2
contradict → star; 2 display, 1 contradict → split kept).

### 4. Path-query correctness

**Correct.**

The left-path `L[minL]` and right-path `R[maxL]` arrays store, for each Day-label
`k`, all nodes whose subtree starts at (respectively ends at) label `k`, ordered by
increasing depth. The query for a cluster [a, b]:

1. Takes the deepest node from `L[a]` and the deepest from `R[b]`.
2. Compares their depths to decide which path to walk.
3. Steps up the shorter path to the LCA candidate, indexing into the other path by the
   depth difference.

This is Day's O(1)-per-cluster compatibility test. The port is a direct transcription
of the FACT code; the left/right symmetry is handled by the `DEPTH[L[a].back()] >=
DEPTH[R[b].back()]` branch that selects which path to anchor first — matching FACT
exactly.

### 5. n > 60 multi-word path

**Not applicable to this method; confirmed at scale.**

The file header and `check-oracle.R` (lines 191-194) both document that
`cons_majorityplus.cpp` contains **no bit-packing** (no BUCKET_SIZE). Compatibility
is decided purely by Day's leaf relabelling plus integer arithmetic on `POS`, `DEPTH`,
`minL`, `maxL` — all of which are `int` and scale to any n without word-boundary
arithmetic. The oracle run at n=80 and n=137 exercise the same code paths as n <= 60,
and both matched FACT exactly with non-trivial output (51 and 95 splits respectively).

### 6. Edge cases

All four cases verified:

| Case | Result | Expected |
|---|---|---|
| All trees identical | NSplits = 5 (= input tree) | idempotent |
| Two trees contradicting every split (independent random trees) | NSplits = 0 | star |
| Two trees: one shared split, each contradicts the other | NSplits = 0 | star (2 display, 2 contradict: not >) |
| Majority split (4/5 display, 1 contradicts) | NSplits >= 1 | split survives |

The majority-survival guarantee holds because a split in more than n/2 trees has
display > n/2 and contradict < n/2 (a tree can only contradict OR not-display, never
both), so display - contradict > 0 always.

### 7. Rooting

**Correct.**

The R wrapper `MajorityPlus()` calls `.FactEdges()` which roots every input at
`labels[[1]]` (taxon 1) via `RootTree(tr, labels[[1]])` and returns the preorder edge
matrix. The C++ entry point receives consistently rooted inputs and produces a Newick
string with integer labels. `.RootLikeFirst()` then roots the consensus to match the
first input tree's root group, the same convention used by all other fast-path
methods and `TreeTools::Consensus()`.

The comment in the C++ (lines 39-42) correctly warns callers not to pass arbitrarily
rooted trees; the R wrapper satisfies this precondition for all calls through the
public API.

## Minor observations (no action required)

**`majContract` leaf-safety invariant.** The contraction loop at lines 61-64 uses
`label[i]` (initialised to -1) inside the `X.leaf[i] > 0` branch without a separate
`goodLabel > 0` guard. This is safe because after the second-pass `updateCounter`,
every leaf node always receives `++goodLabel` from every tree (a singleton cluster
[a, a] always satisfies `a == b`). Confirmed by examining `precompute()` in
`src/fact_tree.cpp`: `size[leaf_node] = 1`, so the contiguous-range check always
passes for leaves, and `++goodLabel[leaf]` fires for every tree. The same latent
pattern exists identically in FACT. Both rely on this invariant silently.

**Dead `sum` array in FACT.** FACT's `majorityPlusConsensus` declares and zeroes a
`sum[numTaxas+5]` array that is never read. The port correctly omits it.

**`majContract` root safety.** If `label[X.root] == -1` (root dropped), `ret.root`
would be -1, causing downstream UB. This cannot happen: the root spans all N taxa and
is always displayed by every tree, so its goodLabel after the recount is numTrees > 0.
Both FACT and the port rely on this invariant silently.

## Test coverage

`tests/testthat/test-selection.R` covers: identical-tree idempotence,
majority-contained-in-majorityPlus lattice invariant, single-tree/fewer-than-4-leaves
trivial path, bare `phylo` short-circuit, and non-list input rejection. The strict
boundary case is covered by the repro script in
`dev/red-team/reviews/majorityplus-fact/repro-01-strict-boundary.R`.

**Gap:** no testthat test for the star-collapse case (all input trees contradicting
every non-trivial split). Low risk given oracle coverage, but worth adding to
`tests/testthat/test-selection.R` for regression protection in CI.

## Summary

The port is a faithful, structurally correct translation of FACT's
`majorityPlusConsensus` / `updateCounter` / `majorityPlusMerge` / `majContract`. Every
algorithmic invariant checked by static analysis is confirmed by the oracle. No bugs found.
