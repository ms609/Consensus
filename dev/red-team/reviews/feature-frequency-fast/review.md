# Red-Team Review: Frequency Consensus (FDCT port)

**Slot:** 3 of 6  
**Date:** 2026-06-02  
**Reviewer:** adversarial (automated red-team)  
**Verdict:** CONDITIONAL SHIP

## Scope

Files reviewed:
- `src/cons_frequency.cpp` — FDCT port, ~50 KB
- `R/selection.R` — `Frequency()` R wrapper and `.PrepareTrees()` helper
- `tests/testthat/test-frequency.R` — shipped test suite
- `dev/red-team/reviews/feature-frequency-fast/repro-01-bounds-battery.R`
- `dev/red-team/reviews/feature-frequency-fast/repro-02-stress.R`

---

## FINDING F1: Missing tip-label validation — silent UB on mismatched taxa sets

**Severity:** MAJOR  
**Location:** `R/selection.R:24–30` (`.PrepareTrees`)

`.PrepareTrees()` records `labels <- TipLabels(trees[[1]])` and passes them to
`RenumberTips(tr, labels)` for every subsequent tree without checking that the
sets match.  If `trees[[i]]` has different, missing, or extra labels,
`RenumberTips` silently maps them to NAs or out-of-range integers.  The C++
`filter_clusters_nlogn` then indexes `taxa_to_leaf_map` with an out-of-range
value, dereferencing a NULL pointer.

**Reproduction:**
```r
t1 <- ape::read.tree(text = "((a,b),(c,d));")
t2 <- ape::read.tree(text = "((a,b),(c,e));")   # 'e' instead of 'd'
Frequency(list(t1, t2))                          # undefined behaviour in C++
```

**Fix:** Add after recording `labels` in `.PrepareTrees()`:
```r
for (i in seq_along(trees)) {
  if (!setequal(TipLabels(trees[[i]]), labels))
    stop("all trees must have the same tip labels")
}
```

---

## FINDING F2: Stale radix state drives the k>>n compression branch

**Severity:** MINOR  
**Location:** `src/cons_frequency.cpp:1204`

`filter()` tests `if (radix->k > radix->n)` to choose counting-sort vs.
quicksort for weight compression.  `radix->k` and `radix->n` are **not reset**
before this check — they hold values from the previous call to a radix
operation.  The intended guard is "current max weight > current node count" but
the code compares stale counts from the last iteration.

The failure mode is a performance regression (wrong sort chosen), not a
correctness error, because the counting sort gracefully handles
`tails[0..k_observed]` even when `k` is an overestimate.  However, if
`radix->k` is 0 on the first iteration, the large-key path is never taken even
when `k >> n`.

**Fix:** Replace `if (radix->k > radix->n)` with an explicit comparison of
the local `k` parameter (the FreqDiff weight bound) against
`tree1->get_nodes_num() + tree2->get_nodes_num()`, or call `radix->clear()`
immediately before the check.

---

## FINDING F3: O(tree-depth) recursion — stack overflow on large caterpillars

**Severity:** MAJOR  
**Location:** `src/cons_frequency.cpp:225` (`fix_tree_supp`), `:429`
(`eulerian_walk`), `:911` (`compute_m`), `:1370` (`newickInto`)

All four functions recurse to a depth equal to the input tree depth.  For a
caterpillar on *n* taxa the depth is *n*.  At *n* ≈ 10 000–50 000 (common in
empirical phylogenomics) this exhausts the default system stack (1–8 MB,
~10 000–50 000 frames).

The shipped tests exercise *n* ≤ 256; the stress battery (`repro-02`) caps at
*n* = 256 as well.  No test ever reaches a depth that would trigger this.

**Reproduction:**
```r
n <- 20000
t <- ape::stree(n, type = "left",  tip.label = paste0("t", seq_len(n)))
u <- ape::stree(n, type = "right", tip.label = paste0("t", seq_len(n)))
Frequency(list(t, u))   # stack overflow in fix_tree_supp or eulerian_walk
```

**Fix:** Convert the four functions to iterative implementations using an
explicit `std::stack`.  `fix_tree_supp` can be replaced with a preorder BFS
loop.  `eulerian_walk` can iterate the preorder node list directly (node IDs
already assigned in preorder by `fix_tree`).  `newickInto` can be replaced with
an iterative string builder.

---

## FINDING F4: NA entries in tree list bypass NULL filter

**Severity:** MINOR  
**Location:** `R/selection.R:19`

`.PrepareTrees()` removes `NULL` elements with
`trees[!vapply(trees, is.null, logical(1))]` but lets `NA` pass through.
An `NA` element is forwarded to `.FactEdges()` → `RenumberTips(NA, labels)`,
producing nonsensical output or a crash.

**Reproduction:**
```r
t <- ape::read.tree(text = "((a,b),(c,d));")
Frequency(list(t, NA, t))
```

**Fix:**
```r
trees <- Filter(function(x) inherits(x, "phylo"), trees)
```
(This also catches `NULL` and any other non-tree object in one pass.)

---

## FINDING F5: Stale scratch-array bounds — unproven but empirically validated

**Severity:** INFO  
**Location:** `src/cons_frequency.cpp:838–843`

The code explicitly acknowledges that the `2n`/`3n` scratch-array sizes and the
radix capacities (`5n` adds, key ≤ `k`) are carried verbatim from the upstream
raw-array allocations and "are not proven closed-form here, but hold empirically
over a wide adversarial battery."  The centroid-path invariant (heavy child ≥
half the leaves) does bound contracted-tree sizes to ≤ 2n − 3, supporting the
`2n` buffers, but no proof comment exists in the code.

**Fix:** Add an `assert(contracted->get_nodes_num() < 3 * n)` guarded by
`#ifndef NDEBUG`, and a comment deriving the bound from the centroid-path
property.

---

## Findings confirmed safe (no action needed)

| Finding | Description | Outcome |
|---------|-------------|---------|
| F6 | `to_del`/`to_del_ti` buffer overflow after `merge_trees` | **SAFE** — merged tree ≤ 2n−1 nodes |
| F7 | Cascade delete: parent+child both marked, use-after-free | **SAFE** — preorder processing updates child parent pointer before child is visited |
| F8 | Strict-plurality deletion condition at ties | **CORRECT** — `>=` at line 1185 correctly deletes on ties, yielding the star topology as required |
| F9 | `filter()` weight-compression / `origw_to_w` interaction | **FRAGILE but safe** — `orig_w` always pre-compression, map populated before sub-contracted queries; needs a clarifying comment |

---

## Coverage Gaps

The shipped test suite covers the lattice invariant (Majority ≤ Frequency ≤
Greedy), all-identical trees, a 2-vs-1 frequency scenario, the 2-tree
all-incompatible case, and a coarse n=40 scalability check.

**Missing test cases:**

1. **k ≥ 3 explicit frequency counts** — split with freq 3 vs. conflicting
   split with freq 2 in k=5; 3-way incompatibility where one wins strictly.
2. **Star topology consensus** — all-unresolved inputs → all-unresolved output.
3. **Caterpillar vs. caterpillar (opposite)** — n=8; exercises the centroid path
   at maximum depth.
4. **Mixed resolution** — some fully-resolved, some star topologies in the same
   run.
5. **Tie at k/2** — k=4 trees, split in exactly 2; must *not* appear (strict
   plurality, tied = deleted).
6. **k=2, one conflicting clade** — minimal frequency-difference case.
7. **k >> n (weight-compression branch)** — the quicksort path in `filter()` is
   exercised by the battery scripts but not by any shipped test.
8. **n=4 minimum** — smallest input that does not short-circuit.
9. **Non-symmetric incompatibility** — split A in k−1 of k trees; k−1 distinct
   single-occurrence conflicts; A must survive.
10. **Mismatched tip labels (error path)** — exercises the missing F1 guard.
11. **NA in tree list (error path)** — exercises the missing F4 guard.
12. **Large-n caterpillar stress** — n=5 000 caterpillar to catch the F3 stack
    overflow before it reaches production users.

---

## Verdict

**CONDITIONAL SHIP**

The algorithmic core is a faithful port of FDCT; the frequency-difference
property, rooting convention, lattice invariants, and memory layout are all
correct for typical inputs.  Two issues require resolution before a general
release:

1. **F1 (MAJOR, R-only fix):** Add taxa-set validation to `.PrepareTrees()`.
   One-line fix; without it a common user mistake causes undefined behaviour in
   C++.

2. **F3 (MAJOR, C++ fix):** Iterativise the four O(depth)-recursive functions.
   Without this, Frequency crashes on caterpillar/pectinate inputs with
   *n* > ~10 000, a real risk for users with large empirical datasets.

F2 and F4 are robustness nits with one-line fixes.  F5 is a documentation gap
that becomes a latent risk only if the contraction/merge logic is modified.
