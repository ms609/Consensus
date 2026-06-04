# Frequency() depth-robustness fix — diagnosis, fix, validation

Follow-up to `audit-efficiency.md` finding **F1** (deep-tree blow-up). Branch
`feature/frequency-depth`. Touches `src/cons_frequency.cpp` only (algorithm
output — the unique frequency-difference split set — is unchanged).

## Diagnosis (measured, not assumed)

The reported symptom was a *catchable* `std::bad_alloc` at n≈30000 on opposite
caterpillars — heap, not a stack overflow. A clean-build sweep (each n in a
fresh process, peak working set via the Windows API; `diag-depth.R`) showed the
cost growing **quadratically in both time and memory**, ruling out the audit's
"O(depth) recursion → stack overflow" hypothesis:

| n | time | peak WS |
|------|--------|---------|
| 2000 | 6.7 s  | 642 MB  |
| 4000 | 27.6 s | 2297 MB (3.6×) |
| 8000 | 118 s  | 8850 MB (3.85×) |

Doubling n quadruples both ⇒ **O(n²)** (≈3.6 GB of ints projected at 30000 ⇒
`bad_alloc`). This is the audit's failure mode **(b)** — per-level O(n) scratch
in an O(depth) recursion — *not* **(a)**, so iterativising the helper recursions
would not have helped.

Phase timing localised all of it to `filter()`. Recursion counters inside
`filter_clusters_nlogn` (per top-level `filter`, at n=1000/2000) pinned the
mechanism exactly:

```
n=1000 : calls=1000 maxdepth=1000 sumNN=500500(≈n²/2) sumContract=998001(≈n²)
n=2000 : calls=2000 maxdepth=2000 sumNN=2001000      sumContract=3996001
```

`maxdepth = n` and `sumPN ≈ 2` per call: the recursion is a **linear chain of
depth n**, and tree2 shrinks by only one leaf per level (ΣNN ≈ n²/2). The
"centroid path" the filter descends is `t1_root → children[0] → … → leaf`, but
on a caterpillar `children[0]` is a **lone leaf**, so the whole (n−1)-leaf clade
becomes a single side branch that is recursed — the exact opposite of the
heavy-path decomposition the O(kn log n) bound requires.

### Root cause

`filter_clusters_nlogn` assumes `children[0]` is the heaviest child (the upstream
comment literally says "construct heavy path of this subtree"). That invariant is
meant to be established by `Tree::reorder()` (heaviest child → position 0), called
on the inputs in `run()`. **But `reorder()` only permutes child pointers and never
renumbers ids, and both the deep-copy `Tree(Tree*)` and `merge_trees` re-add
children in id / `m` order — silently undoing it.** So the heavy-child invariant is
never actually in effect on the trees passed to the filter. This is a latent bug in
the upstream FDCT_new too (its `Tree(Tree*)` re-adds in id order identically); the
reference `freqdiff.exe` is itself O(n²) on caterpillars (n=5000 already overflows
its own radix bound).

## Fix

Re-establish the heavy-child-first invariant on **both** trees at the top of
`filter()` by following `reorder()` with `fix_tree()` (which renumbers preorder so
the heaviest child genuinely lands at id-position 0 and survives the copy):

```cpp
tree1->reorder(); tree1->fix_tree();
tree2->reorder(); tree2->fix_tree();
```

- Done in `filter()` (not once in `run()`) so it also covers the merged
  accumulator `T`, which `merge_trees` rebuilds in `m`-order.
- tree1 heavy-first ⇒ `children[0]` descent is the heavy path ⇒ O(log n) recursion
  depth. tree2 heavy-first ⇒ its `pos_in_parent==0` centroid paths are heavy paths
  ⇒ O(log n) `max_subpath_query`.
- Idempotent on an already heavy-first tree; the run()-level `reorder()` (now dead)
  was removed.

### Output preservation

The path/centroid decomposition is a **pure performance device**: every cluster is
tested exactly once against a contraction that preserves its conflict information,
so the marked-for-deletion set — the unique frequency-difference split set — is
independent of the path choice. Reordering does change the **Newick child-order**
(serialisation), but not the consensus. Confirmed:
- `capture-newick.R`: same Newick **length** for all 16 battery cases, different
  md5 (same topology, reordered children).
- `check-fd-worktree.R`: split sets **exactly** match the boost `freqdiff.exe`
  oracle on all fixtures incl. n80 congruent/incongruent, and on caterpillars.
- (User confirmed split-set identity is the required invariant, consistent with the
  oracle's own `setequal` comparison.)

## Validation (clean `--preclean` build, identity asserted)

1. `build-identity` → **Build identity OK** (ConsTree 0.0.0.9008).
2. Oracle (`check-fd-worktree.R`, worktree-correct — the shipped `check-freqdiff.R`
   prepends the *main* `.agent-cons` and would mask the worktree build): **ALL
   MATCH**, incl. n80 congruent/incongruent and caterpillar n=500/2000.
3. Full `testthat`: **FAIL 0 | PASS 537 | SKIP 1** (skip is an unrelated on-CRAN
   Quartet test). The n=30000 deep guard is flipped on and passes; added
   intermediate n=15000/50000 caterpillars and a deep non-caterpillar (cherry
   ladder, ~16000 tips). One pre-existing fragile assertion (test-frequency.R:69)
   was switched from raw `as.character` to the polarised `splitSet()` helper —
   the consensus is correct (oracle-verified); only its un-polarised string
   representation flipped with the new child-order.
4. Speed (`compare.R Frequency` vs `baseline-2026-06-02.csv`): **no regression
   anywhere**; 4.75×–635× faster (50 trees × 50 leaves independent: 19 s → 0.03 s),
   two cells that timed out now finish; common shallow case unchanged.
5. Depth: opposite-caterpillar n=30000 → **0.5 s / 176 MB** (was `bad_alloc`);
   **n=100000 → 2.0 s / 388 MB**; k=4/k=6 caterpillars and balanced n=50000 fast.

## Out of scope (flagged separately)

`Frequency()` throws `radix sort n overflow` on highly-**incongruent** large inputs
(e.g. two independent `RandomTree`s at n≥5000). Confirmed **pre-existing** (pristine
`main` fails identically) and orthogonal to depth (an incongruence/scratch-bound
issue). Flagged as its own task.
