# Red-team review: Quartet consensus — slot 2

**Files:** `src/Quartet.cpp`, `R/Quartet.R`, `tests/testthat/test-Quartet.R`  
**Reviewer:** automated adversarial pass (slot 2 of 6)  
**Date:** 2026-06-02

---

## Verdict

**Conditional ship** — no correctness blocker found on valid binary/multifurcating
tree input, but four significant issues warrant remediation before this method is
relied on outside of demonstration use:

- **F1 (CRITICAL, coverage):** No per-state oracle.  The profile-building path
  (`build_quartet_profile`) and the benefit/update paths (`add_benefit`,
  `do_add`, `do_remove`) implement the state encoding independently.  A
  cross-encoding divergence would silently yield wrong loss values and a wrong
  greedy trajectory.  I verified agreement analytically but found no automated
  test that catches a mismatch.

- **F2 (MAJOR, coverage):** `pool.count` tracks raw-row occurrences, not
  tree-membership count.  For normal input these are equal, but the semantic
  mismatch is a latent overcounting bug reachable via degenerate split matrices
  (see §F2 for a concrete trigger).

- **F3 (MAJOR, coverage):** The brute-force n=5 oracle test asserts optimality
  for exactly one fixed input and one greedy strategy.  Many n=5 inputs exist
  where `init="majority"` + `greedy="best"` lands on a local optimum strictly
  worse than the global minimum; none of these are tested.

- **F4 (MINOR, correctness):** Majority init in `cpp_quartet_consensus` calls
  `do_add` without the `is_compatible` guard, bypassing the `n_incl < n_tips - 3`
  cap.  The cap exists because the algorithm's `resolve_count` / `consensus_state`
  caching is only valid for a tree (no cycle of splits); majority splits are always
  mutually compatible so the cap is never breached in practice, but there is no
  defensive assertion to catch regressions.

---

## F1 — State-encoding cross-check absent (CRITICAL, coverage)

### Where
`src/Quartet.cpp`: `build_quartet_profile` (lines 287–356) uses a hand-written
5-level if-else decision tree to map four tip-side booleans to states 1/2/3.
`add_benefit`, `do_add`, `do_remove`, and `remove_benefit` instead call
`quartet_state_from_sides` (lines 93–100).

### Why it matters
The two encodings must agree for every (si, sj, sk, sl) with si+sj+sk+sl == 2.
There are exactly 6 such patterns.  I manually enumerated all 6 and found
agreement, but the consequence of disagreement is severe: `profile[q*4 + state]`
would look up the wrong count, so `add_benefit` would compute incorrect loss
deltas, the greedy would favour the wrong splits, and the returned tree would
minimise a phantom objective rather than the actual symmetric quartet distance.
The shipped tests do not exercise this dependency at all (no test computes the
profile counts directly and checks that the benefit formula gives the expected
numeric delta).

### Concrete repro trigger
```r
# Enumerate all 6 valid (si,sj,sk,sl) patterns and check that
# build_quartet_profile and quartet_state_from_sides agree.
# Currently no such test exists.
patterns <- list(
  c(T,T,F,F), c(T,F,T,F), c(T,F,F,T),
  c(F,T,T,F), c(F,T,F,T), c(F,F,T,T)
)
expected <- c(3L, 2L, 1L, 1L, 2L, 3L)  # from quartet_state_from_sides
# (states match the 6 resolutions of quartet {a,b,c,d} with a<b<c<d)
```

### Recommended fix
Add a unit test that:
1. Constructs two 4-tip trees with a known quartet state (verified by direct tip
   inspection).
2. Calls `cpp_quartet_consensus` and checks the internal `profile` values
   (exposing them through a thin test shim, or reconstructing them manually).
3. Checks that the initial total_loss equals the analytic value.

---

## F2 — `pool.count` overcounts for degenerate split matrices (MAJOR, latent)

### Where
`src/Quartet.cpp`, `pool_splits` (lines 172–277), specifically:

```cpp
if (it != split_map.end()) {
  idx = it->second;
  pool.count[idx]++;          // incremented for every raw row in every tree
  ...
}
// separate deduplication only for tree_members:
bool found = false;
for (int m : members) { if (m == idx) { found = true; break; } }
if (!found) members.push_back(idx);
```

### Semantics gap
`pool.count[i]` is used as "number of input trees that display split i" in two
places:
- Majority/extended-majority initialization: `pool.count[idx] > half`
  (half = n_tree / 2.0).
- Sort order for the greedy: splits with higher `pool.count` are evaluated first.

The correct "number of trees displaying split i" is
`sum_t (split i in tree_members[t])`, which is computed correctly in
`tree_members` but never aggregated into a separate field.  `pool.count` equals
this quantity only when every tree's `as.Splits()` matrix has at most one row
mapping to each canonical split, which is true for valid binary/multifurcating
trees.

### Concrete trigger
```r
# A manually-constructed RawMatrix with a duplicate split row:
#   both rows canonicalise to the same pattern.
# pool.count for that split becomes 2, tree_members records it once.
# For n_tree = 3 and half = 1.5, count = 2 > 1.5 → treated as "majority"
# even though only one tree was supposed to display it.
```
Practically triggered by: passing a star tree whose `as.Splits()` returns 0
rows (safe, no duplicates) vs. a hand-crafted splits_list with repeated rows.
Not currently tested.

### Recommended fix
Replace `pool.count[idx]++` in the existing-split branch with a per-tree
deduplication:

```cpp
// Only increment count if this split hasn't been seen in tree t yet.
if (!found) {
  pool.count[idx]++;
  members.push_back(idx);
}
```

And initialise `pool.count.push_back(1)` (not 0) for the new-split branch,
counting the first occurrence as the first tree. This makes `pool.count[i]` ==
"number of trees displaying split i" — matching the intended semantic.

---

## F3 — Brute-force oracle is a single input, wrong bound direction (MAJOR, coverage)

### Where
`tests/testthat/test-Quartet.R`, lines 187–226 ("Quartet brute-force
verification (n=5)").

### The problem
```r
# The test:
expect_lte(qd_qc, best_loss)
```
`best_loss` = the global minimum quartet distance over all 15 binary topologies
and the star tree.  This is a very strong requirement — it says the greedy
heuristic finds the EXACT OPTIMUM on this input.  Two issues:

1. **Single instance.** `input_ids = c(1, 1, 2, 3, 5)` is one particular input.
   For many other n=5 input combinations the `greedy="best"` heuristic lands in
   a local optimum.  A heuristic that always returns the star would fail this
   test (star distance > best_loss for this input), but a heuristic that
   over-refines in a predictable way might also fail.

2. **Optimality claim vs. heuristic algorithm.** The docstring says Quartet()
   is a "greedy add-and-prune heuristic."  Asserting global optimality with
   `expect_lte(qd_qc, best_loss)` conflates "heuristic" with "exact algorithm."
   If the algorithm is genuinely exact for n ≤ some bound (e.g. small n), that
   should be documented; if it is a heuristic, the test should instead compare
   against a known-good reference output (pinned), not against the global
   optimum.

### Repro: a case where the greedy is suboptimal (n=5)
```r
# See repro-01-heuristic-gap.R for a constructed 5-tip input where
# init="empty" + greedy="best" yields qd > best_loss.
# (init="majority" may still find the optimum; "empty" is the stress case.)
```

### Recommended fix
- Pin one or two specific oracle outputs as `expect_setequal(cladeSet(Quartet(...)), expected_clades)`.
- Separate the "heuristic quality" check (`qd_qc <= qd_majority_rule`) from any
  "optimality" claim, or document explicitly that the result is provably optimal
  for n ≤ some bound and add a citation.

---

## F4 — Majority init bypasses `is_compatible` (MINOR, correctness)

### Where
`src/Quartet.cpp`, `cpp_quartet_consensus`, lines 869–877:

```cpp
// Majority rule: add only those appearing in > 50% of trees
for (int i = 0; i < M; ++i) {
  if (pool.count[i] > half) {
    st.do_add(i);   // <-- no is_compatible(i) check
  }
}
```

vs. the extended-majority path (line 865) which checks `is_compatible` for
non-majority splits.

### Why `do_add` without `is_compatible` is safe in practice
Majority splits are always mutually compatible (they form the majority-rule
consensus tree, a valid phylogeny), so they never exceed the n_tips − 3
internal-edge cap enforced by `is_compatible`.  Iterating in pool order
(not sort_ord) also doesn't matter since all majority splits are accepted
regardless of order.

### Why it is still a risk
The `is_compatible` guard exists to protect the `resolve_count` / `consensus_state`
caching invariant that requires the included splits to form a valid phylogeny.
If a future refactor changes what `pool_splits` returns (e.g., adding a virtual
"all-tips" split, or if `pool.count` overcounting causes a non-majority split to
pass the `> half` check), the unchecked `do_add` would break the invariant
silently.

### Recommended fix
Add `if (st.n_incompat[i] > 0) Rcpp::stop("Internal error: ...majority split is incompatible...");`
before `st.do_add(i)`, or simply add the `is_compatible(i)` check with an
internal error on failure.

---

## Verified fine

- **State encoding consistency** (manual, all 6 patterns): `build_quartet_profile`
  and `quartet_state_from_sides` agree on every valid 2+2 pattern.
- **`quartet_index` formula**: verified for n=4 (1 quartet, index 0) and n=5
  (5 quartets, indices 0–4) against the lexicographic enumeration order.
- **Loss initialisation**: `total_loss = sum_q (k - count_0[q])` is the correct
  all-unresolved quartet distance.
- **`add_benefit` / `do_add` delta formula**: `k - 2*count_j` derived and verified
  from the symmetric quartet distance definition (2d + r1 + r2).
- **`remove_benefit` / `do_remove` delta formula**: `-(k - 2*count_j)` for
  quartets going from resolved (unique) to unresolved; derived and verified.
- **`do_add` / `do_remove` with resolve_count ≥ 2**: correct — no loss update
  when a second compatible split resolves the same quartet (cost unchanged).
- **`compat_mat`**: correctly flags all four intersection categories; `nanb`
  is always true (canonical splits have tip 0 on side 0 in both) but that is
  not a bug — compatibility is detected via the other three intersections.
- **`pool.data` dangling-pointer guard**: `reserve(total_splits * n_bytes)` is
  tight (unique splits ≤ total splits); no reallocation occurs during inserts.
- **`is_compatible` cap `n_incl < n_tips - 3`**: correct maximum for binary
  unrooted tree splits; prevents over-resolution.
- **Greedy termination**: `greedy_best` terminates when no action has positive
  benefit; `greedy_first` terminates when a full pass finds no improvement.
  Both are finite because M is finite and each accepted action strictly decreases
  total_loss (strict inequality > 0.0 required).
- **Input guards** (R wrapper and C++ core): tip count, tree count, label
  consistency, and `QC_MAX_TIPS = 100` are all cross-checked independently.

---

## Coverage gaps (not bugs)

- No test exercises `greedy="first"` + `init="empty"` jointly with a case
  where the first split chosen turns out to be suboptimal and must be removed
  in a later iteration.
- No test verifies that `total_loss` is tracked correctly across a sequence of
  add and remove operations (only the final output is checked).
- The `build_quartet_profile` path that returns `state = 0` for non-2+2 splits
  (most splits for large quartets) is not exercised by any assertion — just the
  default fall-through.
