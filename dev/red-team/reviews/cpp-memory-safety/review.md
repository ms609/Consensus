# Red-team review: C++ memory safety
**Date:** 2026-06-03
**Scope:** All `src/*.cpp` except `RcppExports.cpp`
**Prior reviews excluded:** algorithm correctness findings from
`feature-frequency-fast` (F1–F9) and `feature-loose-fast` reviews.
**Verdict:** CONDITIONAL SHIP

Files examined:
- `src/fact_tree.cpp` / `fact_tree.h`
- `src/cons_greedy.cpp`
- `src/cons_loose.cpp`
- `src/cons_majorityplus.cpp`
- `src/local_consensus.cpp`
- `src/rstar.cpp`
- `src/bhv.cpp`
- `src/Quartet.cpp`
- `src/cons_adams.cpp` (stub only — no code paths)

---

## F1 — Unchecked depth-index into right/left path vectors (MAJOR)

**Files/lines:**
- `src/cons_loose.cpp:210` — `b = R[b][DEPTH[B.idx[b]] - DEPTH[a]]`
- `src/cons_loose.cpp:218` — `a = L[a][DEPTH[B.idx[a]] - DEPTH[b]]`
- `src/cons_majorityplus.cpp:202` — `b = R[b][DEPTH[B.idx[b]] - DEPTH[a]]`
- `src/cons_majorityplus.cpp:210` — `a = L[a][DEPTH[B.idx[a]] - DEPTH[b]]`
- `src/cons_majorityplus.cpp:338` — same pattern in `majorityPlusMerge`
- `src/cons_majorityplus.cpp:342` — same pattern in `majorityPlusMerge`

The left/right path vectors `L[x]` and `R[x]` are populated during the ordered-tree
construction phase: `L[minL[v]]` and `R[maxL[v]]` collect nodes in DFS postorder.
The query phase then indexes these vectors by a depth difference, relying on the
FACT invariant that the depth difference never exceeds the number of nodes on the
path.  There is no bounds check.

**Failure mode:** If an input tree is malformed (e.g., a node with a `minL`/`maxL`
value inconsistent with its ancestors), the computed index `DEPTH[B.idx[b]] -
DEPTH[a]` can exceed `R[b].size() - 1`, producing an out-of-bounds `std::vector`
read (undefined behaviour; in debug builds, a `std::out_of_range` exception from
`at()`; in release, a silent wrong-pointer dereference or SIGSEGV).

The same pattern appears identically in three functions across two files:
`looseMerge` (`cons_loose.cpp`), `updateCounter` and `majorityPlusMerge`
(`cons_majorityplus.cpp`).

**Trigger:** A pathological input such as a tree passed as already-ordered (the
BEFORE/AFTER invariant) but with incorrect FACT node metadata — for example, an
R-side tree whose edge matrix encodes a non-preorder traversal, or has a gap in
node IDs.

**Fix (guard):** Replace bare indexing with a bounds-checked form:
```cpp
// before: b = R[b][DEPTH[B.idx[b]] - DEPTH[a]];
int ridx = DEPTH[B.idx[b]] - DEPTH[a];
if (ridx < 0 || ridx >= (int)R[b].size()) {
  // incompatible cluster — treat as non-compatible
  can[t1.first] = 0; skip = true;
} else {
  b = R[b][ridx];
}
```
Apply the same pattern for the symmetric `L[a]` access and for all three
functions.

**Fix (root-cause):** Add input validation in `buildTreeFromEdge()` that verifies
the edge matrix is preorder (parent always appears before child) and rejects it
with `Rcpp::stop()` otherwise.  This prevents the malformed-input scenario.

---

## F2 — Leaf not unconditionally marked as kept in `majContract` (MINOR)

**File/line:** `src/cons_majorityplus.cpp:57–64`

```cpp
for (int i = 0; i < cnt; ++i) {
    if (X.goodLabel[i] > 0) { label[i] = tmp++; ... }
    if (X.leaf[i] > 0) {
        ret.leaf[label[i]] = X.leaf[i];   // label[i] == -1 if goodLabel <= 0
        ret.idx[X.leaf[i]] = label[i];
    }
}
```

Leaf nodes (`X.leaf[i] > 0`) are assigned into `ret` using `label[i]`, but
`label[i]` is only set if `X.goodLabel[i] > 0`.  If a leaf node somehow had
`goodLabel <= 0`, `label[i]` would be `-1` and `ret.leaf[-1]` would be an
out-of-bounds write (UB; on Linux/Windows the write lands 4 bytes before the
heap-allocated vector buffer — silent data corruption).

**Why it does not fire today:** `updateCounter` increments `goodLabel` for leaf
nodes in every call (a leaf's `a == b` branch always triggers `++goodLabel`), so
after `numTrees` calls all leaf goodLabels are `numTrees > 0`.  The invariant is
maintained implicitly by caller contracts and is not documented or asserted.

**Risk:** Any refactoring that resets `goodLabel` without calling `updateCounter`
first (e.g., a partial-tree extension) would silently corrupt the heap at this
line.

**Fix:** Separate the label-assignment and leaf-copy loops with an assertion
between them:
```cpp
for (int i = 0; i < cnt; ++i) if (X.goodLabel[i] > 0) label[i] = tmp++;
// assert: every leaf node must have been labelled
for (int i = 0; i < cnt; ++i) {
    if (X.leaf[i] > 0) {
        assert(label[i] >= 0);   // catches caller invariant violation early
        ret.leaf[label[i]] = X.leaf[i];
        ret.idx[X.leaf[i]] = label[i];
    }
}
```

---

## F3 — `anc[-1]` double-dereference in `contract()` for non-kept root (MINOR)

**File/line:** `src/cons_loose.cpp:70`

```cpp
if (keep[t1.first]) anc[t1.first] = t1.first;
else anc[t1.first] = anc[anc[t1.first]];
```

`anc` is initialised to `-1` everywhere.  When a non-root, non-kept node `v` is
first visited and its parent slot `anc[v]` is also `-1` (because no kept ancestor
has been seen yet on this path), the expression `anc[anc[t1.first]] = anc[-1]`
reads `anc` at index `-1` — an out-of-bounds `std::vector` read (UB).

**When it can fire:** DFS begins at `A.root`.  If the root is not in the `keep`
set (line 67 guards `t1.first != A.root`, preventing the safe branch for the
root), `anc[root]` stays `-1`.  When the root's first child `v` is also not kept,
the code executes `anc[v] = anc[anc[v]] = anc[-1]`.

**Is the root always kept in practice?** In the current code, `precompute()`
assigns `H[N]` to the root node and `looseConsensus` initialises `good[i] =
(H[i] >= 0)`, so `good[N] = 1` and the root is always in `keep`.  The safety
relies on this indirect chain; there is no explicit `keep[A.root] = 1`.

**Fix:** Add one unconditional line after the keep-flagging loop:
```cpp
keep[A.root] = 1;   // root is always retained; prevents anc[-1] on non-kept root
```

---

## F4 — `lcaDepth` walks past root `-1` on disconnected input (MINOR)

**Files/lines:**
- `src/rstar.cpp:79–82`
- `src/local_consensus.cpp:48–52`

```cpp
while (depth[u] > depth[v]) u = parent[u];
```

`parent[root]` is set to `-1`.  On a disconnected edge matrix (two separate
connected components, missing root edge, or a cycle), `u` and `v` can be in
different subtrees.  The equalisation loop then walks `u` all the way to the root
and executes `u = parent[u] = -1`, after which `depth[-1]` is read — an
out-of-bounds `std::vector` read.

The R wrappers for RLC/ILC and R* do not validate that the edge matrix is a
connected tree.

**Fix:** Add a root guard inside each equalisation loop:
```cpp
while (depth[u] > depth[v]) {
    if (parent[u] < 0)
        Rcpp::stop("lcaDepth: malformed tree (reached root before finding LCA)");
    u = parent[u];
}
```
Or validate connectivity of the edge matrix in the R wrapper before dispatching
to C++.

---

## F5 — `gtp_no_common` front-of-vector queue is O(n^2) (INFO)

**File/line:** `src/bhv.cpp:167,197`

```cpp
QR r = queue.front(); queue.erase(queue.begin());
...
queue.insert(queue.begin(), r2);
queue.insert(queue.begin(), r1);
```

`std::vector::erase(begin())` and `insert(begin(), ...)` are O(n) operations,
making the GTP queue processing O(queue_length^2).  No memory-safety concern: the
element is copied before `erase`, and no out-of-bounds access occurs.  Flagged
because the queue length equals the number of splits in the incompatibility group,
which for large unresolved trees could be substantial.

**Fix:** Replace `std::vector` with `std::deque` for O(1) front push/pop.

---

## F6 — `pool.data` pointer stability (INFO, confirmed safe)

**File/line:** `src/Quartet.cpp:206–212`

`pool.data` is reserved with `total_splits * n_bytes` before any insertions
(line 212).  `split_map` stores raw `unsigned char*` pointers into `pool.data`.
If `data.resize()` ever exceeded the reserved capacity it would move the buffer,
invalidating all stored keys — dangling pointers in `split_map`, causing
`find()` to read freed memory.

**Why it is safe today:** `total_splits` is an exact upper bound (sum of all row
counts across all input matrices, counted at line 209).  Each insertion adds
exactly `n_bytes` bytes.  So `data.size()` can never exceed `total_splits *
n_bytes`, and capacity is never exhausted.

**Residual risk:** A future refactor that processes extra splits outside
`splits_list` would silently invalidate this reserve.

**Fix (defensive):** After all insertions, add:
```cpp
assert(pool.data.size() <= pool.data.capacity());
```

---

## Verified safe (per file)

| File | Issue examined | Outcome |
|------|----------------|---------|
| `fact_tree.cpp:buildTreeFromEdge` | `G[p]`, `parent[c]` indexed by 0-based ape node; array size nNode+5; p,c ≤ nNode-1 | SAFE |
| `fact_tree.cpp:precompute` | iterative via `std::stack`; no recursion | SAFE |
| `fact_tree.cpp:newick` | iterative via `std::stack`; no recursion | SAFE |
| `cons_greedy.cpp:39` | `cluster` sized `(2*numTaxas+5)*LEN`; `Ti.cnt ≤ 2*numTaxas-1` for valid ape trees | SAFE |
| `cons_greedy.cpp:98` | `CountingSort` indexed by run-length ≤ numTrees; sized numTrees+5 | SAFE |
| `cons_greedy.cpp:153` | `Tree newT(numTaxas, ret.cnt+1)` — new node `ret.cnt` within array bound `ret.cnt+6` | SAFE |
| `cons_loose.cpp:57` | `label[i]` always set before use: all leaf nodes have `keep[A.idx[taxon]]=1` | SAFE |
| `cons_loose.cpp:POS[numTaxas+2]` | vector sized `numTaxas+5`; index `numTaxas+2 < numTaxas+5` | SAFE |
| `cons_majorityplus.cpp:goodLabel leaf invariant` | after numTrees `updateCounter` calls, leaf goodLabel = numTrees > 0 | SAFE (fragile, see F2) |
| `rstar.cpp:triIdx` | `(size_t)a * n + b` arithmetic; n≤200; max result 7.9M; no overflow | SAFE |
| `rstar.cpp:D array` | `(size_t)i * n + j`; max index 39999; no overflow | SAFE |
| `local_consensus.cpp:comb2` | n≤20; cc≤20; 20\*19/2=190; no int overflow | SAFE |
| `local_consensus.cpp:1u<<n` | n≤20 enforced at line 265; 1u<<20=0x100000; within `uint` range | SAFE |
| `local_consensus.cpp:dpBT[L][(1<<m)-1]` | `dpBT[L]` sized `1<<m` at line 302; index `(1<<m)-1` is last element | SAFE |
| `local_consensus.cpp:dfs` | recursive on Aho graph with n≤20 vertices; depth ≤ 20 | SAFE |
| `local_consensus.cpp:buildNewick` | recursive; bounded by n≤20 | SAFE |
| `bhv.cpp:bit-indexing` | `c[tip>>6] \|= 1ULL << (tip & 63)`; `tip & 63` always 0..63 | SAFE |
| `bhv.cpp:cap/flow arrays` | sized `N*N`; N = na+nb+2 where na,nb = split counts (typically small) | SAFE for typical inputs |
| `Quartet.cpp:compat[i*M+j]` | M ≤ n*(n-3)/2 ≈ 4850 for n=100; M^2 ≈ 23.5M < INT_MAX | SAFE |
| `Quartet.cpp:quartet_index` | int_fast32_t; max ≈ C(97,4) ≈ 3.5M; fits int32 | SAFE |
| `Quartet.cpp:pool.data stability` | reserved to `total_splits * n_bytes` before insertions; no reallocation | SAFE (see F6) |
| `cons_adams.cpp` | stub body only; no data access | N/A |

---

## Summary by severity

| Finding | File(s) | Severity | Fix effort |
|---------|---------|----------|-----------|
| F1: unchecked depth-index into R[]/L[] path vectors | `cons_loose.cpp`, `cons_majorityplus.cpp` | MAJOR | 6-line guard per call site (6 sites) |
| F2: leaf goodLabel invariant undocumented; `label[-1]` latent | `cons_majorityplus.cpp:62` | MINOR | 2-line loop reorder + assert |
| F3: `anc[-1]` dereference for non-kept root in `contract()` | `cons_loose.cpp:70` | MINOR | 1-line unconditional `keep[root]=1` |
| F4: `lcaDepth` walks past root `-1` on disconnected input | `rstar.cpp:79`, `local_consensus.cpp:48` | MINOR | 1-line root-check per loop |
| F5: O(n^2) queue in `gtp_no_common` | `bhv.cpp:167,197` | INFO | Replace `vector` with `deque` |
| F6: `pool.data` pointer stability | `Quartet.cpp:212` | INFO | Add assertion |

---

## Verdict

**CONDITIONAL SHIP**

F1 is the only finding with a realistic trigger from R-side misuse: an edge
matrix that is not strictly preorder (which ape's `Preorder()` guarantees but
which could be bypassed by a user calling the C++ entry point directly, or by an
edge matrix produced by a non-standard source) would produce FACT node metadata
that violates the depth-index invariant, causing an out-of-bounds vector read in
`looseMerge`, `updateCounter`, or `majorityPlusMerge`.  Add the six bounds-check
guards (F1) before a general release.

F2, F3, and F4 are latent risks that do not fire under any currently constructable
R input but are correctness-by-luck patterns.  Add the three assertions/guards
as a one-time hardening pass.

F5 and F6 are informational; no action required before ship.
