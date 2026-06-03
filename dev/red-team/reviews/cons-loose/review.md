# Red-team review: cons_loose.cpp
**Date:** 2026-06-03
**Verdict:** SHIP

Port is structurally faithful to FACT `looseConsensusFast` / `looseMerge` /
`contract`.  Every concern raised in the slot-1 review (feature-loose-fast) has
been re-examined at source level and against live oracle runs.  No correctness
blockers found.

---

## Scope

`src/cons_loose.cpp` (the full file) + `R/selection.R` (`Loose()` wrapper,
`.FactEdges()`, `.PrepareTrees()`, `.RootLikeFirst()`).  The FACT reference
source (`dev/oracle/fact-src/loose.cpp`, `strict.cpp`) was read line-by-line
for the faithfulness diff.  The oracle (`dev/oracle/check-oracle.R`) was run
fresh this session; all comparisons matched.

---

## 1. looseMerge op==1 — BEFORE/AFTER/pid path

### Are BEFORE/AFTER pointers updated atomically?

`BEFORE[a][i]` and `AFTER[a][i]` are incremented at FACT loose.cpp:216-217
and cons_loose.cpp:229-230 only when `B.parent[a] == B.parent[b]`, i.e. when
the two bounding nodes of the A-cluster's leaf range are siblings in B.
Both increments happen in the same basic block (no control flow between them),
and there is no aliasing risk (`a = B.parent[a]` is the same node for both).
Atomic in the only sense that matters here: both or neither.

### Is `pid` set correctly for the newly inserted node?

`pid[v]` is the child-index of node `v` within its parent in B, assigned at
cons_loose.cpp:176 (`pid[B.G[t1.first][c]] = cur++`) and FACT loose.cpp:172-173.
In the op==1 insertion loop (cons_loose.cpp:258-271), the index `i` iterates
over `B.G[a]`'s children and `BEFORE[a][i]` / `AFTER[a][i]` are keyed to that
same child-index.  The roles of pid in the query phase and `i` in the insertion
phase are consistent: both count children in B's (reordered) adjacency list.

There is one subtlety: pid is assigned AFTER the `cur = 0` reset at line 175,
inside the same postorder DFS pass that fills BEFORE/AFTER.  Because the
assignment loop uses `cur` starting at 0 and increments for each child, and the
insertion loop uses `i` starting at 0 and increments for each child, the two
numbering schemes are identical and the two loops are synchronized — this is
verified by inspection and matches FACT exactly.

### Does insertion handle multiple A-children among B's children?

When multiple A-clusters are each compatible with B and each maps to a pair of
B-children `(a_j, b_j)` that share the same parent, multiple BEFORE/AFTER
increments accumulate at that parent's child slots.  The insertion loop at
cons_loose.cpp:255-271 handles this by iterating `a` over ALL of B's nodes and
for each, iterating `i` over ALL children; `cur` tracks the most recently
active intermediate node.  The AFTER-step `cur = ret.parent[cur]` climbs back
to the correct ancestor before the next child is spliced in.  This is a faithful
translation of FACT's inner loop (loose.cpp:239-256): the same `cnt`/`cur`
semantics, the same BEFORE-descend / AFTER-ascend traversal.

The post-condition guard (`next != ret.cnt` -> Rcpp::stop) detects any mismatch
in the bookkeeping at O(1) cost.  It is the assert at FACT loose.cpp:258 promoted
from a compile-out `assert` to a loud runtime check.  (Prior review F5 called
this "dropped" — it was not dropped, only upgraded.)

**Empirical check (new this review):** a three-tree input with two compatible
but differently-resolved flat trees forces `newNodes > 0` in the first
`looseMerge` call (op==1):

    T1: ((t1,t2),(t3,t4),t5,t6)
    T2: ((t3,t4),(t5,t6),t1,t2)
    T3: ((t1,t2),(t5,t6),t3,t4)

`mine` = 3 splits, `fact` = 3 splits, match = TRUE (oracle run, this session).
This is the first time the BEFORE/AFTER insertion path has been verified against
the oracle rather than by static reading alone (closing prior F4).

---

## 2. Contract step

The `contract()` function in `cons_loose.cpp` is a verbatim port of
`strict.cpp:53`.  Differences from the original are purely cosmetic (VLA ->
`std::vector`, `memset` -> `std::fill`/constructor zero-init).

### (a) Is the compatibility check symmetric?

In the two-phase loose algorithm, op==0 (`looseMerge(ret, T[i], 0, ...)`) marks
`ret.good[i] &= can[A.H[i]]` for each non-trivial A-cluster.  This marks A's
clusters as good iff compatible with B.  This is NOT symmetric (it answers "is
A's cluster compatible with B?", not "is B compatible with A?"), and it should
not be: the algorithm keeps only clusters of the one-way-compatible tree that are
compatible with ALL input trees, which is the correct definition of loose
consensus (Jansson, Shen & Sung 2016 Theorem 1).  The asymmetry is intentional
and correct.

### (b) Does contraction create dangling references?

`contract()` rebuilds the tree from scratch (a new `Tree ret` of the correct
size), rewiring each kept node to its nearest kept ancestor via a DFS that tracks
`anc[]`.  The DFS assignment `anc[A.G[t1.first][t1.second]] = anc[t1.first]`
propagates the nearest kept ancestor BEFORE the child is pushed to the stack —
matching FACT strict.cpp:89.  The final `parent` array is rebuilt from scratch
at cons_loose.cpp:79-80.  No dangling references: removed nodes are simply
absent from `ret`, and `ret.G` is only populated for kept nodes.

### (c) Does contraction correctly handle the root?

`ret.root = label[A.root]` (cons_loose.cpp:59).  The root is always kept
(every leaf is kept; the root is kept because `keep[A.idx[i]] = 1` covers all
leaves, and `tmp` starts at `A.N`).  The DFS guard `if (t1.first != A.root &&
keep[t1.first])` skips adding the root to any parent's child list (since the
root has no parent), which is correct.  The star case (all non-trivial clusters
contracted) leaves `ret` with `root` and all leaves directly as children: 0
internal non-root nodes, 0 splits — verified by the star-output oracle check
this session (6 tips, 0 splits, MATCH).

---

## 3. DEPTH array and consecutive-range query

### Initialisation

`DEPTH[B.root] = 0` is set at cons_loose.cpp:145, and
`DEPTH[B.G[t1.first][t1.second]] = DEPTH[t1.first] + 1` propagates it in the
DFS at cons_loose.cpp:160.  This matches FACT loose.cpp:143,154.  The initial
fill of `DEPTH` is `std::vector<int> DEPTH(B.cnt + 5, 0)` at cons_loose.cpp:97 —
all entries zero.  After the DFS the root has DEPTH=0 and every other node has
been set via its parent's DEPTH+1, so only the root carries the initial zero.
This is correct because `DEPTH[B.root] = 0` is explicitly set before the DFS,
and no child of root is processed until root is pushed.

### Half-open vs. closed intervals

`POS[numTaxas+1] = numTaxas+1` and `POS[numTaxas+2] = -1` are sentinel values
for the `B.minL`/`B.maxL` of internal nodes that haven't yet acquired real leaf
bounds.  The query `POS[b] - POS[a] + 1 != A.size[t1.first]` (cons_loose.cpp:199)
uses a CLOSED interval `[a, b]` of the sequential leaf labels.  This matches
FACT loose.cpp:191.  The sentinel values prevent non-leaf nodes from passing
this check spuriously: `POS[numTaxas+2] - POS[numTaxas+1] + 1 = -1 - (numTaxas+1) + 1 < 0`,
which is never equal to size >= 1.

The L/R path construction (cons_loose.cpp:164-165) and the depth-indexed lookup
(cons_loose.cpp:209-218) match FACT verbatim.  The vector indexing
`R[b][DEPTH[B.idx[b]] - DEPTH[a]]` and `L[a][DEPTH[B.idx[a]] - DEPTH[b]]`
are safe because the depth difference is bounded by the number of nodes in the
path, which equals the vector size: `L[x]` / `R[x]` are exactly the left/right
path from the innermost node to the root, so the index cannot exceed the vector
length.

---

## 4. Star-tree output

When no A-cluster survives the op==0 pass, `ret.good[i] = 0` for every i, so
`contract()` keeps only the root and the leaves.  The Newick serialiser
(`newick()` in `fact_tree.cpp`) emits `(c1,c2,...,cN)` for the root: one open
paren, leaf labels comma-separated, one close paren — a valid star topology.
The R side appends `;` and calls `ape::read.tree()`.

**Empirical check (new this review):**

    T1: ((t1,t2),(t3,t4),(t5,t6))
    T2: ((t1,t3),(t2,t5),(t4,t6))

Every split in T1 is incompatible with a split in T2.  `Loose()` returned:
nSplits=0, NTip=6, oracle match=TRUE.  `ape::read.tree()` accepted the Newick
string without error.  The output tree has all 6 tips as direct children of the
root and no internal nodes besides the root.

---

## 5. Polytomous inputs

`buildTreeFromEdge` (fact_tree.cpp:110) builds the Tree from an ape preorder
edge matrix: it calls `t.G[p].push_back(c)` for every (parent, child) edge,
so a node with k children appears with `G[p].size() == k`.  No binary-only
assumption.

In `looseMerge`, the reordering DFS (`B.G[t1.first].clear()` then rebuild from
`L[]`) operates on `B.G` which may have k > 2 children — it clears and
rebuilds the child list, preserving the count.  The BEFORE/AFTER arrays are
sized `B.G[t1.first].size()` (cons_loose.cpp:168-169), one slot per child
regardless of degree.  The pid assignment iterates all children.  The insertion
loop (op==1) iterates all B-children by index.  No binary assumption anywhere
in these paths.

**Empirical check:** the oracle script (`check-oracle.R` lines 151-178) runs two
polytomy cases (trichotomies n=8 and nested polytomies n=9) under both rooted
flags.  All four cases MATCH this session.  The prior review noted this
empirically; this review confirms the mechanism is also correct by inspection.

---

## 6. Large n (n > 60) — no BUCKET_SIZE dependency

`cons_loose.cpp` has no `BUCKET_SIZE` reference.  `BUCKET_SIZE` appears only
in FACT's `looseConsensusSlow` (the O(s^2) bit-packing slow path), which this
port replaces with the structural consecutive-range algorithm.  The fast path is
purely structural: Day's label relabelling + DFS index arithmetic.  All storage
is `std::vector` with runtime sizes proportional to `numTaxas` or `B.cnt`.

**Empirical check:** oracle run at n=80 and n=137, both FACT-exact (check-oracle.R
lines 115-136), confirmed this session.

---

## Faithfulness diff — the six deliberate departures

The translation from FACT is faithful.  The only intentional departures from the
original are all safe:

1. **VLA -> `std::vector`** — identical semantics, portable to MSVC/R CMD check.
2. **Globals -> `numTaxas` parameter** — de-globalised; no semantic change.
3. **`goto end` -> `skip` bool** — the mutations `a = L[a].back()` etc. are in
   `else` arms that don't execute when `skip` is set; the parent min/max
   propagation runs unconditionally after, matching FACT's post-`end:` code.
   Logically equivalent.
4. **Shadowed inner `int t1` -> `q1`/`q2`** — FACT's inner loop reuses the name
   `t1` (loose.cpp:214) for a temporary; the port renames to `q1`/`q2` to avoid
   shadowing the outer loop variable.  Same values.
5. **`assert(cnt==ret.cnt)` -> `Rcpp::stop`** — promoted from a compile-out
   assertion to a loud runtime check (cons_loose.cpp:280-284).  Strictly
   stronger.  (Prior finding F5 miscalled this as "dropped".)
6. **Dead `memset(A.size,...)` dropped** — FACT loose.cpp:111 writes byte-short
   to `A.size` then immediately calls `A.precompute()` which re-zeros `size`
   entirely; the dead write is omitted, relying on `precompute()`.  Correct.

---

## Prior findings re-examined

| # | Prior finding | Status |
|---|---|---|
| F1 | test-loose.R:75-90 vacuous (star output) | Open; test hardening only |
| F2 | test-loose.R:8-38 loose==strict (no real insertion tested) | Open; test hardening only |
| F3 | `if(NSplits==0) next` latent vacuity | Open; cosmetic |
| F4 | op==1 insertion path untested in shipped suite | **Closed by this review** — oracle-verified with newNodes>0 case |
| F5 | Dropped assert | **Closed** — present as `Rcpp::stop` at cons_loose.cpp:280-284 |
| F6 | check-oracle.R print-only (no stop()) | **Closed** — `failCount` + `quit(status=1L)` present at lines 186-189 |

---

## Findings (this review)

Correctness findings: all PASS.  One test-hardening recommendation remains open:

- **T1 (minor, no code change required):** The op==1 newNodes>0 insertion path
  is exercised only by the oracle script, not by `tests/testthat/test-selection.R`.
  Promote one compatible-but-differently-resolved case (e.g. the three-tree T1/T2/T3
  example above) into `test-selection.R` with a pinned `splitSet` expectation
  checking that all three splits are recovered.  This closes the last gap in the
  shipped regression suite for the insertion path.  Not a blocker: the oracle
  script is now asserting, and the method is deterministic.

---

## Evidence summary

Oracle run this session (`dev/oracle/check-oracle.R`): 0 divergences across all
comparisons, including:
- Binary inputs, n=8/9/10, random and conflict seeds, rooted=0 and rooted=1
- n=80, n=137 structural large-n check
- Polytomous inputs: n=8 trichotomies and n=9 nested polytomies, rooted=0 and
  rooted=1
- Star-output case (6-tip fully-conflicting pair): 0 splits, MATCH (new)
- newNodes>0 op==1 case (three compatible flat trees): 3 splits, MATCH (new)
