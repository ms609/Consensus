# Red-team review: fact_tree.cpp / fact_tree.h
**Date:** 2026-06-03
**Verdict:** SHIP

Static review of the port against the FACT original (`dev/oracle/fact-src/tree.{h,cpp}`)
**plus** an empirical run: BRANCH source installed clean into `.agent-cons`
(`fact_tree.o` compiled this build, version 0.0.0.9007); all four consuming
methods (`Loose`/`Greedy`/`MajorityPlus`/`Frequency`) driven through the
discriminating tree shapes named in the brief, and the shipped
`tests/testthat/test-selection.R` lattice+idempotence suite re-run green (28/28).

`fact_tree.cpp/.h` is a faithful, correctly de-globalised port of FACT's
`tree::precompute()` and `tree::printNex()` plus a new `buildTreeFromEdge()`
adapter. The precompute body is character-identical to FACT modulo two
deliberate, behaviour-neutral changes (chained-assignment split; dropped
`assert`). No correctness defect found, in code or in execution. Nothing here
blocks any of the four methods.

---

## No findings rising to CRITICAL or MAJOR.

The items below are MINOR / defensive notes only.

## F1 — dropped `assert(num == N)` removes the sole malformed-tree guard (MINOR)
`Tree::precompute()` (fact_tree.cpp:79) drops FACT's terminal
`assert(num == N)` with the comment "on a well-formed tree this always holds."
True — but it was FACT's *only* structural sanity check. The C++ entry points now
trust the R caller completely: that `.FactEdges()` delivers a contiguous
preorder edge matrix with tips 1..nTip and ≥4 leaves. That contract is in fact
enforced upstream (`.PrepareTrees()` guards `< 4` leaves; `Preorder()` +
`RenumberTips()` + `RootTree()` guarantee the shape), so this is not reachable
from a legitimate R call. Note, not fix: if a future caller bypasses
`.FactEdges()` and feeds a malformed edge matrix, `precompute()` will silently
mis-relabel rather than abort. Acceptable for an internal primitive; flagged so
the invisible precondition is on record.

## F2 — `buildTreeFromEdge` root fallback `nRow > 0 ? … : 0` is dead code (MINOR/defensive)
fact_tree.cpp:130. An empty edge matrix (`nRow == 0`) cannot occur: the n=3
short-circuit and the `< 4` leaf guard in `.PrepareTrees()` mean the C++ is only
ever entered with ≥4 tips, hence a non-empty preorder matrix whose first row's
parent is the root. The ternary's `: 0` arm is unreachable defensive code — same
class as the `# nocov` post-condition guards elsewhere in the package. Fine to
keep; will never be covered.

## Verified fine

- **Construction from the preorder edge matrix (`buildTreeFromEdge`).** The
  ape→FACT index shift is correct: ape node `v` → FACT node `v-1`; tips
  `1..nTip` carry taxon `v` via `leaf[v-1]=v`, `idx[v]=v-1`. `nNode` is taken as
  the max id over both edge columns, so internal nodes above `nTip` are sized in.
  Root = `edge(0,0)-1`, valid because a preorder matrix lists the root's edge
  first; `parent[root]` reset to −1. **Empirically:** star (single internal node,
  n children), caterpillar, fully balanced, random rooted, and mixed-shape inputs
  all reconstruct and round-trip to valid Newick under all four methods.

- **Day's leaf relabelling (`precompute`).** Body matches FACT line-for-line. The
  in-order DFS numbering assigns each subtree a contiguous `[minL,maxL]`; the
  invariant "taxon with label i lies in the subtree hashed at i" holds because
  `label`, `minL`, `maxL` are set together at first leaf visit and merged on
  pop. The one source edit — FACT's `size[t]=labelled[t]=1` chained assignment
  split into two statements — is value-identical (both LHS get 1). `idx` is
  rederived from `leaf` at entry, so a relabelled/contracted tree is hashed
  consistently. **Empirically:** idempotence (consensus of identical trees ==
  that tree, compared by label-based canonical split set) holds for every shape.

- **Array bounds.** Per-taxon arrays sized `N+5`, per-node `cnt+5`, matching
  FACT's guard padding. The Day's-hash arrays `minH/maxH/H` are sized `N+5` and
  indexed by *label* (1..N) — correct, since labels never exceed N. `H[label]`
  stores a *node id* (0..cnt-1) and is read back as an index into `cnt+5`-sized
  arrays (`can[A.H[i]]` in cons_loose, `goodLabel`) — consistent. Every entry
  read after `precompute()` is one it wrote (H initialised to −1 and gated on
  `H[i] >= 0` at the read sites), so no stale-value hazard.

- **Integer types.** All indices are `int`. Node count is ~2N (one internal per
  bifurcation); at the documented n≈10 000 ceiling, ids and the
  `POS[b]-POS[a]+1` / `DEPTH` arithmetic in the consumers stay far inside int32.
  No multiplication of two index-scale quantities anywhere in the primitive.
  **Empirically:** n=500 × 5 trees runs clean for all four methods (the suite's
  larger cross-checks go higher).

- **Newick serialisation (`newick`).** Iterative, matches FACT `printNex` exactly
  except it returns a `std::string` with no trailing `;` (the R side appends it).
  Produces balanced parens for star, caterpillar, balanced and random shapes.
  **Single-child root** emits `(X)` — faithful to FACT and parsed without error
  by `ape::read.tree`; not malformed.

- **Reuse / reinitialisation.** This is the whole point of the port and it holds.
  No file-scope globals (FACT's `extern tree *T; extern int numTaxas` are gone);
  state is per-`Tree`, passed by argument. `precompute()` re-zeroes `size`,
  `minL`, `maxL`, `H`, `label` and `labelled` on every call, so a `Tree` reused
  across merges carries no stale hash. Trees are `std::vector`-backed, hence
  copyable/movable/self-freeing; the loose algorithm's by-value `Tree A, Tree B`
  copies are intended (callers' originals are preserved) and were verified by the
  repeated-merge path producing correct output.

## Coverage gaps

- **No direct unit test of the primitive.** `precompute()`, `newick()` and
  `buildTreeFromEdge()` are exercised only transitively through the four methods'
  R entry points and the oracle. There is no C++-level test that, e.g., asserts
  the `[minL,maxL]` contiguity invariant or the H/idx round-trip directly. Low
  risk (the transitive coverage is strong and the oracle is exact for
  Loose/MajorityPlus), but a malformed edge matrix or a precompute regression
  would only surface as a wrong consensus, not a localised failure.

- **`newick()` single-child-root path is faithful but untested in isolation.**
  The `(X)` output is reachable only if a consumer hands `newick()` a root with
  one child; in the current pipeline the consumers never do, so this arm is
  validated by reading FACT, not by a fixture. If a future method can produce a
  one-child root, add a Newick round-trip fixture.

- **`precompute()` is dead weight for Greedy.** Per the header, Greedy does not
  call `precompute()`; it is exercised only by Loose / MajorityPlus / Frequency
  (and the oracle for the first two). Scope any "precompute is verified" claim to
  those three — Greedy's correctness says nothing about the Day's-hash path.

- **F1 / F2 paths (malformed input, empty edge matrix) are by construction
  unreachable**, so intentionally uncovered. Noted above rather than tested.
