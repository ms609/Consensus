# Review: claude/cranky-kare-bcf5a4 @ 6329d96 (base 534b349)

## Verdict
**Block** — two correctness bugs on legal inputs: an infinite-loop hang on a
zero-length incompatible interior edge (critical), and a wrong/too-large distance
on trees with internal singleton (degree-2) nodes (major). Both live in input
regions with no test coverage. The core GTP geodesic on generic positive-length
binary/multifurcating trees is correct and well validated.

## F1. Zero-length incompatible interior edge -> infinite loop — CRITICAL
Where: src/bhv.cpp:68-75 (min_weight_vc normalisation) + :161-195 (gtp_no_common).
cap[S->a_i] = wA[i]/sA with sA = sum of squared A lengths. A GTP group whose
A-side (or B-side) lengths are all 0 gives sA=0 -> 0/0 = NaN capacities -> BFS
finds no augmenting path, cover degenerate, split never shrinks the subproblem ->
queue never empties. RAN repro-01 under `timeout 20` -> killed (EXIT 143).
Reachable by any tree with a length-0 interior branch incompatible with the other
tree (collapsed/soft polytomies are routinely encoded length-0). Fix: in
min_weight_vc guard sA==0/sB==0 (treat zero-weight vertices as free-to-cover, or
drop <TOL interior edges in from_r as tree_at already does on output).

## F2. Internal singleton node -> duplicated split -> wrong distance — MAJOR
Where: R/BHV.R:25-52 (.TreeToBHV) and :59-91 (.BHVToTree match()).
.TreeToBHV suppresses only the degree-2 *root* (ape::unroot), not internal
degree-2 nodes. Such a node has two edges inducing the SAME bipartition, so the
tree gets a duplicated split coordinate and sits at the wrong BHV point.
d(B, collapse.singles(B)) = 2.828427 for two mathematically identical trees
(repro-02, RAN). Also corrupts tree_at endpoints: d(.BHVTreeAt(A,B,1), B)=1.
Reachable via ape::drop.tip(..., collapse.singles = FALSE). Fix: collapse.singles
(or TreeTools equivalent) at the top of .TreeToBHV, next to the unroot call.

## F3. No edge.length -> opaque error — MINOR
Where: R/BHV.R:43. tree$edge.length NULL -> match()/assignment yields
"replacement has length zero". Validate and stop() with a clear message.

## Coverage notes
Exercised: oracle (15√2), midpoint, cone (incompatible), shared-orthant
(compatible), same-topology Euclidean, metric invariants, mean midpoint/
single-orthant/minimisation, variance avg-vs-sum, pairwise. These pin the vertex
cover and ratio sequence well on generic positive-length trees.
NOT exercised (and exactly where F1/F2 live): zero-length interior edges,
internal singleton nodes, polytomies inside a non-trivial geodesic.
No tautological tests found; triangle-inequality test is weak (1e-9 slack) but
not tautological.

## Verified fine
- Oracle 15√2 + midpoint {3,4}:2.5,{2,3,4,5}:2.5 exact.
- C++ == R prototype to 3.6e-15 over 400 random pairs (4-9 tips).
- Geodesic additivity d(A,M)+d(M,B)=d(A,B) to 1e-15 on positive-length trees
  (independent of prototype).
- Canonical (tip-0-excluding) compatibility == true 4-way bipartition
  compatibility: 0/72124 mismatches up to 8 tips.
- unroot sums basal edges (3+5 -> 8); rooted input handled.
- split<->length marshalling and round-trip exact (incl. tied lengths).
- Mean is a true local minimum (0/60 perturbations beat it), reproducible under
  same seed; stepDist = lambda*geo_dist(g) == d(old,new) exactly.
- Common-edge decomposition (nested + disjoint groups) additive to 1e-15.
- Star / polytomy / 3-tip / 4-tip-NNI / zero-length-compatible cases correct.
- n=1 pairwise and single-tree mean/variance handled.
- DESCRIPTION/NAMESPACE/RcppExports consistent; package builds; tests pass.
