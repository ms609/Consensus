# Changelog

## ConsTree 0.0.0.9008 (development)

**Performance overhaul (in progress).** Reimplementing each consensus
method with its fastest available algorithm, then profiling and
optimising (harness in `dev/profiling/`). Landed so far:

- [`Adams()`](https://constree.github.io/reference/Adams.md) now uses a
  C++ implementation of the asymptotically efficient O(kn log n)
  centroid-path algorithm of Jansson, Li & Sung (2017) (their FACT
  toolkit, used with permission), replacing the previous pure-R
  recursion and per-level Newick round-trip. The shared “spine” of the
  consensus (the leaves lying under the heavy child of every input
  tree’s root) is expanded iteratively down the trees’ centroid paths in
  unison, recursing only on the smaller off-spine blocks. The Adams
  consensus tree is unique, so the output is identical to the classical
  definition (validated clade-for-clade against the slow FACT
  reference); it is orders of magnitude faster (e.g. 50 trees on 200
  leaves: ~0.9 s to ~0.02 s).
- [`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md)
  now uses a C++ port of the optimal O(kn) `majorityPlusConsensus`
  algorithm of Jansson, Shen & Sung (2016) (their FACT toolkit, used
  with permission), replacing the previous R pairwise compatibility
  matrix. It is typically two to three orders of magnitude faster and
  clears inputs that previously timed out (e.g. 50 trees on 100–200
  leaves). Majority-rule (+) is a deterministic count rule (a clade is
  kept when displayed by strictly more trees than contradict it, with no
  frequency tie-break), so the output is the FACT reference exactly,
  verified even at n \> 60.
- [`Loose()`](https://constree.github.io/reference/Loose.md) now uses a
  C++ port of the asymptotically efficient `looseConsensusFast`
  algorithm of Jansson, Shen & Sung (2016) (their FACT toolkit, used
  with permission), replacing the previous R pairwise compatibility
  matrix. The input trees are merged into a one-way compatible tree by
  repeated linear-time consecutive-range queries, then the clusters
  compatible with every input are retained; the output is validated to
  match the FACT reference exactly (the loose consensus is unique, so
  there is no tie-break ambiguity).
- [`Greedy()`](https://constree.github.io/reference/Greedy.md) now uses
  a C++ port of the asymptotically efficient `greedyConsensusFast`
  algorithm of Jansson, Shen & Sung (2016) (their FACT toolkit, used
  with permission), replacing the previous R pairwise compatibility
  matrix. It is typically two to three orders of magnitude faster and
  lifts the practical size ceiling (e.g. 50 trees on 200 leaves: ~15 s
  to ~10 ms); the output is validated to match the FACT reference.
  Equally frequent incompatible splits may now be resolved in a
  different but equally valid order.
- Added a reentrant, allocation-safe C++ tree primitive
  (`src/fact_tree.*`) shared by the fast split-selection methods.
- [`RStar()`](https://constree.github.io/reference/RStar.md) no longer
  caps at 200 leaves. The dense `O(n^3)` triplet tensor is replaced by
  per-tree constant-time LCA queries (`O(kn^2)` memory), and the
  strong-cluster assembly is tightened from `O(n^4)` to about `O(n^3)`
  via the single-linkage (Apresjan) construction of Jansson, Sung, Vu &
  Yiu (2016). The R\* tree is unchanged — verified clade-for-clade
  against the previous implementation up to 200 leaves — and large
  inputs are now bounded by running time rather than a memory wall.
  [`RStar()`](https://constree.github.io/reference/RStar.md) now also
  rejects duplicated tip labels and trees on differing leaf sets (which
  would otherwise be scored silently against mismatched taxa), and the
  internal LCA self-check now covers every input tree rather than only
  the first.
- [`Frequency()`](https://constree.github.io/reference/Frequency.md) now
  uses a C++ port of the near-linear *O*(*kn* log *n*)
  frequency-difference algorithm of Jansson, Sung, Tabatabaee & Yang
  (2024, STACS; their FDCT reference implementation, used with
  permission), replacing the previous R *O*(\_s_²) per-split frequency
  comparison. The old approach exceeded a one-minute budget beyond ~100
  leaves on 50 trees, where the port returns in well under a second. The
  port is boost-free (the upstream `dynamic_bitset` is dead code on this
  near-linear path) and its output is validated to match the FDCT
  `freqdiff` reference exactly.
- [`Frequency()`](https://constree.github.io/reference/Frequency.md) is
  now robust to very deep trees. The filter’s centroid-path
  decomposition descends each subtree’s `children[0]` assuming it is the
  heaviest child, but the deep copy and merge re-add children in id /
  `m` order and silently undid the upstream
  [`reorder()`](https://rdrr.io/r/stats/reorder.factor.html); on a deep
  or caterpillar input the path degenerated to a single leaf, making the
  filter *O*(\_n_²) in both time and heap and throwing a catchable
  `std::bad_alloc` around 30 000 leaves. The heavy-child-first invariant
  is now re-established
  ([`reorder()`](https://rdrr.io/r/stats/reorder.factor.html) +
  `fix_tree()`) at the top of
  [`filter()`](https://rdrr.io/r/stats/filter.html), restoring the
  paper’s *O*(*kn* log *n*): an opposite-caterpillar pair on 30 000
  leaves goes from `bad_alloc` to ~0.5 s (100 000 leaves in ~2 s), and
  incongruent random ensembles that previously timed out now finish in
  milliseconds (e.g. 50 trees on 50 leaves: 19 s to 0.03 s). The
  decomposition is a performance device only, so the unique
  frequency-difference split set is unchanged (the change reorders the
  Newick serialisation but not the consensus; still exact against the
  FDCT `freqdiff` reference).

First public release: a consensus-tree toolkit built on
[TreeTools](https://ms609.github.io/TreeTools/).

Split-based consensus methods:

- [`Strict()`](https://constree.github.io/reference/Strict.md),
  [`Majority()`](https://constree.github.io/reference/Majority.md) /
  [`MajorityRule()`](https://constree.github.io/reference/Majority.md):
  wrappers around
  [`TreeTools::Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.html).
- [`Loose()`](https://constree.github.io/reference/Loose.md): loose
  (semi-strict / combinable-component) consensus.
- [`Greedy()`](https://constree.github.io/reference/Greedy.md): greedy
  (extended majority-rule) consensus.
- [`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md):
  majority-rule (+) consensus.
- [`Frequency()`](https://constree.github.io/reference/Frequency.md):
  frequency-difference consensus.

Rooted-tree consensus methods:

- [`Adams()`](https://constree.github.io/reference/Adams.md): Adams
  consensus.
- [`RStar()`](https://constree.github.io/reference/RStar.md): R\*
  consensus, assembled from the strong clusters of the majority
  resolved-triplet set (Jansson, Sung, Vu & Yiu 2016).
- [`Local()`](https://constree.github.io/reference/Local.md): local
  consensus (MinRLC / MinILC).

Distance- and branch-length-based summaries:

- [`Average()`](https://constree.github.io/reference/Average.md):
  distance-based average consensus.
- [`Quartet()`](https://constree.github.io/reference/Quartet.md): tree
  minimizing the summed symmetric quartet distance to the inputs
  (Takazawa et al. 2026).
- Billera-Holmes-Vogtmann treespace summaries:
  [`BHVDistance()`](https://constree.github.io/reference/BHVDistance.md)
  (Owen-Provan geodesic distance), `BHVPairwiseDistances()`,
  [`BHVMean()`](https://constree.github.io/reference/BHVMean.md)
  (Fréchet mean) and
  [`BHVVariance()`](https://constree.github.io/reference/BHVMean.md).

An introductory vignette (`ConsTree`) demonstrates each method family on
worked exemplars.
