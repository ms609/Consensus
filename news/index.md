# Changelog

## ConsTree 0.0.0.9007 (development)

**Performance overhaul (in progress).** Reimplementing each consensus
method with its fastest available algorithm, then profiling and
optimising (harness in `dev/profiling/`). Landed so far:

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
