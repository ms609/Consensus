# ConsTree 0.0.0.9008 (development)

**Performance overhaul (in progress).** Reimplementing each consensus method
with its fastest available algorithm, then profiling and optimising (harness in
`dev/profiling/`). Landed so far:

- `Adams()` now uses a C++ implementation of the asymptotically efficient
  O(kn log n) centroid-path algorithm of Jansson, Li & Sung (2017) (their
  FACT toolkit, used with permission), replacing the previous pure-R recursion
  and per-level Newick round-trip. The shared "spine" of the consensus (the
  leaves lying under the heavy child of every input tree's root) is expanded
  iteratively down the trees' centroid paths in unison, recursing only on the
  smaller off-spine blocks. The Adams consensus tree is unique, so the output is
  identical to the classical definition (validated clade-for-clade against the
  slow FACT reference); it is orders of magnitude faster (e.g. 50 trees on 200
  leaves: ~0.9 s to ~0.02 s).
- `Loose()` now uses a C++ port of the asymptotically efficient
  `looseConsensusFast` algorithm of Jansson, Shen & Sung (2016) (their FACT
  toolkit, used with permission), replacing the previous R pairwise
  compatibility matrix. The input trees are merged into a one-way compatible
  tree by repeated linear-time consecutive-range queries, then the clusters
  compatible with every input are retained; the output is validated to match
  the FACT reference exactly (the loose consensus is unique, so there is no
  tie-break ambiguity).
- `Greedy()` now uses a C++ port of the asymptotically efficient
  `greedyConsensusFast` algorithm of Jansson, Shen & Sung (2016) (their FACT
  toolkit, used with permission), replacing the previous R pairwise
  compatibility matrix. It is typically two to three orders of magnitude faster
  and lifts the practical size ceiling (e.g. 50 trees on 200 leaves: ~15 s to
  ~10 ms); the output is validated to match the FACT reference. Equally frequent
  incompatible splits may now be resolved in a different but equally valid order.
- Added a reentrant, allocation-safe C++ tree primitive (`src/fact_tree.*`)
  shared by the fast split-selection methods.

First public release: a consensus-tree toolkit built on
[TreeTools](https://ms609.github.io/TreeTools/).

Split-based consensus methods:

- `Strict()`, `Majority()` / `MajorityRule()`: wrappers around
  `TreeTools::Consensus()`.
- `Loose()`: loose (semi-strict / combinable-component) consensus.
- `Greedy()`: greedy (extended majority-rule) consensus.
- `MajorityPlus()`: majority-rule (+) consensus.
- `Frequency()`: frequency-difference consensus.

Rooted-tree consensus methods:

- `Adams()`: Adams consensus.
- `RStar()`: R* consensus, assembled from the strong clusters of the majority
  resolved-triplet set (Jansson, Sung, Vu & Yiu 2016).
- `Local()`: local consensus (MinRLC / MinILC).

Distance- and branch-length-based summaries:

- `Average()`: distance-based average consensus.
- `Quartet()`: tree minimizing the summed symmetric quartet distance to the
  inputs (Takazawa et al. 2026).
- Billera-Holmes-Vogtmann treespace summaries: `BHVDistance()` (Owen-Provan
  geodesic distance), `BHVPairwiseDistances()`, `BHVMean()` (Fréchet mean) and
  `BHVVariance()`.

An introductory vignette (`ConsTree`) demonstrates each method family on worked
exemplars.
