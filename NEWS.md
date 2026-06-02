# ConsTree 0.0.0.9007 (development)

**Performance overhaul (in progress).** Reimplementing each consensus method
with its fastest available algorithm, then profiling and optimising (harness in
`dev/profiling/`). Landed so far:

- `MajorityPlus()` now uses a C++ port of the optimal O(kn)
  `majorityPlusConsensus` algorithm of Jansson, Shen & Sung (2016) (their FACT
  toolkit, used with permission), replacing the previous R pairwise
  compatibility matrix. It is typically two to three orders of magnitude faster
  and clears inputs that previously timed out (e.g. 50 trees on 100–200 leaves).
  Majority-rule (+) is a deterministic count rule (a clade is kept when displayed
  by strictly more trees than contradict it, with no frequency tie-break), so the
  output is the FACT reference exactly, verified even at n > 60.
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
- `RStar()` no longer caps at 200 leaves. The dense `O(n^3)` triplet tensor is
  replaced by per-tree constant-time LCA queries (`O(kn^2)` memory), and the
  strong-cluster assembly is tightened from `O(n^4)` to about `O(n^3)` via the
  single-linkage (Apresjan) construction of Jansson, Sung, Vu & Yiu (2016). The
  R\* tree is unchanged — verified clade-for-clade against the previous
  implementation up to 200 leaves — and large inputs are now bounded by running
  time rather than a memory wall.
- `Frequency()` now uses a C++ port of the near-linear _O_(_kn_ log _n_)
  frequency-difference algorithm of Jansson, Sung, Tabatabaee & Yang (2024,
  STACS; their FDCT reference implementation, used with permission), replacing
  the previous R _O_(_s_²) per-split frequency comparison. The old approach
  exceeded a one-minute budget beyond ~100 leaves on 50 trees, where the port
  returns in well under a second. The port is boost-free (the upstream
  `dynamic_bitset` is dead code on this near-linear path) and its output is
  validated to match the FDCT `freqdiff` reference exactly.

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
