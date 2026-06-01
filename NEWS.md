# ConsTree 1.0.0 (2026-06-01)

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
