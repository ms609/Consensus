# Changelog

## ConsTree 1.0.0 (2026-06-01)

First public release: a consensus-tree toolkit built on
[TreeTools](https://ms609.github.io/TreeTools/).

Split-based consensus methods:

- [`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md),
  [`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md)
  /
  [`MajorityRule()`](https://ms609.github.io/ConsTree/reference/Majority.md):
  wrappers around
  [`TreeTools::Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.html).
- [`Loose()`](https://ms609.github.io/ConsTree/reference/Loose.md):
  loose (semi-strict / combinable-component) consensus.
- [`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md):
  greedy (extended majority-rule) consensus.
- [`MajorityPlus()`](https://ms609.github.io/ConsTree/reference/MajorityPlus.md):
  majority-rule (+) consensus.
- [`Frequency()`](https://ms609.github.io/ConsTree/reference/Frequency.md):
  frequency-difference consensus.

Rooted-tree consensus methods:

- [`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md):
  Adams consensus.
- [`RStar()`](https://ms609.github.io/ConsTree/reference/RStar.md): R\*
  consensus, assembled from the strong clusters of the majority
  resolved-triplet set (Jansson, Sung, Vu & Yiu 2016).
- [`Local()`](https://ms609.github.io/ConsTree/reference/Local.md):
  local consensus (MinRLC / MinILC).

Distance- and branch-length-based summaries:

- [`Average()`](https://ms609.github.io/ConsTree/reference/Average.md):
  distance-based average consensus.
- [`Quartet()`](https://ms609.github.io/ConsTree/reference/Quartet.md):
  tree minimizing the summed symmetric quartet distance to the inputs
  (Takazawa et al. 2026).
- Billera-Holmes-Vogtmann treespace summaries:
  [`BHVDistance()`](https://ms609.github.io/ConsTree/reference/BHVDistance.md)
  (Owen-Provan geodesic distance), `BHVPairwiseDistances()`,
  [`BHVMean()`](https://ms609.github.io/ConsTree/reference/BHVMean.md)
  (Fréchet mean) and
  [`BHVVariance()`](https://ms609.github.io/ConsTree/reference/BHVMean.md).

An introductory vignette (`ConsTree`) demonstrates each method family on
worked exemplars.
