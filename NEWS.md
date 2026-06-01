# Consensus 0.0.0.9000 (development)

- Initial package scaffold, building on
  [TreeTools](https://ms609.github.io/TreeTools/).
- Add `Strict()` and `Majority()` / `MajorityRule()`: thin, consistently named
  wrappers around `TreeTools::Consensus()`.
- Add `QuartetConsensus()`: a consensus tree minimizing the summed symmetric
  quartet distance to the input trees, found with a greedy add-and-prune
  heuristic (Takazawa et al. 2026).  The C++ core is self-contained.
- Add Billera-Holmes-Vogtmann (BHV) treespace summaries with branch lengths:
  `BHVDistance()` (Owen-Provan geodesic distance), `BHVPairwiseDistances()`,
  `BHVMean()` (iterative Fréchet mean) and `BHVVariance()`.  The geodesic core
  is implemented in C++.
