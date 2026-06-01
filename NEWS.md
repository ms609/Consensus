# Consensus 0.0.0.9000 (development)

- Initial package scaffold, building on
  [TreeTools](https://ms609.github.io/TreeTools/).
- Add `Strict()` and `Majority()` / `MajorityRule()`: thin, consistently named
  wrappers around `TreeTools::Consensus()`.
- Add `Quartet()`: an information-maximizing quartet consensus tree,
  found with a greedy add-and-prune heuristic that can also drop rogue taxa
  via `neverDrop` and tune the misinformation `penalty` (Takazawa et al. 2026).
  The C++ core is self-contained.  (Experimental: rogue-dropping is a
  work in progress; its tests are currently expected to fail.)
- Add Billera-Holmes-Vogtmann (BHV) treespace summaries with branch lengths:
  `BHVDistance()` (Owen-Provan geodesic distance), `BHVPairwiseDistances()`,
  `BHVMean()` (iterative Fréchet mean) and `BHVVariance()`.  The geodesic core
  is implemented in C++.
