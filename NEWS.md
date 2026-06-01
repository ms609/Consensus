# Consensus 0.0.0.9005 (development)

- `RStar()`: assemble the consensus from the **strong clusters** of the majority
  resolved-triplet set (Jansson, Sung, Vu & Yiu 2016, Lemma 1.1) — the exact R*
  definition.  This replaces the earlier Aho-BUILD assembly, which could over- or
  under-resolve on conflicting input (it did not guarantee `r(τ) ⊆ R_maj`).
  Validated by a brute-force strong-cluster oracle and by majority-rule
  refinement (every rooted majority clade now appears in the result).

# Consensus 0.0.0.9004 (development)

- `RStar()`: R* consensus of rooted trees, via compiled C++ (Rcpp).  Selects the
  uniquely-favoured (strict plurality) rooted triplet per taxon triple;
  polynomial, so not restricted to small leaf counts (capped at 200 leaves as a
  memory safeguard).

# Consensus 0.0.0.9003 (development)

- `Local()`: local consensus of rooted trees (MinRLC / MinILC) via compiled
  C++ (Rcpp); exact-exponential, limited to 20 leaves.

# Consensus 0.0.0.9002 (development)

- `Adams()`: Adams consensus of rooted trees (validated against the reference
  FACT implementation).

# Consensus 0.0.0.9001 (development)

- Initial package scaffold, building on
  [TreeTools](https://ms609.github.io/TreeTools/).
- `Strict()` and `Majority()` / `MajorityRule()`: thin, consistently named
  wrappers around `TreeTools::Consensus()`.
- Add `QuartetConsensus()`: a consensus tree minimizing the summed symmetric
  quartet distance to the input trees, found with a greedy add-and-prune
  heuristic (Takazawa et al. 2026).  The C++ core is self-contained.
- `Loose()`: loose (semi-strict / combinable-component) consensus.
- `Greedy()`: greedy (extended majority-rule) consensus.
- `MajorityPlus()`: majority-rule (+) consensus.
- `Frequency()`: frequency-difference consensus.
- `Average()`: distance-based average consensus.
- Add Billera-Holmes-Vogtmann (BHV) treespace summaries with branch lengths:
  `BHVDistance()` (Owen-Provan geodesic distance), `BHVPairwiseDistances()`,
  `BHVMean()` (iterative Fréchet mean) and `BHVVariance()`.  The geodesic core
  is implemented in C++.
- Add `AGENTS.md`, `README.md`, and source-paper references.
