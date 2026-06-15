# Red-Team Focus Areas — ConsTree

Rotation table for the `/red-team` skill. One area per invocation, next in
sequence after `last_focus:` (bottom of `log.md`). `start_tier` defaults to
`sonnet` for every area — **maturity is measured, not assumed**; let recorded
yield earn an escalation (`sonnet` → `opus` → `fable`).

This rotation was (re)built 2026-06-15 to consolidate ~18 pre-tier-system
reviews that lived only under `reviews/<area>/review.md`. Those reviews are
folded into `log.md` as prior history; this table is the live rotation.

| # | Area | Files | start_tier | Key questions |
|---|------|-------|------------|---------------|
| 1 | **BHV geodesic & Fréchet mean** | `src/bhv.cpp`, `R/BHV.R` | sonnet | GTP min-weight vertex cover & ratio sequence correct? Zero-length interior / internal-singleton edges (cranky-kare repros — **confirmed handled on main 2026-06-15**, keep covered). `TOL=1e-10` scale-dependence (bhv-numerical F1, open). Fréchet `converged` flag tracks step size not gradient (F2). |
| 2 | **Quartet consensus** | `src/Quartet.cpp`, `R/Quartet.R` | sonnet | `pool.count` raw-row vs tree-membership semantics (fixed e07f1d6). State/oracle pinning. Heuristic global-optimality assertions in tests (T13). `pool.data` pointer stability under reallocation (cpp-mem F6). |
| 3 | **Frequency consensus (FDCT)** | `src/cons_frequency.cpp` | sonnet | Radix capacity overflow (fixed e150d69, 6n). Heavy-path heap O(n²) blow-up (fixed 76b0803). Tip-label validation via `.PrepareTrees` (fixed 1d77249). Deep-recursion / scratch-bounds. |
| 4 | **Loose / Greedy / Strict** | `src/cons_loose.cpp`, `src/cons_greedy.cpp` | sonnet | `anc[-1]` deref when root not kept (cpp-mem F3, open). Insertion (`op==1`) path. Multi-word bitset path (n>60/64). FACT-exactness incl. polytomies. |
| 5 | **MajorityPlus consensus** | `src/cons_majorityplus.cpp` | sonnet | Strict `>` (not `>=`) boundary is FACT-exact (repro confirmed, slot 5). Leaf `goodLabel` invariant undocumented (cpp-mem F2, open). Path-query left/right-deeper arms. |
| 6 | **Adams consensus** | `src/cons_adams.cpp`, `R/adams.R` | sonnet | Centroid-path O(kn log n) via radix sort (7fc9128). `buildNub` m≥2 invariant (6e6d1e7). Degree-1 empty guards. Shared-leaf-set validation (2cfe45f). |
| 7 | **R\* & Local consensus** | `src/rstar.cpp`, `src/local_consensus.cpp`, `R/rstar.R`, `R/local.R` | sonnet | Cross-tree tip-label validation (fixed 1d77249/9d62a3d). `lcaDepth` walking past root −1 on disconnected input (cpp-mem F4, open). `buildParentDepth` preorder assumption. Triplet tally / plurality / strong-cluster assembly. Local DP. |
| 8 | **Average consensus** | `R/Average.R` | sonnet | n<3 leaves guard (fixed b9e2cf4). `edgeLengths` non-logical fall-through (F3, open). NA branch lengths → opaque ape error (F4, info). `check.labels=FALSE` superset (fixed b9e2cf4). |
| 9 | **Transfer consensus** | `src/transfer_consensus.cpp`, `R/transfer.R` | sonnet | **UNREVIEWED.** Ported from TreeDist 2026-06-09/10 (9390bcd, 0a0a356) with OpenMP (1bbc14f). Correctness vs TreeDist oracle? Tip-label validation? Edge cases (identical trees, stars, n small)? OpenMP data races / reduction correctness? |
| 10 | **FACT primitive & C++ memory safety** | `src/fact_tree.cpp/.h`, cross-cutting `src/*.cpp` bounds, OpenMP | sonnet | Dropped `assert(num==N)` malformed-tree guard (fact F1, open). `buildTreeFromEdge` dead root fallback (F2). Depth-index bounds into `R[]`/`L[]` (hardened b396a6f). OpenMP (1bbc14f, NEW) thread-safety across consensus kernels. |
| 11 | **R wrappers & package infrastructure** | `R/wrappers.R`, `R/selection.R`, `R/RcppExports.R`, `NAMESPACE`, `DESCRIPTION` | sonnet | Rcpp arg-count ↔ registration match. Stale `@importFrom` (CompatibleSplits / write.tree both fixed 86711f8). Dispatch (`MR`/`Strict`/`Majority` aliases) tested on bare `phylo` vs list. |
| 12 | **Test suite quality** | `tests/testthat/` | sonnet | ~20 carried-forward hardening gaps (T2–T26 in reviews/test-suite-quality). Vacuous lattice tests (fixed 35869cb). Discriminating power (inclusion vs exclusion). RNG drift (shared seed, T19). Duplicated `cladeSet` helper (T18). |

## Notes
- **Area 9 (Transfer) is the only entirely unreviewed area** and the freshest
  code (with new OpenMP) — highest expected yield. `last_focus` is seeded at 8
  so a bare `/red-team` resumes there.
- Areas 1–8, 10–12 each carry at least one pre-tier-system review (see `log.md`);
  re-visits run at `sonnet` first under the new tier system to *measure* maturity.
- OpenMP (added 1bbc14f, 2026-06-09) is cross-cutting — concurrency questions
  belong to area 10 but touch any kernel that parallelises.
