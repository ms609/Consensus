# Red-Team Findings — ConsTree (open items)

Non-trivial bugs / perf / hardening issues, filed **after verification**.
Trivial issues are fixed inline and not listed here.

Rows below are **carried forward** from the pre-tier-system reviews under
`reviews/<area>/` (each was verified in its originating review) and reconciled
against `main` on **2026-06-15**. Confidence varies by row — the reconciliation
pass (a single read-only sweep) was **wrong on 2 of 2 rows that were spot-checked**
(package-infra F1/F2, both already fixed), so its verdicts are not authoritative:

- **Spot-verified present on main 2026-06-15** (code read directly): C-001, C-002,
  C-003, C-006, C-007, B-001. These are real, currently-open (minor/info) items.
- **Flagged by reconciliation, not yet re-verified**: C-004, C-005 (INFO),
  A-001, A-002, and the T-001 test-quality cluster. Treat as candidates until the
  area's next `/red-team` round confirms with the finder/verifier pair.

Headline "MAJOR" items from the old reviews were re-checked and found **already
fixed** (see *Resolved on reconciliation*), so **no MAJOR correctness finding
currently survives**.

| id | severity | status | area | title | file:line — detail | source |
|----|----------|--------|------|-------|--------------------|--------|
| C-001 | MINOR | OPEN | 4 Loose | `anc[-1]` deref when root not kept | `src/cons_loose.cpp:~70` — `keep[A.root]=1` guard absent; out-of-range read on a non-kept root | reviews/cpp-memory-safety |
| C-002 | MINOR | OPEN | 7 R*/Local | `lcaDepth` walks past root −1 | `src/rstar.cpp:~170`, `src/local_consensus.cpp:~45` — no root sentinel check on disconnected/degenerate input | reviews/cpp-memory-safety |
| C-003 | MINOR | OPEN | 5 MajorityPlus | Leaf `goodLabel` invariant undocumented | `src/cons_majorityplus.cpp:~57` — invariant relied on but not asserted/separated | reviews/cpp-memory-safety |
| C-004 | INFO | OPEN | 1/10 | `pool.data` pointer stability | `src/Quartet.cpp:~212` — no post-insertion assertion that backing store didn't reallocate | reviews/cpp-memory-safety |
| C-005 | INFO/perf | OPEN | 1 BHV | O(n²) queue in `gtp_no_common` | `src/bhv.cpp:~167,197` — `vector::erase(begin())` in a loop; fine for current tree sizes | reviews/cpp-memory-safety |
| C-006 | MINOR | OPEN | 10 FACT | Dropped `assert(num==N)` malformed-tree guard | `src/fact_tree.cpp:~79` — silent acceptance of malformed edge tables | reviews/fact-tree-primitive |
| C-007 | MINOR | OPEN | 10 FACT | `buildTreeFromEdge` dead root-fallback arm | `src/fact_tree.cpp:~130` — unreachable ternary `0` arm; clarity/dead-code | reviews/fact-tree-primitive |
| B-001 | MINOR | OPEN | 1 BHV | `TOL=1e-10` drop-tolerance is scale-dependent | `src/bhv.cpp:24` — absolute tolerance; very large/small branch lengths could mis-drop edges | reviews/bhv-numerical-precision |
| B-002 | INFO | OPEN | 1 BHV | Fréchet `converged` flag tracks step size, not gradient | `R/BHV.R:~293` (documented) — may report convergence on a slow crawl | reviews/bhv-numerical-precision |
| A-001 | MINOR | OPEN | 8 Average | Non-logical `edgeLengths` falls through silently | `R/Average.R` — no `isTRUE()` validation of the flag | reviews/average-consensus |
| A-002 | INFO | OPEN | 8 Average | NA branch lengths → confusing ape error | `R/Average.R` — no pre-check; deferred ape error is opaque | reviews/average-consensus |
| T-001 | MINOR/quality | OPEN | 12 Tests | Test-suite hardening cluster (~20 gaps) | `tests/testthat/` — discriminating power (inclusion-only assertions T6/T9), heuristic-optimality claims undocumented (T13/T22), RNG drift on shared seed (T19), duplicated `cladeSet` helper (T18), and ~14 more | reviews/test-suite-quality |

## Resolved on reconciliation (2026-06-15) — recorded so they are not re-opened
- **package-infrastructure F1** (was MAJOR, "stale `TreeTools::CompatibleSplits`
  import blocks CRAN") — **FIXED** (86711f8); not present in `NAMESPACE`/`R/`.
- **package-infrastructure F2** (stale `ape::write.tree` in `R/local.R`) — **FIXED** (86711f8).
- **cranky-kare F1** (CRITICAL infinite loop on zero-length incompatible interior
  edge) — **does not reproduce on main**; `BHVDistance` returns √3 without
  hanging. The BLOCK was against the abandoned, never-merged `claude/cranky-kare`
  branch. Repro driver kept at `reviews/claude-cranky-kare/run-against-main-f1.R`.
- **cranky-kare F2** (MAJOR wrong distance on internal singleton) — **does not
  reproduce on main**; `d(B, collapse.singles(B)) = 0`. Driver: `run-against-main-f2.R`.
  (`R/BHV.R:30` `collapse.singles()` at the top of `.TreeToBHV` is the reason.)
- **cranky-kare F3** (MINOR opaque error on missing `edge.length`) — **FIXED**;
  `R/BHV.R:28` now `stop("Trees must have edge.length to compute BHV distances.")`.
- **cpp-memory-safety F1** (MAJOR unchecked depth-index) — **FIXED** (b396a6f).
- **test-suite-quality T1/T5/T8** (MAJOR vacuous lattice tests) — **FIXED** (35869cb).
- **R\*/Local, Quartet, Frequency tip-label validation** — **FIXED** (1d77249, e07f1d6).
