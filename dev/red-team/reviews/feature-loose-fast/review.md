# Review: feature/loose-fast @ b0c2f93
Base e7e79a8. Built BRANCH source into temp lib (clean install), ran oracle + shipped tests.

## Verdict: Ship after addressing nits (TEST/DOC hardening only — NO code change required).
Port is empirically FACT-exact incl. op==1/polytomy under both rooted flags. No correctness blocker.

## Findings (all test/doc)
- F1 minor: test-loose.R:75-90 vacuous — Loose() is a STAR there; any(x %in% character(0))==FALSE.
- F2 minor: test-loose.R:8-22 & 24-38 — loose==strict (3==3) every spec; can't tell Loose from Strict.
- F3 minor: test-loose.R:24-38 `if(NSplits==0) next` latent vacuity (not triggered today).
- F4 minor: only exact-topology test (40-50, idempotence) uses IDENTICAL inputs -> newNodes==0,
  insertion loops never run. Insertion path's exact-topology check lives ONLY in print-only oracle.
- F5 nit: dropped assert(cnt==ret.cnt) (loose.cpp:258) at cons_loose.cpp:274.
- F6 nit: check-oracle.R print-only (no stop()). Highest-value fix: promote one polytomy oracle
  case into tests/testthat with a pinned PolarizeSplits set (closes F4+F6).

## Empirically FACT-exact this run: binary n9/n10/n8 r0&1; polytomy trichotomy+nested r0&1; n=80,137. 0 fail.
cmpPol sound (plain cmp DIFFERS on polytomy = real orientation artifact). Faithfulness claims 1-7 verified.
