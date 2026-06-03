# Red-Team Review Session — 2026-06-02/03 (22:00–03:00 BST)

Automated adversarial review series, six slots.  Each slot picked the next
unreviewed functional area or open branch and wrote findings to
`dev/red-team/reviews/<area>/`.

---

## Slot summary

| Slot | Area | Verdict | Key findings |
|------|------|---------|--------------|
| 1 | **Loose consensus** (`feature/loose-fast`) | SHIP (nits only) | No correctness blocker; port is FACT-exact incl. polytomies. F1–F6 all tests/docs hardening: vacuous test cases, insertion path untested, dropped assert. Review: `reviews/feature-loose-fast/review.md` |
| 2 | **Quartet consensus** (`feature-quartet`) | CONDITIONAL SHIP | F1 (CRITICAL coverage): no per-state oracle; F2 (MAJOR): `pool.count` raw-row vs tree-membership semantic mismatch; F3 (MAJOR): brute-force oracle too narrow. Review: `reviews/feature-quartet/review.md` |
| 3 | **Frequency consensus / FDCT port** (`feature-frequency-fast`) | CONDITIONAL SHIP | F1 (MAJOR): missing tip-label validation in `.PrepareTrees()` — mismatched taxa cause silent wrong results or NULL deref in C++. Review: `reviews/feature-frequency-fast/review.md` |
| 4 | **BHV geodesic distance** (`claude/cranky-kare` branch) | **BLOCK** | F1 (CRITICAL): zero-length incompatible interior edge → infinite loop (confirmed repro, killed under timeout). F2 (MAJOR): wrong distance on degree-2 singleton internal nodes. Review: `reviews/claude-cranky-kare/review.md` |
| 5 | **MajorityPlus / FACT boundary** (`majorityplus-fact`) | REPRO ONLY — no full review | Repro script `repro-01-strict-boundary.R` confirmed that strict `>` (not `>=`) boundary is FACT-exact: twoToTwo drops, twoToOne keeps. No full `review.md` written this slot. |
| 6 | **R\* consensus + Local consensus** (`rstar-local-consensus`) | CONDITIONAL SHIP | F1 (MAJOR): no cross-tree taxon-label validation in `RStar()` and `Local()` — mismatched taxa silently corrupt triplet tally. F2 (MEDIUM): `buildParentDepth` preorder assumption undocumented/unguarded in `rstar.cpp`. F3 (MINOR): dead `seen` vector. G1/G2: no n=3 test, no mismatched-taxa test. Algorithm correctness (tally, plurality, strong-cluster assembly, Local DP) verified and passed. Review: `reviews/rstar-local-consensus/review.md` |

---

## Cross-cutting pattern

F1 in slots 3 and 6 (and implicitly in slot 2) are variants of the same
root cause: **no wrapper validates that all input trees share the same tip
label set**.  The `.PrepareTrees()` helper (used by Loose, Greedy, MajorityPlus,
Frequency) and the bespoke RStar/Local wrappers all derive the canonical label
vector from the first tree only, then call `RenumberTips(tr, labels)` on each
subsequent tree without a set-equality check.  A single fix in `.PrepareTrees()`
and parallel fixes in `RStar()` and `Local()` would close all three instances.
