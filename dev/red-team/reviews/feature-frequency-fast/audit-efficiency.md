# Skeptical Efficiency Audit: Frequency() / FDCT freqdiff port @ main 6b6f0e3

Reviewer: external (source-only; package NOT installed/run). Base = HEAD (already merged to main).

## One-line verdict
Efficiency claim TRUE (genuine O(kn log n) calc_w_knlogn, boost-free, no O(n^2) pass).
Correctness/rigor claim WEAKER THAN ADVERTISED but NOT a hidden divergence.
Block-class item: F3 (O(depth) recursion -> stack overflow on in-domain large caterpillars), shipped unfixed behind a tautological n=2000 guard.

## The "property-validated not FACT-exact" mystery -- SOLVED, user's hypothesis NOT supported
- The user feared "couldn't match the reference, so weakened the check." NO source evidence supports this.
- The label "property-validated; Boost oracle pending" was set by the human (commit 6cbc0d7, 2026-06-01) BEFORE the C++ port + oracle existed. It meant "lattice-property validated, oracle still pending."
- The port (e981693, 2026-06-02) built the oracle and its commit message claims "exact MATCH vs freqdiff.exe on all fixtures, incl. two new n>60 cases (n80/k31 congruent, n80/k20 incongruent)." The "Boost oracle pending" suffix was dropped; the "property-validated" tag was NEVER upgraded to "oracle-exact".
- So the downgrade is a STALE / UNDER-CLAIMED label relative to the author's own recorded evidence -- not a rationalized divergence. The tie-break-diagnosis machinery in check-freqdiff.R is a category error for a unique-output method but predates nothing suspicious; it is defensive scaffolding, not an escape hatch.

## Why "property-validated" is nonetheless HONEST (the real, defensible reason)
NOT "manual vs CI" -- check-oracle.R (Loose/MajorityPlus) is equally manual-only yet those are "FACT-exact". The true distinction is SHIPPED-TEST STRENGTH:
- Loose CI test pins the EXACT consensus on a nontrivial mixed fixture (test-loose.R:86-103, expect_setequal vs a hand-computed tree).
- Frequency CI tests assert full split-set identity ONLY on the two degenerate cases (all-agree->input test-frequency.R:50-54; all-conflict->star NSplits==0). EVERY nontrivial mixed case gets only the lattice envelope majority<=f<=greedy (lines 35-36, 115-116) or single-split membership (s12 %in% fs). For a UNIQUE-output method the full split set is well-defined and cheap to pin at small n -- and no CI test does it. A regression that gets a split wrong while staying inside the envelope passes CI silently. The only exact-vs-reference check (the n80 oracle match) lives in a non-CI dev script whose pass is recorded solely in a commit message.

## Findings
F1 (BLOCKER, efficiency-in-domain): O(depth) recursion unfixed. fix_tree_supp (225/229), build_taxas_ranges_supp (249/257), eulerian_walk (424), compute_m (909), newickInto (1362). Depth = n on caterpillars; ~n>1e4 overflows the stack. Test at test-frequency.R:307-320 caps n=2000 and its own comment says it "should complete today" / guards "once the iterative-rewrite fix lands" -- the fix did NOT land; tautological guard (passes pre- and post-"fix"). Crashes on in-domain empirical phylogenomic inputs.
F2 (MAJOR, rigor): No CI exact-identity test on a nontrivial fixture for a unique-output method. Add a hand-computed expect_setequal mirroring test-loose.R:86-103. Upgrade the roster label to "oracle-exact (n<=80 confirmed)" OR keep "property-validated" but say WHY (exact match is a non-CI manual spot-check).
F3 (MINOR, rigor): exact-match claim rests solely on a one-time manual commit-message assertion; check-freqdiff.R is invoked by nothing (not CI/make/tests). n80 match is unverifiable from source -- carve out as the separate empirical step (run check-freqdiff.R at n80+, plus larger n).
F4 (INFO): scratch-bound 2n/3n + radix 5n "not proven closed-form" (cons_frequency.cpp:838-843) -- documented, std::vector[] unchecked. Pre-existing review F5; unchanged.

## Checked and FINE (no algorithmic shortcut)
- LIVE weight pass IS calc_w_knlogn (run() line 1292 calls it; the bitset/node_bitvec_t path is genuinely omitted). Radix/divide-and-conquer label_nodes + radix_t counting sort present.
- alloc_int_matrix is NOT an n^2 LCA table: sole call site rmq_preprocess:357 passes block_size = int_log2(size) (line 321), <= sqrt(n) distinct block types -> sub-linear O(log^2 n) each. The advisor-flagged "n^2 LCA table" trap is a non-issue.
- Sparse-table RMQ (M, O(n)), Euler-tour LCA O(n)/tree, centroid-path decomposition (subpath_query_info_t / general_rmq), union-find skyline sweep -- all near-linear; no quadratic cluster comparison, no O(kn^2) weight calc.
- BOOST-FREE confirmed: git grep boost/dynamic_bitset over src/ returns only the explanatory comment (line 29). Dead struct truly dead, not reintroduced.
- The boost-compiled freqdiff.exe oracle does NOT ship: .Rbuildignore line 9 ^dev$ excludes dev/. Artifact-level boost-free claim holds.
- No TODO/FIXME/approximate/deferred/stub markers (only the documented scratch-bounds note).
- F1(tip-label) and F4(NA) from the prior review.md ARE fixed in merged source (.PrepareTrees selection.R:19,25-28).
