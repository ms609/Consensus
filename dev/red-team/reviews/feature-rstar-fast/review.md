# Review: feature/rstar-fast @ 9f81e91 (base e7e79a8)

Scope: the R* consensus rewrite — `src/rstar.cpp` (tensor-free O(kn²)-memory
construction: per-tree O(1) LCA → similarity tally → MST single-linkage
candidates → strong-cluster filter/build), the cap removal in `R/rstar.R`, and
the supporting oracles.

Started as a background cynical review (agent `aef91099`), which hung on an
account session limit after completing the analytical pass and one large
empirical probe but **before** writing a verdict or finishing the brute-force
oracle. This document completes it: the analytical findings are the original
reviewer's (independently re-checked against the source); the empirical results
and the wrapper findings (W1–W3) are the continuation.

## Verdict
**Pass on the optimised algorithm; one MAJOR pre-existing wrapper bug to fix.**
The cluster-correctness core is exhaustively validated — analysis plus five
independent empirical oracles agree, and the output tree is provably unchanged
from the previous definition-exact build to n=200. No defect was found in the
new algorithm or its memory/cap changes. The residual issues all live in the
**R wrapper's input contract** (unchanged by this rewrite) and one defensive
self-check; W1 (silent wrong answer on ragged taxon sets) is worth fixing before
this is used on real multi-locus data.

## Status — all findings resolved in the commit that adds this review
- **W1** — `RStar()` now rejects trees on differing leaf sets, so the silent
  ragged-input case errors loudly.  (`main` independently added the same
  `setequal` tip-label check via commit `b9e2cf4`, merged in `ae6ca3b`; the
  duplicate-label half, W2, is the part that fix missed.)
- **W2** — duplicate tip labels now error (`anyDuplicated` guard; `setequal`
  alone ignores multiplicity).
- **W3** — the `@return` doc now explains that an unresolved-root consensus is a
  basal polytomy that `ape::is.rooted()` reports as unrooted.
- **C1** — the C++ LCA self-check now runs over every input tree, not just
  tree 0.
All R* gates were re-run green after the fixes (test-rstar 79/79; four pillars;
brute-force strong-cluster oracle 180/180; new-vs-legacy 216/216; full suite
clean).  The findings below are retained as the record of why each fix exists.

## Build provenance
Validated against a **forced clean recompile** of `rstar.cpp` from the committed
9f81e91 source (`rm rstar.o ConsTree.dll; R CMD INSTALL`), because an earlier
`INSTALL` reported `make: Nothing to be done` — i.e. it would have reused
whatever object the 21:26 profiling run left behind. All results below are from
the provably-matching binary (version 0.0.0.9008, cap absent).

## W1. Trees on differing leaf sets → silently wrong clusters — MAJOR
Where: `R/rstar.R:98-112` + `src/rstar.cpp:212`.
`RStar()` takes `labels <- TipLabels(trees[[1]])`, `RenumberTips(tr, labels)` on
every other tree, then passes tree[0]'s `n` to the C++ for **all** trees
(`buildTreeFromEdge(edge, n)`). RenumberTips errors when a tree carries an
*extra* label (`"... missing from tipOrder"`), but a tree whose labels are a
strict *subset* of `labels` passes silently. In the C++ the absent leaves' ids
then alias that tree's internal nodes, so triples involving a missing taxon are
scored from those nodes' LCA depths instead of abstaining.
**Proven, not theoretical.** `probe-w1-oracle.R` compares the C++ output on
ragged input against the natural "each tree abstains on triples touching a taxon
it lacks" semantics (independent brute force): they **diverge in 129 / 190
ragged trials (68%)**, with no error or warning. Divergence runs both ways —
spurious clusters the data do not support (RStar emits `{t1,t3,t5}` where only
`{t1,t5}` holds; `{t1,t5}` where the abstain answer is a star) *and* dropped
clusters (RStar under-resolves where the abstain answer has three). A single
hand case (`probe-w1-divergence.R`) can look correct by luck — when the missing
taxon happens to be outermost where it is mentioned the fabricated reading
agrees — which is exactly why the bug is dangerous; the oracle sweep removes the
luck. Whether one regards R* as *undefined* on unequal leaf sets or as obliged
to abstain, silently returning aliasing-dependent clusters is wrong: the method
must **fail loud** here, as it already does for supersets. Realistic trigger:
multi-locus gene trees with incomplete, non-identical taxon sampling. Fix:
before the C++ call require every tree to have the identical leaf set as
`trees[[1]]` (equal length AND `setequal` labels), `stop()` otherwise (also
forecloses a theoretical out-of-bounds read at extreme size gaps). Pre-existing
— the relabel logic predates this rewrite; the only `R/rstar.R` change here was
deleting the `n > 200` cap.

## W2. Duplicate tip labels → accepted, meaningless tree — MINOR
Where: `R/rstar.R:98`. `RStar(list((((a,a),b),c), ...))` returns a tree with two
tips labelled `a`; R* is undefined on duplicated labels. RAN. Fix: reject
`anyDuplicated(labels)` with a clear message. Pre-existing.

## W3. `@return` "rooted by construction" overstated for unresolved roots — MINOR
Where: `R/rstar.R:61-62`. When the triplets do not resolve the root, the correct
R* output is a star or basal polytomy; `ape::is.rooted()` then returns FALSE
(root degree > 2). RAN: independent random trees, a star input, and a
basal-polytomy input all return `is.rooted = FALSE` trees. The trees are
*correct*; the doc claim is not literally true for these (common) cases.
Fix: soften the wording (note the consensus may be a basal polytomy that ape
reports as unrooted when the root is unresolved), or trivially root the output.

## C1. Stage-0 LCA self-check covers only tree 0 — MINOR (reviewer's note)
Where: `src/rstar.cpp:222-239`. The O(1)-LCA-vs-ancestor-walk self-check (n≤50)
runs on `edgeList[0]` only. All k trees share the identical `fillLcaDepths` path,
so a structural RMQ bug would surface on tree 0 with high probability across the
suite — a reasonable defensive check, not a guarantee. Optional: rotate the
checked tree by trial index, or document the rationale.

## Verified sound (the optimised core — analysis, re-checked against source)
- **Candidate generation.** MST maximum-spanning-tree single-linkage produces a
  laminar **superset** of the old per-threshold connected components; spurious
  binary intermediates are filtered by the strong test. Empirically a superset
  with **zero** missed threshold components over 20000 random + 40000
  randomized-tie-order similarity matrices (`probe-candidates.R`,
  `probe-tieorder.R`, both reproduced). Output is therefore tie-order invariant.
- **Min-side verification** (`sab - inside == need`). `s(a,b)` excludes a,b (the
  tally is over distinct i<j<l), and the Stage-1 increment rule and `inRmaj` use
  the **identical** strict-plurality decision `cab>0 && cab>cax && cab>cbx`, so
  the derived outside-count is exact. Airtight.
- **`closePair` vs the old triplet logic** — bit-identical by case; the final
  `return -1` is genuinely unreachable on any realizable LCA-depth triple (a
  (5,3,5)-type triple cannot occur in a tree). `# nocov` justified.
- **Integer overflow** — every `D[...]` / `s[...]` index casts the first factor
  to `size_t` before multiplying; `D` and `s` allocations cast `k`/`i` first.
  No pre-promotion `int` overflow in the hot indices.
- **`outChildren` sized 2n** — a laminar family on n leaves has ≤ n−1 non-trivial
  members, +root ≤ 2n−1 < 2n. Tight but correct; `memo` sized to match.
- **Prefilter `sab ≥ need`** — a true necessary condition for a strong cluster;
  cannot reject a genuine one.
- **Child-certification amortization** — within-block pairs are skipped safely:
  blocks form only on acceptance, so each accUF block is a previously-accepted
  strong cluster C ⊆ A and {x∉A} ⊆ {x∉C}; single-linkage guarantees children are
  processed before parents.
- **k=0** (`rStarConsensus(list(), n)`) — `s` is n×n, independent of k, zero-init;
  the D/closePair/inRmaj paths never run; result is a star. No OOB.

## Empirical validation completed (all green, on the clean-rebuilt binary)
- Four property pillars (`check-rstar.R`): identity 0/40; congruent == aho-build
  (n=4..16); strict-refinement 0/150; majority-refinement (Lemma 1.1) 0/400.
- New-vs-legacy exact clade diff (`check-vs-legacy.R`): **216/216** to n=200;
  runs at n>200 (cap gone).
- In-repo brute-force strong-cluster gate (`check-strong-clusters.R`, independent
  `ape::mrca` path): **180/180** — binary, partly-resolved, mixed-set batteries.
- Reviewer's second independent brute-force oracle (`probe-oracle-fast.R`, seed
  42, n=8..10): **250/250**.
- Uncapped-regime refinement (`probe-bign.R`, NEW): R* refines strict AND
  majority at n=250/300/500 (k=12); timings 0.19 / 0.22 / 0.76 s — the practical
  win of the lifted cap, confirmed.

## Wrapper behaviour confirmed correct
Extra-label / superset trees error cleanly; same-set-but-permuted tip vectors are
realigned correctly; `multiPhylo`, a list with a `NULL` element, a non-list, k=0
(→NULL), k=1 (→the tree), and n<3 (→the tree) are all handled sanely.
