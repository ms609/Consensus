# Red-team review: test suite quality
**Date:** 2026-06-03
**Scope:** all files in `tests/testthat/`
**Prior reviews consulted:** feature-loose-fast, feature-frequency-fast, feature-quartet, rstar-local-consensus (majorityplus-fact has no review.md yet)
**Focus:** test quality — tautologies, weak fixtures, missing negative tests, structural problems. Coverage gaps already identified in prior reviews are not re-listed; this review concentrates on the quality of the *existing* tests.

---

## test-selection.R

### T1 — Lattice-ordering test is vacuous on the star (MAJOR)
`tests/testthat/test-selection.R:12-23`

The test asserts `all(s %in% l)`, `all(l %in% g)`, etc. across four methods on
`ape::as.phylo(0:20, 9)`. If every method collapses to the star (zero splits),
all four split-string sets are `character(0)` and `all(x %in% character(0))`
returns `TRUE` vacuously. The test would pass even if Loose, Greedy, and
MajorityPlus all returned the star topology when the real consensus should be
non-trivial.

**What a real bug would slip past:** An off-by-one in the frequency threshold
that makes Greedy drop every split (returning the star) would still satisfy
`all(s %in% g)` because `s` is `character(0)`. Neither lattice direction would
catch it.

**How to harden:** add `expect_gt(length(g), 0L)` (or `length(l) > 0`) to
assert at least one of the sets is non-empty, so the ordering checks are not
trivially satisfied. The `test-loose.R:58-74` block already demonstrates the
right pattern (`tested <- 0L; ... expect_gt(tested, 0L)`).

### T2 — `expect_s3_class` assertion cannot distinguish a wrong-topology phylo (MINOR)
`tests/testthat/test-selection.R:43`

```r
expect_s3_class(MajorityPlus(trees), "phylo")
```
Any function that returns *a* `phylo` — even the star or the wrong tree — would
pass. The assertion lives inside a block that also checks `all(m %in% mp)`, so
the topology is partially pinned, but the class check alone adds zero
discriminating power and could mislead future readers into thinking it is a
meaningful property test.

**Recommendation:** delete the class check or replace it with the split-string
subset assertion that is already present one line above.

---

## test-loose.R

### T3 — `NSplits` count check cannot verify *which* splits are present (MINOR)
`tests/testthat/test-loose.R:82-83`

```r
expect_equal(TreeTools::NSplits(l), TreeTools::NSplits(base))
```
The identical-input idempotence test checks that the split *count* equals the
input's count, but a wrong tree with the same number of splits would pass. The
companion `expect_setequal(canonSplits(...))` on line 83 does pin the exact
splits, so the count check is redundant. However, if the `expect_setequal` were
removed in a future refactor, the count check would give false confidence.

**Recommendation:** delete the redundant `NSplits` count check and rely solely
on the `expect_setequal` for idempotence.

### T4 — `"Every loose split is compatible"` skips specs where the result is the star, but does not flag when all specs are skipped (MINOR)
`tests/testthat/test-loose.R:55-74`

```r
if (TreeTools::NSplits(lSplits) == 0L) {
  next
}
tested <- tested + 1L
```
The `tested > 0` guard at line 73 confirms at least one spec produced a
non-star consensus. However, if exactly one spec passes the guard and the others
all skip, only one spec's compatibility is checked. The comment says "three
specs" but two of `c(0:20, c(0,0,0,1,2,53,99), 1:15)` may produce stars for
this n=9 ensemble, leaving the property effectively verified on one case.

**Recommendation (low priority):** replace `tested > 0` with
`expect_gte(tested, 2L)` to ensure at least two specs exercise the compatibility
check.

---

## test-greedy.R

### T5 — Lattice ordering check is vacuous when majority is the star (MAJOR)
`tests/testthat/test-greedy.R:7-23`

```r
expect_true(all(mSplits %in% gSplits))
```
For four out of four `spec` values, the majority-rule consensus may be the star
on these random inputs (no split appears in >50% of the `as.phylo` ensemble
for that seed). `mSplits` is then `character(0)` and the assertion is vacuously
`TRUE` for any `g`, including the star. The companion `expect_s3_class(g, "phylo")`
and label check similarly cannot catch a wrong-topology result.

**What a real bug would slip past:** A Greedy implementation that returns the
star unconditionally would pass all three assertions in the loop body.

**How to harden:** mirror the pattern in `test-greedy.R:37-50` where a specific
fixture pins a majority-present clade. Add `expect_gt(length(gSplits), 0L)` for
at least one spec, or verify that the majority is non-empty by construction
before entering the loop body.

### T6 — `Greedy keeps a clade` test checks inclusion but not exclusion (MINOR)
`tests/testthat/test-greedy.R:37-50`

```r
expect_true(all(wanted %in% gSplits))
```
This checks that the two `wanted` splits appear in the greedy output, but does
not check that *no wrong* splits appear. If the greedy added both `(t1,t2)` and
`(t3,t4)` (correctly) but also added incompatible garbage splits, the test would
pass. The algorithm is sound in practice, but the test does not verify the full
topology.

**Recommendation (low priority):** add `expect_setequal(gSplits, wanted)` to
pin the full result for this specific fixture, as is done in test-loose.R:51-52.

---

## test-majorityplus.R

### T7 — Random-input scale test only checks the lattice relation, not the defining count rule (MINOR)
`tests/testthat/test-majorityplus.R:85-115`

The `"MajorityPlus handles larger perturbed input"` block asserts:
1. Result is a `phylo`.
2. Tip labels are correct.
3. Every majority split is in the MajorityPlus output.

Assertion 3 is exactly the lattice relation `majority ⊆ majorityPlus`. It would
be satisfied by any implementation that returns `Greedy(trees)` (since
`majority ⊆ greedy` always). The defining discriminator for MajorityPlus —
that a split displayed by more trees than contradict it is *kept* even when
below majority — is only exercised by the small exact fixtures (lines 49-83).
The larger perturbed fixture adds no discriminating power beyond a crash/class
check.

**What a real bug would slip past:** An implementation that always returns
`Greedy` would pass all assertions in this block.

**Recommendation:** for the `n=20` case, also verify that the result has *fewer*
splits than `Greedy(trees)` on some random seed where the two are known to
differ, or add a second exact fixture at `n ≈ 15` that pins the full split set
against an independently computed reference.

---

## test-frequency.R

### T8 — Lattice-ordering loop is vacuous when both bounds are the star (MAJOR)
`tests/testthat/test-frequency.R:10-24`

Same structural problem as T1 and T5: `all(mSplits %in% fSplits)` and
`all(fSplits %in% gSplits)` are both vacuously `TRUE` when the majority is a
star and the frequency-diff collapses to the star. The loop iterates four
`spec` values; for some of them the majority is the star, making the
`majority ⊆ frequency` check trivially pass.

**What a real bug would slip past:** A Frequency implementation that returns the
star for all inputs would pass every assertion in this loop (class, labels,
majority ⊆ frequency, frequency ⊆ greedy) for any spec where the real frequency
consensus also happens to be a star.

**How to harden:** add `expect_gt(length(fSplits), 0L)` for at least one spec,
or pick a spec where the frequency consensus is known to be non-trivial and
assert it explicitly.

### T9 — `"keeps a clade"` test checks inclusion but not exclusion of the conflicting split (MINOR)
`tests/testthat/test-frequency.R:39-52`

```r
expect_true(all(wanted %in% fSplits))
```
This verifies that the two splits from the majority-agreement trees survive, but
does not check that the conflicting split from the third tree is *absent*. The
frequency-difference rule should drop `(t1,t3)` (freq 1) because it conflicts
with `(t1,t2)` (freq 2). If a bug caused both to be retained (making an
incompatible tree), the test would pass.

**How to harden:** compute `conflicting_splits` from
`"((t1, t3), (t2, t4), t5);"` and add:
```r
expect_false(any(conflicting_splits %in% fSplits))
```

### T10 — Congruent ensemble in scalability test uses tip-label swaps that may produce a trivially star result (MINOR)
`tests/testthat/test-frequency.R:83-87`

```r
congruent = structure(lapply(seq_len(18), function(i) {
  tr <- base
  ij <- sample.int(40L, 2L)
  tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
  tr
}), class = "multiPhylo")
```
Swapping two tip *labels* changes which tip is called what, but every relabelled
tree has the same topology as `base` from the algorithm's perspective (the
splits are the same bipartitions, just with different taxa names). If any two
trees in this ensemble disagree on a split due to a label swap — which they will
— the frequency consensus may collapse in unexpected ways. The test only checks
`majority ⊆ frequency ⊆ greedy`, which passes trivially for the star. The regime
being tested (near-congruent, "rich surviving split pool") is not validated as
actually non-trivial.

**How to harden:** assert `expect_gt(length(fSplits), 0L)` for the congruent
ensemble to confirm the "rich surviving split pool" claim in the comment.

---

## test-BHV.R

### T11 — `BHVDistance(t1, t1)` tests identity of indiscernibles via the same object (MINOR)
`tests/testthat/test-BHV.R:52`

```r
expect_equal(BHVDistance(t1, t1), 0)
```
`t1` is the same R object on both sides. The implementation may have a fast-path
`if (identical(t1, t2)) return(0)` that bypasses the actual distance computation.
A bug in the geodesic computation for same-topology trees with different internal
representations would not be caught here.

**Recommendation (low priority):** add a test with two independently parsed
trees with identical Newick strings to confirm the distance is 0 without
relying on object identity. The existing `test-BHV.R:29-38` round-trip test
partially covers this via `BHVDistance(tree, back) == 0`, which is stronger.

### T12 — `BHVVariance` self-mean equality uses tolerance `1e-2` without documentation (MINOR)
`tests/testthat/test-BHV.R:150-153`

```r
expect_equal(v, BHVVariance(trees, mean = BHVMean(trees)), tolerance = 1e-2)
```
`BHVMean` is stochastic (documented), so tolerance is necessary. However, `1e-2`
is large relative to typical variance values for 6-tip trees with random branch
lengths in `[0.1, 1]`. A systematic 1% error in the variance computation would
never be detected.

**Recommendation:** document why `1e-2` is the right tolerance, or run
`BHVMean` with a fixed seed on both calls to allow a tighter tolerance.

---

## test-Quartet.R

### T13 — Brute-force n=5 oracle asserts *global optimality* for a heuristic algorithm (MAJOR — known)
`tests/testthat/test-Quartet.R:187-226`

Prior review (feature-quartet F3) identified this. The assertion
`expect_lte(qd_qc, best_loss)` claims the heuristic finds the global optimum.
This is a valid property to test *if* the algorithm is provably exact for small
n, but neither the code nor the test documents this. The real risk: a future
refactor slightly degrades the heuristic on this particular input, causing a
test failure that looks like a correctness regression but is actually just
heuristic variation — or a heuristic bug that degrades quality is not caught
because the optimum happens to be trivially achievable.

**Status:** known from prior review; re-flagged because it has not been
addressed. Fix: either (a) document (and cite) why the algorithm is exact for
n=5, or (b) pin the output as `expect_setequal(splitSet(qc), pinned_splits)`.

### T14 — `expect_lte(NSplits(qc), 0)` uses a weak bound that admits negative counts (MINOR)
`tests/testthat/test-Quartet.R:26`

```r
expect_lte(NSplits(qc), 0)
```
`NSplits` returns a non-negative integer; `<= 0` is equivalent to `== 0` but
looks like a weak bound. A negative return value (indicative of an internal
error) would pass this test. Use `expect_equal(NSplits(qc), 0L)` to make the
intent unambiguous and catch anomalous negative counts.

### T15 — `Quartet different init strategies` tests only class and tip count — no topology content (MINOR)
`tests/testthat/test-Quartet.R:41-60`

```r
expect_s3_class(qc_empty, "phylo")
expect_equal(length(qc_empty$tip.label), 8)
```
Both assertions would pass if all three strategies returned the star tree or any
other valid phylo. The test gives zero information about whether the three init
strategies produce plausible results.

**What a real bug would slip past:** An init strategy that resets the greedy to
the star at every call would pass both assertions.

**How to harden:** add `expect_gt(NSplits(qc_maj), 0L)` for the majority-rule
init (the majority rule on 15 random 8-tip trees is very unlikely to be a star),
and verify that `init = "majority"` yields at least as many splits as `Majority(trees)`.

### T16 — `Quartet minimizes quartet distance` compares against majority-rule only (MINOR)
`tests/testthat/test-Quartet.R:73-95`

```r
expect_lte(qd_qc, qd_mr)
```
This asserts that the quartet consensus is no *worse* than majority-rule under
the quartet distance. A result equal to majority-rule (e.g., because the greedy
gets stuck) would pass. The test does not verify that the quartet objective is
actually minimised relative to the star.

**Recommendation:** also assert `qd_qc < qd_star` where `qd_star` is the star
tree's quartet distance, to confirm the greedy added at least one useful split.

---

## test-adams.R

### T17 — `"Adams can recover a clade present in no input tree"` has only one fixed fixture (MINOR)
`tests/testthat/test-adams.R:34-43`

The test verifies a single pair `(t1, t2)` and checks that `"b,c"` appears.
The Adams consensus algorithm's novel-clade behaviour depends on the intersection
structure of the path sets; a regression in a different intersection case (three
conflicting paths, or deeper nesting) would not be caught.

**Recommendation:** add one additional hand-crafted fixture with a different
topology structure where a different novel clade is expected.

### T18 — `cladeSet` helper is duplicated verbatim across three test files (STRUCTURAL)
`tests/testthat/test-adams.R:3-15`, `test-rstar.R:13-25`, `test-local.R:3-15`

All three test files define a `cladeSet` function with identical logic. A bug
in the helper would need to be fixed in three places. More importantly, if
`DescendantEdges` behaviour changes upstream, all three test files would
silently compute wrong clade sets — the `expect_setequal` assertions would fail
in hard-to-diagnose ways.

**Recommendation:** consolidate `cladeSet` into a shared `tests/testthat/helper-clades.R`
(testthat auto-sources `helper-*.R` files), so the logic lives in one audited place.

---

## test-rstar.R

### T19 — Brute-force oracle loop shares one `set.seed` across 12 trials — order-dependent (MODERATE)
`tests/testthat/test-rstar.R:133-142`

```r
set.seed(2024)
for (trial in seq_len(12)) {
  n <- sample(4:7, 1); k <- sample(2:7, 1)
  trees <- .alignTrees(lapply(seq_len(k), function(i) {
    tr <- ape::rtree(n, rooted = TRUE); tr[["edge.length"]] <- NULL; tr
  }))
  expect_setequal(cladeSet(RStar(trees)), .strongClusters(trees))
}
```
A single `set.seed(2024)` governs all 12 trials. If any preceding test consumes
a different amount of RNG (e.g., after a future test is inserted before this
block), all 12 trial inputs shift silently. A bug that only manifests for
specific `(n, k)` combinations could be hidden or revealed purely by RNG
position at test entry.

**How to harden:** seed each trial independently:
```r
for (trial in seq_len(12)) {
  set.seed(2024 + trial)
  n <- sample(4:7, 1); k <- sample(2:7, 1)
  ...
}
```

### T20 — `.strongClusters` oracle itself has no self-test (INFO)
`tests/testthat/test-rstar.R:57-68`

The brute-force oracle used to validate `RStar` is not itself tested against a
known worked example. A bug in `.tripletClose` or `.favoured` could cause both
the oracle and the production code to silently agree on a wrong answer.

**Recommendation (low priority):** add one test that calls `.strongClusters`
directly on the Degnan et al. example (lines 79-86) and checks that it returns
`c("a,b", "a,b,c")`. This pins the oracle's correctness independently.

### T21 — No test for `RStar` with mismatched taxon sets (MAJOR — known)
`tests/testthat/test-rstar.R` (absent)

Prior review (rstar-local-consensus F5) identified this gap. After the F1 fix
(taxon validation) is applied to `R/rstar.R`, a regression test is essential to
prevent silent reversion to the broken behaviour.

---

## test-local.R

### T22 — `"Local rooted and induced produce different clades"` relies on a runtime search for a discriminating fixture (MAJOR)
`tests/testthat/test-local.R:144-169`

```r
set.seed(42)
found_trees <- NULL
for (trial in seq_len(200)) {
  ...
  if (!is.null(cr) && !is.null(ci) && !setequal(cr, ci)) {
    found_trees <- trees; break
  }
}
expect_false(is.null(found_trees), ...)
```
The test *searches at runtime* for an input where MinRLC ≠ MinILC. Three problems:

1. **Fragile failure message:** if a future commit makes both variants agree on
   all 200 tried inputs (e.g., a systematic bug that makes both variants
   identical), `expect_false(is.null(found_trees))` fails with "Expected to
   find trees where rooted != induced." This reads like a fixture problem but
   actually indicates an algorithm regression. The failure gives no information
   about which commit caused it.

2. **RNG drift breaks the fixture silently:** the comment says "seed=42, trial
   52 in random scan", but if a new package version changes `ape::rtree`'s RNG
   consumption, a different trial wins and the pinned "discriminating" case
   changes without warning. The final `expect_false(setequal(cr, ci))` would
   still pass (a different input still discriminates), but the property being
   tested has quietly shifted.

3. **Expensive:** 200 iterations of `Local(..., "rooted")` + `Local(...,
   "induced")` inside a unit test loop is the most expensive block in the suite.

**How to harden:** extract the discriminating trees found at seed=42, trial ~52
and hard-code them as the fixture. This makes the test deterministic, fast,
self-documenting, and immune to RNG drift. The comment "seed=42, trial 52 in
random scan" already points to the specific input; it just needs to be extracted.

### T23 — No test for `Local` with mismatched taxon sets (MAJOR — known)
`tests/testthat/test-local.R` (absent)

Same gap as T21, identified in prior review (F5). After the F1 fix to `R/local.R`,
a regression test is needed.

---

## test-Average.R

### T24 — Multiple fixtures use `ape::rtree` instead of `TreeTools::RandomTree` (STRUCTURAL)
`tests/testthat/test-Average.R:46, 63, 80, 92, 120, 152`

Per project convention (memory note "Prefer RandomTree over rtree"),
`ape::rtree` introduces a degree-2 root node. In `Average()` the trees are
unrooted internally, so the degree-2 root is collapsed, but the fixture itself
is non-standard and inconsistent with conventions used elsewhere in the suite.
If `Average` ever changed to process trees without unrooting, the degree-2 root
could silently affect results.

**Recommendation:** use `TreeTools::RandomTree(n, root = FALSE)` or
`ape::rmtree` (already used at line 103) for consistency.

### T25 — `"A single tree is its own average"` does not test `Average(list(tree), outgroup = ...)` (MINOR)
`tests/testthat/test-Average.R:92-98`

```r
expect_false(ape::is.rooted(Average(list(tree))))
expect_true(ape::is.rooted(Average(tree, outgroup = tree[["tip.label"]][[1]])))
```
The second assertion combines a bare `phylo` (not a list) with an `outgroup` —
a different code path than `Average(list(tree), outgroup = ...)`. If the
outgroup-rooting logic has a bug that only affects the single-tree-list path,
neither assertion catches it.

**Recommendation:** add `expect_true(ape::is.rooted(Average(list(tree), outgroup = tree[["tip.label"]][[1]])))`.

---

## test-wrappers.R

### T26 — `expect_equal` on `phylo` objects is fragile to internal representation changes (MODERATE)
`tests/testthat/test-wrappers.R:3-4, 7-9`

```r
expect_equal(Strict(trees), TreeTools::Consensus(trees, p = 1))
expect_equal(Majority(trees), TreeTools::Consensus(trees, p = 0.5))
```
`expect_equal` on `phylo` objects compares the full underlying list structure.
Two trees with the same topology can differ in internal node ordering, `edge`
matrix row order, or attribute presence, causing `expect_equal` to fail on a
*correct* wrapper or pass despite a topological difference if the representations
happen to coincide.

**How to harden:** use split-string comparison (as done throughout the rest of
the suite) rather than structural equality:
```r
labels <- TreeTools::TipLabels(trees[[1]])
expect_setequal(splitSet(Strict(trees), labels),
                splitSet(TreeTools::Consensus(trees, p = 1), labels))
```

---

## Verified strong

These patterns are well-designed and should be preserved or replicated:

- **Exact-topology pinning in test-loose.R:86-103** — pins the full split set
  against an independently computed reference; this is the gold standard for
  the suite. Every method should have at least one test at this level.
- **Non-emptiness guard in test-loose.R:58-74** — `tested <- 0L; ...; expect_gt(tested, 0L)`
  ensures the compatibility property check is not vacuously satisfied.
- **Defining-property boundary in test-majorityplus.R:49-64** — the 2-vs-1 and
  2-vs-2 exact boundary tests would catch an off-by-one in the `>` vs `>=`
  comparison; the most adversarial fixture in the suite.
- **`Frequency drops splits merely tied`** (`test-frequency.R:54-67`) — tests
  the exact property that distinguishes Frequency from Greedy and includes
  `expect_gt(NSplits(Greedy(...)), NSplits(f))` to confirm Greedy resolves where
  Frequency does not.
- **RStar brute-force oracle** (`test-rstar.R:132-142`) — comparing to an
  independent O(k n 2^n) implementation of the strong-cluster definition is
  excellent adversarial coverage for a method without a reference binary.
- **Local smoke case with exact pinned clades** (`test-local.R:92-109`) — exact
  topology pin for both method variants on a hand-designed conflicting input.
- **BHV Owen-Provan oracle** (`test-BHV.R:40-42`) — numeric oracle against a
  published exact value; would catch any regression in the geodesic length
  computation.
- **Mismatched leaf labels rejection** (`test-BHV.R:239-243`) — the only method
  currently testing the mismatched-label error path; other methods (RStar, Local,
  Frequency, Greedy, Loose) should follow this pattern once their F1 fixes are
  applied.
- **`Loose drops a split` non-emptiness guard** (`test-loose.R:139`) — load-bearing;
  see T3 for why removing it would silently re-introduce vacuity.
