# Regression tests for the C++ frequency-difference consensus
# (`frequencyConsensusCpp`, ported from the FDCT_new freqdiff2.h near-linear
# algorithm of Jansson, Sung, Tabatabaee & Yang, 2024).  The rigorous
# correctness check is the dev oracle (`dev/oracle/freqdiff/check-freqdiff.R`,
# against `freqdiff.exe`); these shipped tests assert the definitional
# properties of the frequency-difference consensus, which -- unlike greedy -- has
# a unique output (a split is kept iff it is STRICTLY more frequent than every
# split that conflicts with it, so there is no tie-break freedom).

test_that("Frequency preserves the leaf set and lies between majority and greedy", {
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    f <- Frequency(trees)
    expect_s3_class(f, "phylo")
    expect_setequal(TreeTools::TipLabels(f), labels)

    mSplits <- as.character(TreeTools::as.Splits(Majority(trees), tipLabels = labels))
    fSplits <- as.character(TreeTools::as.Splits(f, tipLabels = labels))
    gSplits <- as.character(TreeTools::as.Splits(Greedy(trees), tipLabels = labels))
    expect_true(all(mSplits %in% fSplits))   # majority        <= frequency
    expect_true(all(fSplits %in% gSplits))   # frequency-diff   <= greedy
  }
})

test_that("Frequency returns the input when all trees agree", {
  base <- ape::as.phylo(42, 12)
  trees <- structure(list(base, base, base, base), class = "multiPhylo")
  f <- Frequency(trees)
  # Identical inputs: every cluster has weight k with no conflicts, so the
  # (binary) input itself is returned.
  expect_equal(TreeTools::NSplits(f), TreeTools::NSplits(base))
  expect_setequal(
    unname(as.character(TreeTools::as.Splits(f, tipLabels = TreeTools::TipLabels(base)))),
    unname(as.character(TreeTools::as.Splits(base)))
  )
})

test_that("Frequency keeps a clade strictly more frequent than its rivals", {
  # (t1, t2) occurs in two of three trees; the conflicting (t1, t3) in one.  The
  # frequency-difference rule retains (t1, t2): 2 > 1 over every conflicting split.
  trees <- c(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
  )
  labels <- TreeTools::TipLabels(trees[[1]])
  fSplits <- as.character(TreeTools::as.Splits(Frequency(trees), tipLabels = labels))
  wanted <- as.character(TreeTools::as.Splits(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"), tipLabels = labels))
  expect_true(all(wanted %in% fSplits))
})

test_that("Frequency drops splits merely tied with a conflicting rival", {
  # Two fully incompatible trees: every non-trivial split occurs exactly once, so
  # none is STRICTLY more frequent than its conflicting rivals -> frequency keeps
  # nothing (an unresolved star), whereas greedy's tie-break would retain a
  # resolution.  This is the property that distinguishes Frequency from Greedy.
  trees <- c(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
  )
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_equal(TreeTools::NSplits(f), 0L)
  expect_gt(TreeTools::NSplits(Greedy(trees)), TreeTools::NSplits(f))
})

test_that("Frequency scales to larger ensembles (centroid-path recursion)", {
  # n = 40 exercises code paths the small fixtures never reach -- the centroid-
  # path recursion, special-node insertion in the tree contraction, and (k near
  # n) the weight-compression branch -- guarding them against crashes/regressions
  # in CI (the dev oracle checks exact values).  Both regimes: an incongruent
  # ensemble (the harder filter/merge path) and a near-congruent one (a rich
  # surviving split pool that drives the cluster grafting).  Assertions are the
  # lattice invariants, which hold for any inputs (robust to RNG drift).
  set.seed(20L)
  base <- TreeTools::RandomTree(40L)
  ensembles <- list(
    incongruent = structure(
      lapply(seq_len(18), function(i) TreeTools::RandomTree(40L)), class = "multiPhylo"),
    congruent = structure(lapply(seq_len(18), function(i) {
      tr <- base
      ij <- sample.int(40L, 2L)
      tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
      tr
    }), class = "multiPhylo")
  )
  for (trees in ensembles) {
    labels <- TreeTools::TipLabels(trees[[1]])
    f <- Frequency(trees)
    expect_s3_class(f, "phylo")
    expect_setequal(TreeTools::TipLabels(f), labels)
    mSplits <- as.character(TreeTools::as.Splits(Majority(trees), tipLabels = labels))
    fSplits <- as.character(TreeTools::as.Splits(f, tipLabels = labels))
    gSplits <- as.character(TreeTools::as.Splits(Greedy(trees), tipLabels = labels))
    expect_true(all(mSplits %in% fSplits))   # majority        <= frequency
    expect_true(all(fSplits %in% gSplits))   # frequency-diff   <= greedy
  }
})

test_that("Frequency short-circuits the degenerate cases", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Frequency(tree), tree)           # a bare phylo is returned as-is
  expect_equal(Frequency(list(tree)), tree)     # a single-tree list, likewise
  triplet <- ape::as.phylo(0:1, 3)              # < 4 leaves: no consensus to take
  expect_equal(Frequency(triplet), triplet[[1]])
  expect_error(Frequency(5), "list of trees")
})
