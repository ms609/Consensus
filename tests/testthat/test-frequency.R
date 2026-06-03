# Regression tests for the C++ frequency-difference consensus
# (`frequencyConsensusCpp`, ported from the FDCT_new freqdiff2.h near-linear
# algorithm of Jansson, Sung, Tabatabaee & Yang, 2024).  The rigorous
# correctness check is the dev oracle (`dev/oracle/freqdiff/check-freqdiff.R`,
# against `freqdiff.exe`); these shipped tests assert the definitional
# properties of the frequency-difference consensus, which -- unlike greedy -- has
# a unique output (a split is kept iff it is STRICTLY more frequent than every
# split that conflicts with it, so there is no tie-break freedom).

# Canonical, order-independent split-string set (polarised for comparison).
splitSet <- function(tree, labels) {
  if (TreeTools::NSplits(tree) == 0L) {
    character(0)
  } else {
    unname(as.character(TreeTools::PolarizeSplits(
      TreeTools::as.Splits(tree, tipLabels = labels)
    )))
  }
}

test_that("Frequency preserves the leaf set and lies between majority and greedy", {
  # All four specs on n=9 produce non-empty majority and frequency splits (verified
  # empirically), so the lattice assertions are non-vacuous across all iterations.
  non_vacuous_checked <- FALSE
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
    if (length(fSplits) > 0L) non_vacuous_checked <- TRUE
  }
  # At least one spec must produce a non-star frequency consensus, confirming
  # that the lattice assertions above are not all trivially satisfied.
  expect_true(non_vacuous_checked)
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

# ---------------------------------------------------------------------------
# Coverage-gap tests added from red-team review slot 3
# ---------------------------------------------------------------------------

test_that("Frequency retains the majority clade in a k=5 three-way vote", {
  # Coverage gap 1: k >= 3 explicit frequency counts.
  # Three trees share clade {t1,t2}; two trees share the incompatible clade
  # {t1,t3}.  Frequency-difference rule: 3 > 2, so {t1,t2} must appear; 2 is not
  # strictly greater than 3, so {t1,t3} must not appear.
  labels <- paste0("t", 1:7)
  trees <- list(
    ape::read.tree(text = "((t1,t2),(t3,t4),(t5,t6,t7));"),
    ape::read.tree(text = "((t1,t2),(t3,t5),(t4,t6,t7));"),
    ape::read.tree(text = "((t1,t2),(t4,t6),(t3,t5,t7));"),
    ape::read.tree(text = "((t1,t3),(t2,t4),(t5,t6,t7));"),
    ape::read.tree(text = "((t1,t3),(t2,t5),(t4,t6,t7));")
  )
  f <- Frequency(trees)
  fs <- splitSet(f, labels)
  # The freq-3 clade {t1,t2} must be present.
  s12 <- splitSet(ape::read.tree(text = "((t1,t2),(t3,t4,t5,t6,t7));"), labels)
  expect_true(s12 %in% fs)
  # The freq-2 clade {t1,t3} must be absent.
  s13 <- splitSet(ape::read.tree(text = "((t1,t3),(t2,t4,t5,t6,t7));"), labels)
  expect_false(s13 %in% fs)
})

test_that("Frequency of all-star trees is unresolved", {
  # Coverage gap 2: every input is a star (zero non-trivial splits) so every
  # split has frequency 0 -- no split can strictly exceed any rival, and the
  # output must also be unresolved.
  n <- 7L
  tips <- paste0("t", seq_len(n))
  stars <- structure(
    replicate(5, ape::stree(n, tip.label = tips), simplify = FALSE),
    class = "multiPhylo"
  )
  f <- Frequency(stars)
  expect_s3_class(f, "phylo")
  expect_equal(TreeTools::NSplits(f), 0L)
})

test_that("Frequency of opposite caterpillars is a valid tree", {
  # Coverage gap 3: left- vs. right-leaning caterpillar on n=8 identical labels.
  # The two caterpillars are the same unrooted tree, so the result should be
  # fully resolved.  Assert the lattice invariants and structural validity;
  # do not assume a specific topology (unrooted identity is the key constraint).
  n <- 8L
  tips <- paste0("t", seq_len(n))
  t_left  <- ape::stree(n, type = "left",  tip.label = tips)
  t_right <- ape::stree(n, type = "right", tip.label = tips)
  trees <- list(t_left, t_right)
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_setequal(TreeTools::TipLabels(f), tips)
  # Lattice invariants must hold.
  fs <- splitSet(f, tips)
  ms <- splitSet(Majority(trees),  tips)
  gs <- splitSet(Greedy(trees),    tips)
  expect_true(all(ms %in% fs))
  expect_true(all(fs %in% gs))
})

test_that("Frequency handles mixed fully-resolved and star inputs", {
  # Coverage gap 4: an ensemble containing both binary trees and star topologies.
  # The star trees contribute no splits; the binary trees' splits may still be
  # retained if no input conflicts with them.  Result must be a valid phylo.
  labels <- paste0("t", 1:6)
  binary <- ape::read.tree(text = "((t1,t2),(t3,(t4,(t5,t6))));")
  star   <- ape::stree(6L, tip.label = labels)
  trees  <- list(binary, binary, star, binary, star)
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_setequal(TreeTools::TipLabels(f), labels)
  # Lattice: majority <= frequency <= greedy.
  fs <- splitSet(f,              labels)
  ms <- splitSet(Majority(trees), labels)
  gs <- splitSet(Greedy(trees),   labels)
  expect_true(all(ms %in% fs))
  expect_true(all(fs %in% gs))
})

test_that("Frequency produces a star when conflicting splits are exactly tied (k/2)", {
  # Coverage gap 5: k=4 trees, split A in trees 1-2, incompatible split B in
  # trees 3-4.  freq(A) = freq(B) = 2; neither is STRICTLY greater than its
  # rival, so both are deleted and the result is unresolved.
  trees <- list(
    ape::read.tree(text = "((t1,t2),(t3,t4),(t5,t6,t7));"),
    ape::read.tree(text = "((t1,t2),(t3,t5),(t4,t6,t7));"),
    ape::read.tree(text = "((t1,t3),(t2,t4),(t5,t6,t7));"),
    ape::read.tree(text = "((t1,t3),(t2,t5),(t4,t6,t7));")
  )
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_equal(TreeTools::NSplits(f), 0L)
})

test_that("Frequency drops the sole conflicting pair in a 2-tree input", {
  # Coverage gap 6: k=2 with one pair of incompatible clades ({t1,t2} vs
  # {t1,t3}).  Each appears once -- tied -- so neither can survive.  A clade
  # shared by both trees with no conflict ({t4,t5}) must be retained.
  labels <- paste0("t", 1:5)
  trees <- list(
    ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));"),
    ape::read.tree(text = "((t1,t3),(t2,(t4,t5)));")
  )
  f <- Frequency(trees)
  fs <- splitSet(f, labels)
  # Conflicting split {t1,t2} must not appear.
  s12 <- splitSet(ape::read.tree(text = "((t1,t2),(t3,t4,t5));"), labels)
  expect_false(s12 %in% fs)
  # Uncontested shared split {t4,t5} must appear.
  s45 <- splitSet(ape::read.tree(text = "((t4,t5),(t1,t2,t3));"), labels)
  expect_true(s45 %in% fs)
})

test_that("Frequency of k>>n identical trees matches the single-tree result", {
  # Coverage gap 7: 60 identical trees on n=8 tips drives the weight-compression
  # branch in filter() (k >> n triggers the quicksort path).  All splits have
  # maximum frequency with no conflicts, so the consensus must equal the input.
  base   <- ape::read.tree(text = "((t1,t2),(t3,(t4,(t5,(t6,(t7,t8))))));")
  labels <- TreeTools::TipLabels(base)
  trees  <- structure(replicate(60, base, simplify = FALSE), class = "multiPhylo")
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_equal(TreeTools::NSplits(f), TreeTools::NSplits(base))
  expect_setequal(splitSet(f, labels), splitSet(base, labels))
})

test_that("Frequency works on the minimum n=4 input", {
  # Coverage gap 8: four-tip trees do not short-circuit (the threshold is < 4)
  # and produce a valid result.  Two conflicting binary trees on 4 leaves yield a
  # star (the only non-trivial split in each is incompatible with the other).
  trees <- list(
    ape::read.tree(text = "((t1,t2),(t3,t4));"),
    ape::read.tree(text = "((t1,t3),(t2,t4));")
  )
  f <- Frequency(trees)
  expect_s3_class(f, "phylo")
  expect_setequal(TreeTools::TipLabels(f), paste0("t", 1:4))
  # Both non-trivial splits are tied and incompatible: the result is a star.
  expect_equal(TreeTools::NSplits(f), 0L)
})

test_that("Frequency retains a clade with non-symmetric incompatibility", {
  # Coverage gap 9: split A = {t1,t2,t3} appears in k-1=4 of k=5 trees; four
  # *distinct* splits, each conflicting only with A, appear once each.  Because
  # freq(A) = 4 strictly exceeds every rival's frequency of 1, A must survive.
  labels <- paste0("t", 1:8)
  trees <- list(
    ape::read.tree(text = "(((t1,t2,t3),(t4,t5)),(t6,t7,t8));"),
    ape::read.tree(text = "(((t1,t2,t3),(t4,t6)),(t5,t7,t8));"),
    ape::read.tree(text = "(((t1,t2,t3),(t4,t7)),(t5,t6,t8));"),
    ape::read.tree(text = "(((t1,t2,t3),(t5,t6)),(t4,t7,t8));"),
    ape::read.tree(text = "(((t1,t4),(t2,t5)),(t3,t6,t7,t8));")
  )
  f <- Frequency(trees)
  fs <- splitSet(f, labels)
  # The freq-4 clade {t1,t2,t3} must appear in the consensus.
  sA <- splitSet(ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6,t7,t8));"), labels)
  expect_true(sA %in% fs)
})

test_that("Frequency errors on mismatched tip labels", {
  # Coverage gap 10 (F1 guard): .PrepareTrees() validates the taxa set and
  # stops with a clear message before reaching C++.
  t1 <- ape::read.tree(text = "((a,b),(c,d));")
  t2 <- ape::read.tree(text = "((a,b),(c,e));")   # 'e' instead of 'd'
  expect_error(Frequency(list(t1, t2)), "tip labels")
})

test_that("Frequency silently drops NA entries in the tree list", {
  # Coverage gap 11 (F4 guard): Filter(inherits(., "phylo"), ...) in
  # .PrepareTrees() silently removes NA and other non-phylo objects.
  t <- ape::read.tree(text = "((a,b),(c,d));")
  expect_equal(Frequency(list(t, NA, t)), Frequency(list(t, t)))
})

test_that("Frequency completes on large caterpillar inputs without stack overflow", {
  # Coverage gap 12 (F3 regression guard): O(depth) recursion in fix_tree_supp,
  # eulerian_walk, compute_m, and newickInto can overflow the system stack for
  # caterpillar trees with n >> 10 000.  n = 2 000 is well below the empirical
  # crash threshold of the unfixed code and should complete today; it also serves
  # as a regression guard once the iterative-rewrite fix lands.
  n    <- 2000L
  tips <- paste0("t", seq_len(n))
  t_left  <- ape::stree(n, type = "left",  tip.label = tips)
  t_right <- ape::stree(n, type = "right", tip.label = tips)
  f <- Frequency(list(t_left, t_right))
  expect_s3_class(f, "phylo")
  expect_setequal(TreeTools::TipLabels(f), tips)
})
