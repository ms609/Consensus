# Regression tests for the C++ loose consensus (`looseConsensusCpp`, ported from
# FACT `looseConsensusFast`).  The rigorous correctness check is the dev oracle
# (`dev/oracle/check-oracle.R`, against `fact.exe`); the loose consensus is
# unique (compatible-with-all, no frequency tie-break), so it is FACT-exact at
# every n.  These shipped tests assert the defining property and the lattice
# relationships any correct implementation must satisfy -- and, crucially, pin
# cases that DISTINGUISH loose from strict and that exercise looseMerge's op == 1
# child-insertion path (which binary or identical inputs never reach).

# Canonical, orientation-safe split strings: PolarizeSplits fixes each
# bipartition's orientation, so the string is invariant to how the tree happens
# to be rooted (loose output and a reference tree may root differently).
canonSplits <- function(tr, labels) {
  if (TreeTools::NSplits(tr) == 0L) {
    return(character(0L))
  }
  unname(as.character(TreeTools::PolarizeSplits(
    TreeTools::as.Splits(tr, tipLabels = labels))))
}

test_that("Loose preserves the leaf set and contains the strict consensus", {
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    l <- Loose(trees)
    expect_s3_class(l, "phylo")
    expect_setequal(TreeTools::TipLabels(l), labels)

    # strict <= loose: every strict split survives in loose, for any correct
    # implementation.
    expect_true(all(canonSplits(Strict(trees), labels) %in% canonSplits(l, labels)))
  }
})

test_that("Loose is strictly more resolved than strict on compatible inputs", {
  # The defining advantage of loose (semi-strict) over strict: a split present in
  # SOME but not all inputs is kept when it is compatible with EVERY input.  Here
  # the two inputs are mutually compatible but unequally resolved, so strict
  # collapses to a star while loose keeps both nested clades.  The leaf-set and
  # compatibility tests cannot detect this (there loose == strict on every spec),
  # so this guards against a regression that silently aliases Loose to Strict.
  labels <- paste0("t", 1:6)
  trees <- c(
    ape::read.tree(text = "((t1,t2,t3),t4,t5,t6);"),
    ape::read.tree(text = "((t1,t2),t3,t4,t5,t6);"))
  expect_equal(TreeTools::NSplits(Strict(trees)), 0L)
  expect_gt(TreeTools::NSplits(Loose(trees)),
            TreeTools::NSplits(Strict(trees)))
  # Loose keeps both compatible (nested) clades {t1,t2} and {t1,t2,t3}.
  expect_setequal(
    canonSplits(Loose(trees), labels),
    canonSplits(ape::read.tree(text = "(((t1,t2),t3),t4,t5,t6);"), labels))
})

test_that("Every loose split is compatible with every input tree", {
  # The defining property: loose keeps exactly the splits contradicted by no
  # input, so each retained split must be compatible with each input tree.
  tested <- 0L
  for (spec in list(0:20, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    lSplits <- TreeTools::as.Splits(Loose(trees), tipLabels = labels)
    if (TreeTools::NSplits(lSplits) == 0L) {
      next
    }
    tested <- tested + 1L
    for (tr in trees) {
      tSplits <- TreeTools::as.Splits(tr, tipLabels = labels)
      compat <- as.matrix(TreeTools::CompatibleSplits(lSplits, tSplits))
      expect_true(all(compat))
    }
  }
  expect_gt(tested, 0L)  # at least one spec yielded a non-star consensus
})

test_that("Identical trees are their own loose consensus", {
  base <- ape::as.phylo(42, 12)
  labels <- TreeTools::TipLabels(base)
  trees <- structure(list(base, base, base, base), class = "multiPhylo")
  l <- Loose(trees)
  # Identical binary inputs -> the loose consensus is the (binary) input itself.
  expect_equal(TreeTools::NSplits(l), TreeTools::NSplits(base))
  expect_setequal(canonSplits(l, labels), canonSplits(base, labels))
})

test_that("Loose resolves the op==1 insertion path to the exact topology", {
  # Binary or identical inputs keep newNodes == 0, so looseMerge's BEFORE/AFTER
  # child-insertion loops never run.  These trichotomies force an A-cluster
  # ({t4,t5}) to be inserted among a B node's children, exercising that path --
  # and unlike the polytomy property test below we pin the EXACT loose tree, so
  # an under-resolution bug on the insertion path is caught here in CI, not only
  # by the print-only oracle.  Expected loose consensus, computed independently:
  #   {t1,t2} (in input 3, compatible with {t1,t2,t3}), {t4,t5} (in inputs 1-2,
  #   nested in {t3,t4,t5}), {t6,t7,t8} (in all) ==> ((t1,t2),t3,(t4,t5),(t6,t7,t8)).
  labels <- paste0("t", 1:8)
  trees <- c(
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2),(t3,t4,t5),(t6,t7,t8));"))
  expect_setequal(
    canonSplits(Loose(trees), labels),
    canonSplits(ape::read.tree(text = "((t1,t2),t3,(t4,t5),(t6,t7,t8));"), labels))
})

test_that("Loose handles polytomous inputs", {
  # A second polytomy shape (nested), as a property check on top of the exact
  # check above: each retained split must be compatible with every input, and
  # the result is a real (non-star) consensus.
  trees <- c(
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3),t4,(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),t7,(t8,t9));"))
  labels <- TreeTools::TipLabels(trees[[1]])
  l <- Loose(trees)
  expect_s3_class(l, "phylo")
  expect_setequal(TreeTools::TipLabels(l), labels)
  lSplits <- TreeTools::as.Splits(l, tipLabels = labels)
  expect_gt(TreeTools::NSplits(lSplits), 0L)   # a real (non-star) consensus
  for (tr in trees) {
    tSplits <- TreeTools::as.Splits(tr, tipLabels = labels)
    expect_true(all(as.matrix(TreeTools::CompatibleSplits(lSplits, tSplits))))
  }
})

test_that("Loose drops a split that any input contradicts", {
  # Two inputs show the clade {t3,t4}; the third shows {t3,t5},{t4,t6} instead,
  # incompatible with {t3,t4}.  Loose keeps the clade present in ALL inputs
  # ({t1,t2}) but drops the more-frequent {t3,t4} because an input contradicts
  # it -- exactly where loose parts company with majority/greedy.  The earlier
  # version of this test fed inputs whose loose consensus was a STAR, making the
  # "dropped" assertion vacuously true for any implementation; this set yields a
  # non-star consensus (asserted below), so the drop is a meaningful check.
  labels <- paste0("t", 1:6)
  trees <- c(
    ape::read.tree(text = "((t1,t2),(t3,t4),(t5,t6));"),
    ape::read.tree(text = "((t1,t2),(t3,t4),(t5,t6));"),
    ape::read.tree(text = "((t1,t2),(t3,t5),(t4,t6));"))
  lSplits <- canonSplits(Loose(trees), labels)
  expect_gt(length(lSplits), 0L)  # NON-star: makes the assertions below meaningful
  kept <- canonSplits(ape::read.tree(text = "((t1,t2),t3,t4,t5,t6);"), labels)
  dropped <- canonSplits(ape::read.tree(text = "((t3,t4),t1,t2,t5,t6);"), labels)
  expect_true(all(kept %in% lSplits))      # {t1,t2} retained (present in every input)
  expect_false(any(dropped %in% lSplits))  # {t3,t4} dropped (an input contradicts it)
})
