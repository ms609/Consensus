# Regression tests for the C++ loose consensus (`looseConsensusCpp`, ported from
# FACT `looseConsensusFast`).  The rigorous correctness check is the dev oracle
# (`dev/oracle/check-oracle.R`, against `fact.exe`); the loose consensus is
# unique (compatible-with-all, no frequency tie-break), so it is FACT-exact at
# every n.  These shipped tests assert the defining property and the lattice
# relationships that any correct implementation must satisfy.

test_that("Loose preserves the leaf set and contains the strict consensus", {
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    l <- Loose(trees)
    expect_s3_class(l, "phylo")
    expect_setequal(TreeTools::TipLabels(l), labels)

    # strict <= loose: every strict split (in every tree) survives in loose
    # (contradicted by none), for any correct implementation.
    sSplits <- as.character(TreeTools::as.Splits(Strict(trees), tipLabels = labels))
    lSplits <- as.character(TreeTools::as.Splits(l, tipLabels = labels))
    expect_true(all(sSplits %in% lSplits))
  }
})

test_that("Every loose split is compatible with every input tree", {
  # The defining property: loose keeps exactly the splits contradicted by no
  # input, so each retained split must be compatible with each input tree.
  for (spec in list(0:20, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    lSplits <- TreeTools::as.Splits(Loose(trees), tipLabels = labels)
    if (TreeTools::NSplits(lSplits) == 0L) next
    for (tr in trees) {
      tSplits <- TreeTools::as.Splits(tr, tipLabels = labels)
      compat <- as.matrix(TreeTools::CompatibleSplits(lSplits, tSplits))
      expect_true(all(compat))
    }
  }
})

test_that("Identical trees are their own loose consensus", {
  base <- ape::as.phylo(42, 12)
  trees <- structure(list(base, base, base, base), class = "multiPhylo")
  l <- Loose(trees)
  # Identical binary inputs -> the loose consensus is the (binary) input itself.
  expect_equal(TreeTools::NSplits(l), TreeTools::NSplits(base))
  expect_setequal(
    unname(as.character(TreeTools::as.Splits(l, tipLabels = TreeTools::TipLabels(base)))),
    unname(as.character(TreeTools::as.Splits(base)))
  )
})

test_that("Loose handles polytomous inputs", {
  # Binary inputs never present a node with > 2 children, so they never exercise
  # looseMerge's BEFORE/AFTER child-insertion across a multi-child node.  These
  # polytomous inputs do.  Loose is deterministic, so assert the defining
  # property (each retained split compatible with every input) -- orientation-
  # safe via CompatibleSplits; the exact tree is pinned by the dev oracle.
  trees <- c(
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2),(t3,t4,t5),(t6,t7,t8));")
  )
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
  # Two trees show the clade (t1, t2); the third shows (t1, t3) instead, which
  # is incompatible with it.  Loose excludes (t1, t2) because an input
  # contradicts it -- exactly where loose parts company with majority/greedy,
  # which would keep the more frequent split.
  trees <- c(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
  )
  labels <- TreeTools::TipLabels(trees[[1]])
  lSplits <- as.character(TreeTools::as.Splits(Loose(trees), tipLabels = labels))
  contradicted <- as.character(TreeTools::as.Splits(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"), tipLabels = labels))
  expect_false(any(contradicted %in% lSplits))
})
