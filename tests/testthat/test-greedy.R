# Regression tests for the C++ greedy consensus (`greedyConsensusCpp`, ported
# from FACT `greedyConsensusFast`).  The rigorous correctness check is the dev
# oracle (`dev/oracle/check-oracle.R`, against `fact.exe`); these shipped tests
# assert tie-break-robust properties that must hold regardless of how equally
# frequent incompatible splits are resolved.

test_that("Greedy preserves the leaf set and contains the majority consensus", {
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    g <- Greedy(trees)
    expect_s3_class(g, "phylo")
    expect_setequal(TreeTools::TipLabels(g), labels)

    # Extended majority-rule contains majority-rule, for ANY valid tie-break:
    # the majority splits are most frequent and mutually compatible, so greedy
    # accepts them all before considering anything that could conflict.
    mSplits <- as.character(TreeTools::as.Splits(Majority(trees),
                                                 tipLabels = labels))
    gSplits <- as.character(TreeTools::as.Splits(g, tipLabels = labels))
    expect_true(all(mSplits %in% gSplits))
  }
})

test_that("Greedy returns a fully resolved tree when the inputs agree", {
  base <- ape::as.phylo(42, 12)
  trees <- structure(list(base, base, base, base), class = "multiPhylo")
  g <- Greedy(trees)
  # Identical binary inputs -> the consensus is the (binary) input itself.
  expect_equal(TreeTools::NSplits(g), TreeTools::NSplits(base))
  expect_setequal(
    unname(as.character(TreeTools::as.Splits(g, tipLabels = TreeTools::TipLabels(base)))),
    unname(as.character(TreeTools::as.Splits(base)))
  )
})

test_that("Greedy keeps a clade supported by a clear majority", {
  # Three trees: two share the rooted clade (t1, t2); the third conflicts.  The
  # clade occurs twice vs once, so greedy must retain it (no tie).
  trees <- c(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
  )
  labels <- TreeTools::TipLabels(trees[[1]])
  gSplits <- as.character(TreeTools::as.Splits(Greedy(trees), tipLabels = labels))
  wanted <- as.character(TreeTools::as.Splits(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"), tipLabels = labels))
  expect_true(all(wanted %in% gSplits))
})
