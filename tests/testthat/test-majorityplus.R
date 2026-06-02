# Regression tests for the C++ majority-rule (+) consensus
# (`majorityPlusConsensusCpp`, ported from FACT `majorityPlusConsensus`).  The
# rigorous correctness check is the dev oracle (`dev/oracle/check-oracle.R`,
# against `fact.exe`, including an n > 60 exact block); these shipped tests pin
# the defining count-rule semantics -- a clade is kept iff displayed by strictly
# more trees than contradict it, with compatible-but-absent trees counting
# neither way -- and the lattice relationship majority subseteq majorityPlus.

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

# The polarised split string for an explicit in-group membership vector.
splitOf <- function(members, labels) {
  unname(as.character(TreeTools::PolarizeSplits(
    TreeTools::as.Splits(matrix(members, nrow = 1), tipLabels = labels)
  )))
}

test_that("MajorityPlus preserves the leaf set and contains the majority consensus", {
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99), 1:15)) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    mp <- MajorityPlus(trees)
    expect_s3_class(mp, "phylo")
    expect_setequal(TreeTools::TipLabels(mp), labels)
    # Every majority split is supported by more than half the trees, hence by
    # strictly more than contradict it, so it survives majority-rule (+).
    expect_true(all(splitSet(Majority(trees), labels) %in% splitSet(mp, labels)))
  }
})

test_that("Identical binary inputs are returned unchanged", {
  base <- ape::as.phylo(42, 12)
  trees <- structure(list(base, base, base, base), class = "multiPhylo")
  labels <- TreeTools::TipLabels(base)
  mp <- MajorityPlus(trees)
  expect_equal(TreeTools::NSplits(mp), TreeTools::NSplits(base))
  expect_setequal(splitSet(mp, labels), splitSet(base, labels))
})

test_that("A clade is kept iff displayed by strictly more trees than contradict it", {
  labels <- paste0("t", 1:5)
  display <- ape::read.tree(text = "((t1, t2), (t3, t4), t5);")  # shows {t1,t2}
  conflict <- ape::read.tree(text = "((t1, t3), (t2, t4), t5);") # breaks {t1,t2}
  s12 <- splitOf(c(TRUE, TRUE, FALSE, FALSE, FALSE), labels)     # the {t1,t2} split

  # Displayed by 2, contradicted by 1 (2 > 1): retained.
  twoToOne <- structure(list(display, display, conflict), class = "multiPhylo")
  expect_true(s12 %in% splitSet(MajorityPlus(twoToOne), labels))

  # Displayed by 2, contradicted by 2 (2 is NOT > 2): dropped.  This is the
  # boundary that separates majority-rule (+) from the strict count.
  conflict2 <- ape::read.tree(text = "((t1, t4), (t2, t3), t5);")
  twoToTwo <- structure(list(display, display, conflict, conflict2),
                        class = "multiPhylo")
  expect_false(s12 %in% splitSet(MajorityPlus(twoToTwo), labels))
})

test_that("MajorityPlus can keep a split the majority-rule consensus omits", {
  # {t1,t2}: displayed by 2 trees, contradicted by 1, and merely compatible with
  # 2 more (a basal polytomy leaving t1, t2, t3 unresolved -- it neither groups
  # nor splits t1,t2).  Displayed (2) > contradicted (1) so majority-rule (+)
  # keeps it; but 2 is not > 5/2, so the majority-rule consensus drops it.
  labels <- paste0("t", 1:6)
  display  <- ape::read.tree(text = "((t1, t2), (t3, t4), (t5, t6));")
  conflict <- ape::read.tree(text = "((t1, t3), (t2, t4), (t5, t6));")
  compat   <- ape::read.tree(text = "(t1, t2, t3, (t4, (t5, t6)));")
  trees <- structure(list(display, display, conflict, compat, compat),
                     class = "multiPhylo")
  s12 <- splitOf(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), labels)
  mpSplits <- splitSet(MajorityPlus(trees), labels)
  mSplits <- splitSet(Majority(trees), labels)
  expect_true(s12 %in% mpSplits)
  expect_false(s12 %in% mSplits)
})

test_that("MajorityPlus handles larger perturbed input (path-query at scale)", {
  # A base topology perturbed by random tip swaps produces clusters whose
  # boundaries fall on both sides of B's leaf order, so the merge exercises BOTH
  # arms of the C++ left-right path-query compatibility test -- not just the
  # left-deeper arm the tiny fixtures above happen to hit.  (The dev oracle
  # checks these inputs against fact.exe at n = 80/137; here we pin validity and
  # the lattice relation at a size the small fixtures never reach.)
  perturbed <- function(n, k, seed, nSwap = 3L) {
    set.seed(seed)
    labs <- paste0("t", seq_len(n))
    base <- TreeTools::RootTree(TreeTools::RandomTree(labs, root = TRUE), labs[[1]])
    swap <- function(tr) {
      for (s in seq_len(nSwap)) {
        ij <- sample.int(n, 2L)
        tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
      }
      TreeTools::RootTree(tr, labs[[1]])
    }
    structure(c(list(base), lapply(seq_len(k - 1L), function(i) swap(base))),
              class = "multiPhylo")
  }
  for (n in c(8L, 20L)) {
    trees <- perturbed(n, 15L, n + 7L)
    labels <- TreeTools::TipLabels(trees[[1]])
    mp <- MajorityPlus(trees)
    expect_s3_class(mp, "phylo")
    expect_setequal(TreeTools::TipLabels(mp), labels)
    # Every majority split survives majority-rule (+) (the lattice relation).
    expect_true(all(splitSet(Majority(trees), labels) %in% splitSet(mp, labels)))
  }
})

test_that("MajorityPlus returns a valid tree on conflicting input", {
  # Mutually incompatible equal-frequency splits all sit at the displayed ==
  # contradicted boundary, so none is retained: the result is the (valid) star.
  labels <- paste0("t", 1:5)
  trees <- structure(list(
    ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
    ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
  ), class = "multiPhylo")
  mp <- MajorityPlus(trees)
  expect_s3_class(mp, "phylo")
  expect_setequal(TreeTools::TipLabels(mp), labels)
  expect_equal(TreeTools::NSplits(mp), 0L)
})
