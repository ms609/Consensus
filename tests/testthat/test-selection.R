# Canonical split-string set for a tree, for order-independent comparison.
splitSet <- function(tree, labels) {
  if (TreeTools::NSplits(tree) == 0L) {
    character(0)
  } else {
    unname(as.character(TreeTools::PolarizeSplits(
      TreeTools::as.Splits(tree, tipLabels = labels)
    )))
  }
}

test_that("Loose and Greedy respect the consensus lattice", {
  # Constructed fixture with non-empty split sets at each lattice level.
  #   - (a,b) is in ALL 6 trees  -> strict, loose, majority, and greedy all keep it
  #   - (c,d) is in 4/6 trees    -> loose, majority, and greedy keep it; strict does not
  # This guarantees strict (1 split) is a proper subset of loose/majority (2
  # splits), making the subset assertions non-vacuous even for a buggy
  # implementation that collapses every method to the star.
  t_all <- ape::read.tree(text = "((a,b),(c,d),e);")  # resolves (a,b) and (c,d)
  t_ab  <- ape::read.tree(text = "((a,b),c,d,e);")    # resolves only (a,b)
  trees <- c(t_all, t_all, t_all, t_all, t_ab, t_ab)
  labels <- TreeTools::TipLabels(trees[[1]])
  s <- splitSet(Strict(trees), labels)
  l <- splitSet(Loose(trees), labels)
  m <- splitSet(Majority(trees), labels)
  g <- splitSet(Greedy(trees), labels)
  # Pin the computed split sets so a star-returning implementation cannot pass
  # by vacuity (character(0) %in% character(0) is always TRUE).
  expect_setequal(s, "a b | c d e")
  expect_setequal(m, c("a b | c d e", "a b e | c d"))
  # Lattice ordering: strict <= loose <= greedy, majority <= greedy.
  expect_true(all(s %in% l))   # strict   <= loose
  expect_true(all(s %in% g))   # strict   <= greedy
  expect_true(all(l %in% g))   # loose    <= greedy
  expect_true(all(m %in% g))   # majority <= greedy
})

test_that("Identical trees are their own consensus", {
  tree <- ape::as.phylo(42, 9)
  trees <- structure(list(tree, tree, tree), class = "multiPhylo")
  labels <- TreeTools::TipLabels(tree)
  target <- splitSet(tree, labels)
  expect_setequal(splitSet(Loose(trees), labels), target)
  expect_setequal(splitSet(Greedy(trees), labels), target)
  expect_setequal(splitSet(MajorityPlus(trees), labels), target)
  expect_setequal(splitSet(Frequency(trees), labels), target)
})

test_that("MajorityPlus contains the majority-rule consensus", {
  trees <- ape::as.phylo(0:20, 9)
  labels <- TreeTools::TipLabels(trees[[1]])
  m <- splitSet(Majority(trees), labels)
  mp <- splitSet(MajorityPlus(trees), labels)
  expect_true(all(m %in% mp))
  # The retained splits define a valid (compatible) tree:
  expect_s3_class(MajorityPlus(trees), "phylo")
})

test_that("Frequency-difference lies between majority and greedy", {
  # All three specs on n=9 produce non-empty majority and frequency splits (verified
  # empirically), so the subset assertions below are non-vacuous.
  non_vacuous_checked <- FALSE
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99))) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    m <- splitSet(Majority(trees), labels)
    f <- splitSet(Frequency(trees), labels)
    g <- splitSet(Greedy(trees), labels)
    expect_true(all(m %in% f))   # majority        <= frequency-difference
    expect_true(all(f %in% g))   # frequency-diff   <= greedy
    expect_s3_class(Frequency(trees), "phylo")
    if (length(f) > 0L) non_vacuous_checked <- TRUE
  }
  # At least one spec must produce a non-star frequency consensus, confirming
  # that the lattice assertions above are not all trivially satisfied.
  expect_true(non_vacuous_checked)
})

test_that("A single tree, or fewer than four leaves, returns the input", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Loose(list(tree)), tree)
  expect_equal(Greedy(list(tree)), tree)
  triplet <- ape::as.phylo(0:1, 3)
  expect_equal(Loose(triplet), triplet[[1]])
})

test_that("Every selection method short-circuits a bare single tree", {
  # A bare `phylo` (not a list) triggers the `.PrepareTrees()` trivial-tree path
  # and is returned unchanged; exercised for the methods whose trivial branch is
  # otherwise untested.
  tree <- ape::as.phylo(1, 9)
  expect_equal(MajorityPlus(tree), tree)
  expect_equal(Frequency(tree), tree)
  expect_equal(Loose(tree), tree)
  expect_equal(Greedy(tree), tree)
})

test_that("Selection methods reject non-list input", {
  expect_error(Loose(5), "list of trees")
  expect_error(MajorityPlus("not a tree"), "list of trees")
})

test_that("selection methods error on mismatched tip labels", {
  t1 <- ape::read.tree(text = "((a,b),(c,d));")
  t2 <- ape::read.tree(text = "((a,b),(c,e));")
  expect_error(Loose(list(t1, t2)),     "tip label")
  expect_error(Greedy(list(t1, t2)),    "tip label")
  expect_error(Frequency(list(t1, t2)), "tip label")
})

test_that("NA elements in tree list are silently dropped", {
  t <- TreeTools::RandomTree(8, root = TRUE)
  expect_equal(Loose(list(t, NA, t)), Loose(list(t, t)))
})
