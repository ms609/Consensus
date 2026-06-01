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
  trees <- ape::as.phylo(0:20, 9)
  labels <- TreeTools::TipLabels(trees[[1]])
  s <- splitSet(Strict(trees), labels)
  l <- splitSet(Loose(trees), labels)
  m <- splitSet(Majority(trees), labels)
  g <- splitSet(Greedy(trees), labels)
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
  for (spec in list(0:20, 0:30, c(0, 0, 0, 1, 2, 53, 99))) {
    trees <- ape::as.phylo(spec, 9)
    labels <- TreeTools::TipLabels(trees[[1]])
    m <- splitSet(Majority(trees), labels)
    f <- splitSet(Frequency(trees), labels)
    g <- splitSet(Greedy(trees), labels)
    expect_true(all(m %in% f))   # majority        <= frequency-difference
    expect_true(all(f %in% g))   # frequency-diff   <= greedy
    expect_s3_class(Frequency(trees), "phylo")
  }
})

test_that("A single tree, or fewer than four leaves, returns the input", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Loose(list(tree)), tree)
  expect_equal(Greedy(list(tree)), tree)
  triplet <- ape::as.phylo(0:1, 3)
  expect_equal(Loose(triplet), triplet[[1]])
})

test_that("Greedy matches phangorn's greedy consensus (allCompat)", {
  # `phangorn::allCompat()` triggers an R-devel/phangorn ABI mismatch that is
  # undefined behaviour: it sometimes throws a catchable error and sometimes
  # segfaults (which `tryCatch()` cannot intercept), so it must not run in the
  # default suite.  Greedy is already validated exactly against the FACT
  # reference binary (dev/oracle), so this is an optional redundant cross-check;
  # enable it in a healthy environment with CONSENSUS_PHANGORN_TESTS=1.
  skip_if_not(identical(Sys.getenv("CONSENSUS_PHANGORN_TESTS"), "1"),
              "phangorn allCompat cross-check disabled (R-devel/phangorn ABI)")
  skip_if_not_installed("phangorn")
  trees <- ape::as.phylo(0:30, 10)
  labels <- TreeTools::TipLabels(trees[[1]])
  expect_setequal(splitSet(Greedy(trees), labels),
                  splitSet(phangorn::allCompat(trees), labels))
})
