# Two trees share an unrooted topology iff they induce the same set of splits.
# PolarizeSplits() normalises each split so that tip 1 is always on the
# zero-side, giving a canonical form that makes %in% on Splits polarity-safe.
# Without it, as.Splits() can encode the same bipartition with opposite polarity
# depending on the tree's internal rooting, causing false negatives.
unrootedMatch <- function(a, b) {
  labels <- TreeTools::TipLabels(a)
  sA <- TreeTools::PolarizeSplits(TreeTools::as.Splits(a, tipLabels = labels))
  sB <- TreeTools::PolarizeSplits(TreeTools::as.Splits(b, tipLabels = labels))
  length(sA) == length(sB) && all(sA %in% sB)
}

# Residual sum of squares of a tree's path-length distances against a target,
# computed independently of the package's internal fitter.
pairRss <- function(tree, target, labs) {
  d <- ape::cophenetic.phylo(ape::unroot(tree))[labs, labs]
  sum((d[lower.tri(d)] - target[lower.tri(target)]) ^ 2)
}

# `method = "ls"` needs a TreeSearch build that provides LeastSquaresTree().
skip_without_ls <- function() {
  skip_if_not_installed("TreeSearch")
  skip_if(!exists("LeastSquaresTree", where = asNamespace("TreeSearch"),
                  mode = "function"),
          "Installed TreeSearch lacks LeastSquaresTree()")
}

test_that("Average() returns the least-squares tree on the toy example", {
  skip_without_ls()
  # ((((A,X),B),C),D) + (((A,B),C),(D,X)): averaging the path-length matrices
  # and fitting the closest tree returns input tree T2 (verified by exhaustive
  # least-squares over all 15 unrooted five-leaf topologies), not the apparent
  # compromise ((((A,B),X),C),D).
  t1 <- ape::read.tree(text = "((((A,X),B),C),D);")
  t2 <- ape::read.tree(text = "(((A,B),C),(D,X));")

  expect_true(unrootedMatch(Average(list(t1, t2), method = "ls"), t2))
  expect_true(unrootedMatch(Average(list(t1, t2), method = "nj"), t2))
  expect_true(unrootedMatch(Average(list(t1, t2), method = "fastme.bal"), t2))
})

test_that("Average() fits an additive matrix exactly", {
  skip_without_ls()
  set.seed(1)
  tree <- ape::rtree(8)
  # The average of identical matrices is itself additive, so least squares
  # should recover the generating tree with zero residual.
  avg <- Average(list(tree, tree), method = "ls")
  expect_true(unrootedMatch(avg, tree))

  labs <- tree[["tip.label"]]
  target <- ape::cophenetic.phylo(tree)[labs, labs]
  fitted <- ape::cophenetic.phylo(avg)[labs, labs]
  expect_equal(fitted, target, tolerance = 1e-6)
})

test_that("Average('ls') finds the exhaustive optimum on real branch lengths", {
  skip_without_ls()
  # Two distinct topologies with real branch lengths: the average path-length
  # matrix is non-additive (the genuine posterior-summary case).
  set.seed(1)
  trees <- list(ape::rtree(6), ape::rtree(6))
  labs <- trees[[1]][["tip.label"]]
  target <- Reduce(`+`, lapply(trees, function(tr) {
    ape::cophenetic.phylo(tr)[labs, labs]
  })) / length(trees)

  # Exhaustive least squares: every binary topology on six leaves (enumerated as
  # rooted trees, which redundantly cover all unrooted topologies).
  candidates <- ape::as.phylo(seq_len(TreeTools::NRooted(6)) - 1L, 6L,
                              tipLabels = labs)
  exhaustive <- min(vapply(candidates, function(tr) {
    pairRss(TreeSearch::LeastSquaresFit(tr, target, method = "nnls"), target, labs)
  }, numeric(1)))

  expect_equal(pairRss(Average(trees, method = "ls"), target, labs),
               exhaustive, tolerance = 1e-8)
})

test_that("A single tree is its own average", {
  tree <- ape::rtree(6)               # rooted, with branch lengths
  expect_true(unrootedMatch(Average(list(tree)), tree))
  expect_true(unrootedMatch(Average(tree), tree))
  # honour the unrooted-by-default contract even on the single-tree path
  expect_false(ape::is.rooted(Average(list(tree))))
  expect_true(ape::is.rooted(Average(tree, outgroup = tree[["tip.label"]][[1]])))
})

test_that("Average() honours distance methods, scaling and rooting", {
  trees <- ape::rmtree(4, 7)
  for (m in c("nj", "bionj", "fastme.bal", "fastme.ols")) {
    tr <- Average(trees, method = m)
    expect_s3_class(tr, "phylo")
    expect_setequal(tr[["tip.label"]], trees[[1]][["tip.label"]])
  }

  expect_s3_class(Average(trees, method = "nj", scale = "max"), "phylo")
  expect_s3_class(Average(trees, method = "nj", edgeLengths = FALSE), "phylo")

  rooted <- Average(trees, method = "nj", outgroup = trees[[1]][["tip.label"]][[1]])
  expect_true(ape::is.rooted(rooted))
})

test_that("scale = 'max' makes the average invariant to per-tree rescaling", {
  set.seed(42)
  a <- ape::rtree(7)
  b <- ape::rtree(7)
  b10 <- b
  b10[["edge.length"]] <- b10[["edge.length"]] * 10

  expect_true(unrootedMatch(Average(list(a, b), scale = "max", method = "nj"),
                            Average(list(a, b10), scale = "max", method = "nj")))
  expect_false(unrootedMatch(Average(list(a, b), scale = "none", method = "nj"),
                             Average(list(a, b10), scale = "none", method = "nj")))
})

test_that("Average() validates its input", {
  good <- ape::read.tree(text = "((A,B),(C,D));")
  bad  <- ape::read.tree(text = "((A,B),(C,E));")
  expect_error(Average(list(good, bad)), "same leaf labels")
  expect_error(Average("not a tree"), "list of trees")

  trees <- ape::rmtree(2, 6)
  expect_error(Average(trees, weights = c(1, 2, 3)), "one entry per tree")
  expect_error(Average(trees, weights = c(-1, 1)), "non-negative")
  expect_s3_class(Average(trees, edgeLengths = TRUE, method = "nj"), "phylo")
})

test_that("weights shift the average toward the favoured tree", {
  t1 <- ape::read.tree(text = "((((A,X),B),C),D);")
  t2 <- ape::read.tree(text = "(((A,B),C),(D,X));")
  # all weight on one tree -> that tree's (additive) matrix -> that tree
  expect_true(unrootedMatch(Average(list(t1, t2), method = "nj",
                                    weights = c(1, 0)), t1))
  expect_true(unrootedMatch(Average(list(t1, t2), method = "nj",
                                    weights = c(0, 1)), t2))
})

test_that("edgeLengths auto-detection counts edges when lengths are mixed", {
  withLen <- ape::rtree(6)
  noLen <- withLen
  noLen[["edge.length"]] <- NULL
  # one tree lacks branch lengths, so the default (NA) falls back to edge counts
  expect_s3_class(Average(list(withLen, noLen), method = "nj"), "phylo")
  expect_error(Average(list(withLen, noLen), edgeLengths = TRUE),
               "not every tree has branch lengths")
})
