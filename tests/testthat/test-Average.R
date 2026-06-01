# Two trees share an unrooted topology iff they induce the same set of splits;
# a TreeTools-native stand-in for `phangorn::RF.dist(...) == 0` that keeps this
# helper off the R-devel/phangorn ABI path.  `%in%` on `Splits` is
# polarity-invariant -- it treats A|B and B|A as the same split -- so trees
# built by different methods (with different tip orderings) still compare equal;
# `as.character()` is *not* polarity-canonical, so a naive split-string set
# comparison spuriously reports identical topologies as different.
unrootedMatch <- function(a, b) {
  labels <- TreeTools::TipLabels(a)
  sA <- TreeTools::as.Splits(a, tipLabels = labels)
  sB <- TreeTools::as.Splits(b, tipLabels = labels)
  length(sA) == length(sB) && all(sA %in% sB)
}

test_that("Average() returns the least-squares tree on the toy example", {
  skip_if_not_installed("TreeSearch")
  # ((((A,X),B),C),D) + (((A,B),C),(D,X)): averaging the path-length matrices
  # and fitting the closest tree returns input tree T2 (verified by exhaustive
  # least-squares over all 15 unrooted five-leaf topologies), not the apparent
  # compromise ((((A,B),X),C),D).
  t1 <- ape::read.tree(text = "((((A,X),B),C),D);")
  t2 <- ape::read.tree(text = "(((A,B),C),(D,X));")

  expect_true(unrootedMatch(Average(list(t1, t2), method = "ls"), t2))
  expect_true(unrootedMatch(Average(list(t1, t2), method = "nj"), t2))
  expect_true(unrootedMatch(Average(list(t1, t2), method = "fastme.bal"), t2))

  # the rearrangement search escapes a deliberately poor starting tree
  expect_true(unrootedMatch(Average(list(t1, t2), method = "ls", start = t1), t2))
})

test_that("Average() fits an additive matrix exactly", {
  skip_if_not_installed("TreeSearch")
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
  skip_if_not_installed("TreeSearch")
  # Two distinct topologies with real branch lengths: the average path-length
  # matrix is non-additive (the genuine posterior-summary case).
  set.seed(1)
  trees <- list(ape::rtree(6), ape::rtree(6))
  labs <- trees[[1]][["tip.label"]]
  target <- Reduce(`+`, lapply(trees, function(tr) {
    ape::cophenetic.phylo(tr)[labs, labs]
  })) / length(trees)

  rss <- function(tr) {
    fit <- phangorn::nnls.tree(stats::as.dist(target), ape::unroot(tr),
                               method = "unrooted", trace = 0)
    d <- ape::cophenetic.phylo(fit)[labs, labs]
    sum((d[lower.tri(d)] - target[lower.tri(target)]) ^ 2)
  }
  exhaustive <- min(vapply(
    phangorn::allTrees(6, rooted = FALSE, tip.label = labs), rss, numeric(1)))

  expect_equal(rss(Average(trees, method = "ls")), exhaustive, tolerance = 1e-8)
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

test_that("A single tree is its own average", {
  tree <- ape::rtree(6)
  expect_true(unrootedMatch(Average(list(tree)), tree))
  expect_true(unrootedMatch(Average(tree), tree))
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

test_that("Average() validates its input", {
  good <- ape::read.tree(text = "((A,B),(C,D));")
  bad  <- ape::read.tree(text = "((A,B),(C,E));")
  expect_error(Average(list(good, bad)), "same leaf labels")
  expect_error(Average("not a tree"), "list of trees")

  trees <- ape::rmtree(2, 6)
  expect_error(Average(trees, weights = c(1, 2, 3)), "one entry per tree")
  expect_error(Average(trees, weights = c(-1, 1)), "non-negative")
  # trees from rmtree() have branch lengths, so edgeLengths = TRUE is satisfiable
  expect_s3_class(Average(trees, edgeLengths = TRUE, method = "nj"), "phylo")
})
