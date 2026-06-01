# Owen & Provan (2011) Fig. 1 / Example 2: a geodesic that crosses an
# intermediate orthant (k = 2), for which the paper gives the exact length
# 15 * sqrt(2) and the midpoint tree {3,4}:2.5, {2,3,4,5}:2.5.  This is the
# load-bearing oracle: a broken vertex cover would still pass same-topology
# and cone cases, but not this one.
op_T  <- ape::read.tree(text = "(0:0,(((1:0,2:0):4,(3:0,4:0):10):3,5:0):0);")
op_Tp <- ape::read.tree(text = "(0:0,(1:0,((2:0,3:0):4,(4:0,5:0):3):10):0);")

# Unrooted random tree with positive edge lengths (TreeTools idiom; avoids the
# degree-2 root that `ape::rtree()` introduces).
rtreeBHV <- function(n, labels = letters[seq_len(n)]) {
  tr <- TreeTools::RandomTree(labels, root = FALSE)
  tr[["edge.length"]] <- stats::runif(nrow(tr[["edge"]]), 0.1, 1)
  tr
}

test_that(".TreeToBHV() maps splits to their edge lengths", {
  tree <- ape::read.tree(text = "(a:0,(((b:0,c:0):4,(d:0,e:0):10):3,f:0):0);")
  tl <- sort(tree[["tip.label"]])
  rep <- ConsTree:::.TreeToBHV(tree, tl)
  # three interior splits with lengths 4, 10, 3
  cladeSizes <- rowSums(rep[["membership"]])
  lenBySize <- stats::setNames(rep[["lengths"]], cladeSizes)
  expect_equal(sort(rep[["lengths"]]), c(3, 4, 10))
  expect_equal(unname(lenBySize["4"]), 3)   # the 4-tip clade {b,c,d,e} has length 3
  expect_true(all(rep[["leaf"]] == 0))
})

test_that(".BHVToTree() inverts .TreeToBHV() (round-trip)", {
  set.seed(7)
  tree <- rtreeBHV(9)
  tl <- sort(tree[["tip.label"]])
  rep <- ConsTree:::.TreeToBHV(tree, tl)
  back <- ConsTree:::.BHVToTree(
    list(membership = rep[["membership"]], lengths = rep[["lengths"]],
         leaf = rep[["leaf"]]), tl)
  expect_equal(BHVDistance(tree, back), 0, tolerance = 1e-9)
})

test_that("BHVDistance() reproduces the Owen-Provan 15*sqrt(2) oracle", {
  expect_equal(BHVDistance(op_T, op_Tp), 15 * sqrt(2), tolerance = 1e-9)
})

test_that("BHVDistance() obeys metric invariants", {
  set.seed(1)
  t1 <- rtreeBHV(8)
  t2 <- rtreeBHV(8)
  t3 <- rtreeBHV(8)
  # symmetry
  expect_equal(BHVDistance(t1, t2), BHVDistance(t2, t1))
  # identity of indiscernibles
  expect_equal(BHVDistance(t1, t1), 0)
  # triangle inequality
  expect_lte(BHVDistance(t1, t3), BHVDistance(t1, t2) + BHVDistance(t2, t3) + 1e-9)
  # same topology -> Euclidean distance of branch lengths (one coordinate per
  # edge, so adding 1 to each shifts the point by sqrt(nEdge))
  t1b <- t1
  t1b[["edge.length"]] <- t1[["edge.length"]] + 1
  expect_equal(BHVDistance(t1, t1b), sqrt(length(t1[["edge.length"]])))
  # positive homogeneity: scaling both trees scales the distance
  s1 <- op_T; s1[["edge.length"]] <- op_T[["edge.length"]] * 3.5
  s2 <- op_Tp; s2[["edge.length"]] <- op_Tp[["edge.length"]] * 3.5
  expect_equal(BHVDistance(s1, s2), 3.5 * 15 * sqrt(2), tolerance = 1e-9)
})

test_that("cone path is taken between incompatible single edges", {
  # share clade {1,2,3,4}:5; differ {1,2}:3 (incompatible with) {2,3}:4
  k1 <- ape::read.tree(text = "(0:0,((((1:0,2:0):3,3:0):0,4:0):5,5:0):0);")
  k2 <- ape::read.tree(text = "(0:0,(((1:0,(2:0,3:0):4):0,4:0):5,5:0):0);")
  expect_equal(BHVDistance(k1, k2), 3 + 4)   # cone: ||e|| + ||f||
})

test_that("compatible split changes use the shared orthant, not the cone", {
  # {1,2}:3 (in T1) and {4,5}:4 (in T2) are compatible (disjoint), so the
  # geodesic crosses orthant {12,45}: d = sqrt(3^2 + 4^2) = 5, not the cone 7.
  t1 <- ape::read.tree(text = "((1:0,2:0):3,0:0,3:0,4:0,5:0);")
  t2 <- ape::read.tree(text = "((4:0,5:0):4,0:0,1:0,2:0,3:0);")
  expect_equal(BHVDistance(t1, t2), 5)
})

test_that("geodesic interpolation matches the Owen-Provan midpoint", {
  mid <- ConsTree:::.BHVTreeAt(op_T, op_Tp, 0.5)
  interior <- mid[["edge.length"]][mid[["edge"]][, 2] > length(mid[["tip.label"]])]
  expect_equal(sort(interior), c(2.5, 2.5))
  # midpoint is equidistant, half the geodesic length
  expect_equal(BHVDistance(op_T, mid), 15 * sqrt(2) / 2, tolerance = 1e-9)
  expect_equal(BHVDistance(mid, op_Tp), 15 * sqrt(2) / 2, tolerance = 1e-9)
  # endpoints reproduce the input trees
  expect_equal(BHVDistance(op_T,  ConsTree:::.BHVTreeAt(op_T, op_Tp, 0)), 0,
               tolerance = 1e-9)
  expect_equal(BHVDistance(op_Tp, ConsTree:::.BHVTreeAt(op_T, op_Tp, 1)), 0,
               tolerance = 1e-9)
})

test_that("BHVMean() of two trees converges to their geodesic midpoint", {
  set.seed(1)
  m <- BHVMean(list(op_T, op_Tp))
  mid <- ConsTree:::.BHVTreeAt(op_T, op_Tp, 0.5)
  expect_true(attr(m, "converged"))
  expect_lt(BHVDistance(m, mid), 1e-2)
  expect_equal(BHVDistance(m, op_T), BHVDistance(m, op_Tp), tolerance = 1e-2)
})

test_that("BHVMean() of one orthant is the coordinatewise Euclidean mean", {
  set.seed(2)
  base <- rtreeBHV(7)
  samp <- lapply(1:40, function(i) {
    tr <- base
    tr[["edge.length"]] <- base[["edge.length"]] * stats::runif(1, 0.5, 1.5)
    tr
  })
  m <- BHVMean(samp)
  eu <- base
  eu[["edge.length"]] <- Reduce(`+`, lapply(samp, `[[`, "edge.length")) / length(samp)
  expect_lt(BHVDistance(m, eu), 1e-2)
})

test_that("a single tree is its own mean", {
  expect_equal(BHVDistance(BHVMean(list(op_T)), op_T), 0, tolerance = 1e-9)
})

test_that("BHVMean() minimises the sum of squared geodesic distances", {
  set.seed(8)
  trees <- lapply(1:12, function(i) rtreeBHV(7))
  mn <- BHVMean(trees)
  ss <- function(centre) sum(vapply(trees, function(t) BHVDistance(centre, t)^2,
                                    numeric(1)))
  # the Fréchet mean's sum of squares cannot exceed that of any sample tree
  expect_lte(ss(mn), min(vapply(trees, ss, numeric(1))) + 1e-6)
})

test_that("BHVVariance() relates average and sum correctly", {
  set.seed(3)
  trees <- lapply(1:15, function(i) rtreeBHV(6))
  mn <- BHVMean(trees)
  va <- BHVVariance(trees, mean = mn, type = "average")
  vs <- BHVVariance(trees, mean = mn, type = "sum")
  expect_equal(va, vs / length(trees))
  # explicit definition: average squared geodesic distance from the mean
  d2 <- vapply(trees, function(t) BHVDistance(mn, t)^2, numeric(1))
  expect_equal(va, mean(d2))
})

test_that("BHVVariance() computes its own mean when none is supplied", {
  set.seed(11)
  trees <- lapply(1:8, function(i) rtreeBHV(6))
  v <- BHVVariance(trees)                  # internally calls BHVMean()
  expect_length(v, 1L)
  expect_gte(v, 0)
  # The Fréchet mean is near-optimal, so the self-computed variance closely
  # matches one based on a separately estimated mean (BHVMean() is stochastic).
  expect_equal(v, BHVVariance(trees, mean = BHVMean(trees)), tolerance = 1e-2)
})

test_that("a single-split (four-leaf) tree maps to one BHV coordinate", {
  # A four-leaf unrooted tree has exactly one interior split, the smallest
  # non-trivial .TreeToBHV() case; two such trees differ in that one coordinate.
  q1 <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):1);")
  q2 <- ape::read.tree(text = "((a:1,c:1):3,(b:1,d:1):1);")
  expect_true(is.finite(BHVDistance(q1, q2)))
  expect_equal(BHVDistance(q1, q1), 0)
})

test_that("BHVMean() accepts a bare phylo and an all-star sample", {
  # A bare `phylo` takes the single-tree branch of .BHVTreeList().
  expect_equal(BHVDistance(BHVMean(op_T), op_T), 0, tolerance = 1e-9)
  # A sample of star trees has no interior splits, so the mean is itself a star,
  # exercising the empty-membership branch of .BHVToTree().
  star <- TreeTools::StarTree(letters[1:5])
  star[["edge.length"]] <- rep(1, nrow(star[["edge"]]))
  m <- BHVMean(list(star, star))
  expect_s3_class(m, "phylo")
  expect_equal(TreeTools::NSplits(m), 0)
})

test_that("BHVMean() rejects a collection containing a non-tree", {
  expect_error(BHVMean(list(op_T, "not a tree")),
               "phylo")
})


test_that("BHVDistance() pairwise (list, NULL) returns the same dist as BHVDistance(trees, trees)", {
  set.seed(5)
  trees <- lapply(1:4, function(i) rtreeBHV(7))
  expect_equal(BHVDistance(trees), BHVDistance(trees, trees))
})

test_that("BHVDistance() single-vs-list returns a named vector matching scalar calls", {
  set.seed(6)
  t0 <- rtreeBHV(6)
  trees <- setNames(lapply(1:3, function(i) rtreeBHV(6)), c("A", "B", "C"))
  v <- BHVDistance(t0, trees)
  expect_length(v, 3)
  expect_named(v, c("A", "B", "C"))
  expect_equal(v[["A"]], BHVDistance(t0, trees[["A"]]))
  expect_equal(v[["C"]], BHVDistance(t0, trees[["C"]]))
  # list-vs-single is symmetric
  w <- BHVDistance(trees, t0)
  expect_equal(unname(v), unname(w))
  expect_named(w, c("A", "B", "C"))
})

test_that("BHVDistance() list-vs-list returns a correctly oriented matrix", {
  set.seed(7)
  trees1 <- setNames(lapply(1:2, function(i) rtreeBHV(6)), c("r1", "r2"))
  trees2 <- setNames(lapply(1:3, function(i) rtreeBHV(6)), c("c1", "c2", "c3"))
  m <- BHVDistance(trees1, trees2)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 3L))
  expect_equal(rownames(m), c("r1", "r2"))
  expect_equal(colnames(m), c("c1", "c2", "c3"))
  # entries agree with scalar calls
  expect_equal(m["r1", "c2"], BHVDistance(trees1[["r1"]], trees2[["c2"]]))
  expect_equal(m["r2", "c3"], BHVDistance(trees1[["r2"]], trees2[["c3"]]))
})

test_that("BHVDistance() list-vs-same-list returns dist, not matrix", {
  set.seed(9)
  trees <- lapply(1:4, function(i) rtreeBHV(6))
  d <- BHVDistance(trees, trees)
  expect_s3_class(d, "dist")
  # identical(tree1, tree2) path agrees with the NULL path
  expect_equal(d, BHVDistance(trees))
})

test_that("BHVDistance(single) errors informatively", {
  t1 <- rtreeBHV(5)
  expect_error(BHVDistance(t1), "tree2")
})

test_that("BHVDistance() carries tree names as dist labels", {
  set.seed(10)
  trees <- setNames(lapply(1:3, function(i) rtreeBHV(5)), c("X", "Y", "Z"))
  d <- BHVDistance(trees)
  m <- as.matrix(d)
  expect_equal(rownames(m), c("X", "Y", "Z"))
})

test_that("mismatched leaf labels are rejected", {
  a <- rtreeBHV(6, letters[1:6])
  b <- rtreeBHV(6, letters[2:7])
  expect_error(BHVDistance(a, b), "same")
})

test_that("a zero-length interior edge is treated as absent (no hang)", {
  # {1,2} has length 0 in A and is incompatible with {2,3}:1 in B; a 0-length
  # interior edge is an absent split, so it must not stall the vertex cover.
  a <- ape::read.tree(text = "(0:1,(1:0,2:0):0,3:1,4:1,5:1);")
  b <- ape::read.tree(text = "(0:1,(2:1,3:1):1,1:1,4:1,5:1);")
  d <- BHVDistance(a, b)
  expect_true(is.finite(d))
  # equals the limit as the 0-length edge shrinks: same as dropping it
  aDrop <- ape::read.tree(text = "(0:1,1:0,2:0,3:1,4:1,5:1);")
  expect_equal(d, BHVDistance(aDrop, b))
})

test_that("internal singleton nodes do not change the tree's position", {
  # a tree and its collapse.singles() equivalent are identical in BHV space
  b <- ape::read.tree(
    text = "(0:1,((1:1,(2:1,(3:1,(4:1,5:1):1):1):1):2):1,6:1,7:1);")
  expect_true(any(table(b[["edge"]][, 1]) == 1))   # has a singleton
  expect_equal(BHVDistance(b, ape::collapse.singles(b)), 0, tolerance = 1e-12)
})

test_that("trees without edge lengths are rejected", {
  a <- TreeTools::RandomTree(6, root = FALSE)            # no edge.length
  b <- rtreeBHV(6, a[["tip.label"]])
  expect_error(BHVDistance(a, b), "edge.length")
})
