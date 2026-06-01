# Rooted-clade set (sorted descendant tip-label sets, size 2..n-1),
# order- and tip-numbering-independent: uses each tree's own labels.
cladeSet <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree)
  tl <- tree[["tip.label"]]
  clades <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]
    paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sizes <- lengths(strsplit(clades, ","))
  sort(unique(clades[sizes > 1 & sizes < nTip]))
}

# ---------------------------------------------------------------------------
# Edge-case guard: single tree / single-phylo input returns the input
# ---------------------------------------------------------------------------
test_that("Local returns input for single-phylo input", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Local(tree), tree)
})

test_that("Local returns input for single-element list", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Local(list(tree)), tree)
})

test_that("Local rejects non-list input", {
  expect_error(Local(5), "list of trees")
  expect_error(Local("not a tree"), "list of trees")
})

test_that("the local-consensus scaffold self-check is wired up", {
  expect_identical(ConsTree:::consensus_rcpp_selfcheck(), 42L)
})

# ---------------------------------------------------------------------------
# Fewer than 3 leaves: return the input tree
# ---------------------------------------------------------------------------
test_that("Local returns input for fewer than 3 leaves", {
  twoLeaf <- ape::read.tree(text = "(a, b);")
  trees <- list(twoLeaf, twoLeaf)
  result <- Local(trees)
  expect_setequal(result[["tip.label"]], c("a", "b"))
})

# ---------------------------------------------------------------------------
# Result is a valid phylo with all n tips
# ---------------------------------------------------------------------------
test_that("Local returns a valid phylo with the correct tip labels", {
  trees <- lapply(seq_len(5), function(i) ape::as.phylo(i, 8))
  for (tp in c("rooted", "induced")) {
    result <- Local(trees, tp)
    expect_s3_class(result, "phylo")
    expect_setequal(result[["tip.label"]], TreeTools::TipLabels(trees[[1]]))
  }
})

# ---------------------------------------------------------------------------
# Small-n guard: stop above 20 leaves
# ---------------------------------------------------------------------------
test_that("Local errors for n > 20", {
  big <- list(ape::as.phylo(1, 21), ape::as.phylo(2, 21))
  expect_error(Local(big), "20")
})

# ---------------------------------------------------------------------------
# Idempotence: Local(list(t,t,t)) clades == clades of t, both types
# ---------------------------------------------------------------------------
test_that("Local is idempotent (rooted)", {
  tree <- ape::as.phylo(42, 9)
  trees <- list(tree, tree, tree)
  result <- Local(trees, "rooted")
  expect_setequal(cladeSet(result), cladeSet(tree))
})

test_that("Local is idempotent (induced)", {
  tree <- ape::as.phylo(42, 9)
  trees <- list(tree, tree, tree)
  result <- Local(trees, "induced")
  expect_setequal(cladeSet(result), cladeSet(tree))
})

# ---------------------------------------------------------------------------
# Smoke case from the spec:
#   t1 = (1,((2,3),4))  t2 = (1,((2,4),3))
#   MinRLC (rooted) keeps exactly clade {2,3,4} (cherry (2,3) and (2,4)
#   conflict, but both trees agree on the 3-tip clade {2,3,4}).
# ---------------------------------------------------------------------------
test_that("Local smoke case: two conflicting trees keep {2,3,4} only", {
  t1 <- ape::read.tree(text = "(1,((2,3),4));")
  t2 <- ape::read.tree(text = "(1,((2,4),3));")
  trees <- list(t1, t2)
  result_r <- Local(trees, "rooted")
  expect_s3_class(result_r, "phylo")
  expect_setequal(result_r[["tip.label"]], c("1", "2", "3", "4"))
  expect_setequal(cladeSet(result_r), "2,3,4")
})

test_that("Local smoke case: induced also keeps {2,3,4}", {
  t1 <- ape::read.tree(text = "(1,((2,3),4));")
  t2 <- ape::read.tree(text = "(1,((2,4),3));")
  trees <- list(t1, t2)
  result_i <- Local(trees, "induced")
  expect_s3_class(result_i, "phylo")
  expect_setequal(cladeSet(result_i), "2,3,4")
})

# ---------------------------------------------------------------------------
# Star tree: when no common triplets resolve anything, Local returns a star
# ---------------------------------------------------------------------------
test_that("Local returns a star tree for fully conflicting trees", {
  # These two trees have no common triplet:
  # (1,(2,(3,4))) and (2,(1,(3,4))) share {3,4} but as a cherry both agree,
  # so use a more conflicting pair.
  t1 <- ape::read.tree(text = "((1,2),(3,4));")
  t2 <- ape::read.tree(text = "((1,3),(2,4));")
  trees <- list(t1, t2)
  result <- Local(trees, "rooted")
  expect_s3_class(result, "phylo")
  expect_setequal(result[["tip.label"]], c("1", "2", "3", "4"))
  # Star: no internal clade of size 2..n-1
  expect_length(cladeSet(result), 0L)
})

# ---------------------------------------------------------------------------
# Consistent trees: fully resolved input returned unchanged
# ---------------------------------------------------------------------------
test_that("Local on identical fully-resolved trees returns the same clades", {
  tree <- ape::as.phylo(7, 6)
  trees <- list(tree, tree, tree)
  for (tp in c("rooted", "induced")) {
    result <- Local(trees, tp)
    expect_setequal(cladeSet(result), cladeSet(tree))
  }
})

# ---------------------------------------------------------------------------
# Discriminating fixture: rooted and induced genuinely diverge.
# (Confirmed against the oracle binary in dev/oracle/local/check-local.R)
# ---------------------------------------------------------------------------
test_that("Local rooted and induced produce different clades on a conflict case", {
  # Reproducible trees where MinRLC != MinILC (seed=42, trial 52 in random scan)
  set.seed(42)
  found_trees <- NULL
  for (trial in seq_len(200)) {
    n <- sample(5:8, 1)
    k <- sample(3:8, 1)
    trees <- lapply(seq_len(k), function(i) ape::rtree(n, rooted = TRUE))
    labs  <- trees[[1]][["tip.label"]]
    trees <- lapply(trees, function(tr) TreeTools::RenumberTips(tr, labs))
    cr <- tryCatch(cladeSet(Local(trees, "rooted")),  error = function(e) NULL)
    ci <- tryCatch(cladeSet(Local(trees, "induced")), error = function(e) NULL)
    if (!is.null(cr) && !is.null(ci) && !setequal(cr, ci)) {
      found_trees <- trees
      break
    }
  }
  # The scan must find at least one discriminating case
  expect_false(is.null(found_trees),
    info = "Expected to find trees where rooted != induced")
  cr <- cladeSet(Local(found_trees, "rooted"))
  ci <- cladeSet(Local(found_trees, "induced"))
  # The two results must differ
  expect_false(setequal(cr, ci),
    info = "Rooted and induced clades should differ on this input")
})
