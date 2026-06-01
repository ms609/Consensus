# Rooted-clade set (sorted descendant tip-label sets, size 2..n-1), order- and
# tip-numbering-independent: uses each tree's own labels.
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

test_that("Adams keeps clades shared by every tree, drops conflicting ones", {
  # Two rooted trees that agree on the clade (a, b) but disagree on c vs d
  trees <- c(ape::read.tree(text = "(((a, b), c), d);"),
             ape::read.tree(text = "(((a, b), d), c);"))
  ad <- Adams(trees)
  expect_s3_class(ad, "phylo")
  expect_setequal(ad[["tip.label"]], c("a", "b", "c", "d"))
  # Root resolves to ((a, b), c, d): only the shared clade (a, b) survives
  expect_setequal(cladeSet(ad), "a,b")
})

test_that("Adams returns a star when the trees share no nesting", {
  trees <- c(ape::read.tree(text = "((1, 2), (3, 4));"),
             ape::read.tree(text = "((1, 3), (2, 4));"))
  expect_length(cladeSet(Adams(trees)), 0L)
})

test_that("Adams can recover a clade present in no input tree", {
  # T1 nests (c, d); T2 nests (a, b); neither contains (b, c)
  t1 <- ape::read.tree(text = "(a, (b, (c, d)));")
  t2 <- ape::read.tree(text = "(d, (c, (a, b)));")
  trees <- structure(list(t1, t2), class = "multiPhylo")
  novel <- "b,c"
  expect_false(novel %in% cladeSet(t1))
  expect_false(novel %in% cladeSet(t2))
  expect_setequal(cladeSet(Adams(trees)), novel)
})

test_that("Adams is idempotent", {
  tree <- ape::as.phylo(42, 9)
  trees <- structure(list(tree, tree, tree), class = "multiPhylo")
  expect_setequal(cladeSet(Adams(trees)), cladeSet(tree))
})

test_that("Adams of a single tree, or fewer than three leaves, returns the input", {
  tree <- ape::as.phylo(1, 9)
  expect_equal(Adams(tree), tree)
  expect_equal(Adams(list(tree)), tree)
  twoLeaf <- ape::read.tree(text = "(a, b);")
  expect_equal(Adams(structure(list(twoLeaf, twoLeaf), class = "multiPhylo")),
               twoLeaf)
})

test_that("Adams rejects non-list input", {
  expect_error(Adams(5), "list of trees")
  expect_error(Adams("not a tree"), "list of trees")
})

test_that("Adams preserves the full leaf set", {
  trees <- ape::as.phylo(0:9, 12)
  ad <- Adams(trees)
  expect_s3_class(ad, "phylo")
  expect_setequal(ad[["tip.label"]], TreeTools::TipLabels(trees[[1]]))
})
