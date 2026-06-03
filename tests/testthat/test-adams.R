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

# Independent, deliberately naive reference for the classical Adams recursion --
# the definition the C++ must reproduce.  At each node it partitions the taxa by
# the cross-tree "which child of the LCA" signature (a direct restate of Adams
# 1972), and collects every non-trivial block as a clade.  Slow (KeepTip per
# node), so used only to fuzz Adams() at small n without invoking fact.exe.
refAdamsCladeSet <- function(trees) {
  labels <- TreeTools::TipLabels(trees[[1]])
  treesP <- lapply(trees, function(tr)
    TreeTools::Preorder(TreeTools::RenumberTips(tr, labels)))
  partition <- function(taxa) {
    sig <- vapply(treesP, function(tr) {
      kept <- TreeTools::Preorder(TreeTools::KeepTip(tr, taxa))
      edge <- kept[["edge"]]
      root <- edge[1L, 1L]
      nt <- TreeTools::NTip(kept)
      bid <- integer(nt)
      for (b in which(edge[, 1L] == root)) {
        below <- TreeTools::DescendantEdges(edge[, 1L], edge[, 2L], edge = b)
        ct <- edge[below, 2L]
        bid[ct[ct <= nt]] <- b
      }
      bid[match(taxa, kept[["tip.label"]])]
    }, integer(length(taxa)))
    if (is.null(dim(sig))) sig <- matrix(sig, nrow = length(taxa))
    unname(split(taxa, apply(sig, 1L, paste, collapse = ",")))
  }
  clades <- character(0)
  rec <- function(taxa) {
    if (length(taxa) <= 2L) return(invisible())
    for (B in partition(taxa)) {
      if (length(B) >= 2L) clades[[length(clades) + 1L]] <<- paste(sort(B), collapse = ",")
      rec(B)
    }
  }
  rec(labels)
  sort(unique(clades))
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

test_that("Adams matches the classical reference across random profiles", {
  cfgs <- list(c(n = 5L, k = 2L), c(n = 6L, k = 3L), c(n = 8L, k = 4L),
               c(n = 10L, k = 3L), c(n = 12L, k = 5L), c(n = 9L, k = 8L))
  for (cfg in cfgs) {
    for (seed in 1:6) {
      set.seed(seed * 1000L + cfg[["n"]])
      labs <- paste0("t", seq_len(cfg[["n"]]))
      trees <- structure(
        lapply(seq_len(cfg[["k"]]),
               function(i) TreeTools::RandomTree(labs, root = TRUE)),
        class = "multiPhylo")
      expect_setequal(cladeSet(Adams(trees)), refAdamsCladeSet(trees))
    }
  }
})

test_that("Adams matches the reference on caterpillars and shared structure", {
  catTree <- function(labs) {
    nwk <- labs[length(labs)]
    for (i in (length(labs) - 1L):1L) nwk <- paste0("(", labs[[i]], ",", nwk, ")")
    ape::read.tree(text = paste0(nwk, ";"))
  }
  labs <- paste0("t", seq_len(15L))
  set.seed(99L)
  # distinct caterpillars (long spines: the centroid-path walk's stress case)
  cats <- structure(lapply(1:4, function(i) catTree(sample(labs))),
                    class = "multiPhylo")
  expect_setequal(cladeSet(Adams(cats)), refAdamsCladeSet(cats))
  # a shared backbone perturbed by tip swaps (congruent: many real clades)
  base <- TreeTools::RandomTree(labs, root = TRUE)
  pert <- structure(lapply(1:5, function(i) {
    tr <- base
    for (s in 1:2) {
      ij <- sample.int(15L, 2L)
      tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
    }
    tr
  }), class = "multiPhylo")
  expect_setequal(cladeSet(Adams(pert)), refAdamsCladeSet(pert))
})

test_that("Adams emits no degree-1 nodes (spine suppression)", {
  # A star input: the spine collapses to one leaf; the result must be a clean
  # star, never a leaf wrapped in a degree-1 internal node.
  star <- c(ape::read.tree(text = "((1, 2), (3, 4));"),
            ape::read.tree(text = "((1, 3), (2, 4));"))
  ad <- Adams(star)
  expect_false(any(tabulate(ad[["edge"]][, 1L]) == 1L))
  expect_setequal(ad[["tip.label"]], c("1", "2", "3", "4"))
  expect_length(cladeSet(ad), 0L)
  # Identical caterpillars must reproduce the caterpillar exactly.
  cat <- ape::read.tree(text = "((((a, b), c), d), e);")
  trees <- structure(list(cat, cat, cat), class = "multiPhylo")
  ad2 <- Adams(trees)
  expect_false(any(tabulate(ad2[["edge"]][, 1L]) == 1L))
  expect_setequal(cladeSet(ad2), cladeSet(cat))
})

test_that("Adams resolves three leaves", {
  agree <- c(ape::read.tree(text = "((a, b), c);"),
             ape::read.tree(text = "((a, b), c);"))
  expect_setequal(cladeSet(Adams(agree)), "a,b")
  disagree <- c(ape::read.tree(text = "((a, b), c);"),
                ape::read.tree(text = "((a, c), b);"))
  expect_length(cladeSet(Adams(disagree)), 0L)
})

test_that("Adams suppresses a degree-1 spine node (deepest single-block step)", {
  # A rare n=7/k=2 configuration whose deepest spine step resolves to a single
  # block with no spine-bottom leaf remaining: this exercises the degree-1
  # suppression branch of the chain assembly (pass the lone child through rather
  # than wrap it in a one-child internal node).
  trees <- structure(list(
    ape::read.tree(text = "((t1, t6), (((t2, t5), (t3, t4)), t7));"),
    ape::read.tree(text = "((((t1, (t6, t7)), (t3, t4)), t5), t2);")
  ), class = "multiPhylo")
  ad <- Adams(trees)
  expect_false(any(tabulate(ad[["edge"]][, 1L]) == 1L))
  expect_setequal(cladeSet(ad), refAdamsCladeSet(trees))
})

test_that("Adams requires a common leaf set", {
  t1 <- ape::read.tree(text = "((a, b), c);")
  t2 <- ape::read.tree(text = "((a, b), d);")  # 'd' replaces 'c'
  expect_error(Adams(structure(list(t1, t2), class = "multiPhylo")),
               "same leaves")
})
