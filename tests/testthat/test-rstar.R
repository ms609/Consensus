# Tests for RStar() (R* consensus).
#
# R* has no reference binary, but its definition is exact (Jansson, Sung, Vu &
# Yiu 2016, Lemma 1.1): the R* tree's clades are exactly the STRONG CLUSTERS of
# R_maj.  So beyond the worked examples and guards we assert two provable facts:
#   * the brute-force strong-cluster oracle (RStar clades == strong clusters,
#     computed independently via ape::mrca over all 2^n subsets for small n); and
#   * majority-rule refinement (every rooted majority clade appears in RStar) --
#     which Lemma 1.1 guarantees and which the earlier BUILD implementation
#     violated.

# Rooted clades (size 2..n-1) as sorted, comma-joined tip-label strings.
cladeSet <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree)
  tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]
    paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sizes <- lengths(strsplit(cs, ","))
  sort(unique(cs[sizes > 1 & sizes < nTip]))
}
tt <- function(x) ape::read.tree(text = x)

# --- Independent strong-cluster machinery (ape::mrca; not the C++ path) -------
.tripletClose <- function(tree) {
  tree2 <- tree; tree2[["edge.length"]] <- rep(1, nrow(tree[["edge"]]))
  dep <- ape::node.depth.edgelength(tree2)
  m <- ape::mrca(tree)
  function(x, y, z) {
    dxy <- dep[m[x, y]]; dxz <- dep[m[x, z]]; dyz <- dep[m[y, z]]
    if (dxy == dxz && dxy == dyz) return(NA_character_)            # fan -> abstain
    if (dxy > dxz && dxy > dyz) return(paste(sort(c(x, y)), collapse = ","))
    if (dxz > dxy && dxz > dyz) return(paste(sort(c(x, z)), collapse = ","))
    paste(sort(c(y, z)), collapse = ",")
  }
}
.favoured <- function(staters, x, y, Z) {     # is (xy)|Z uniquely favoured?
  v <- vapply(staters, function(st) st(x, y, Z), character(1))
  v <- v[!is.na(v)]
  if (!length(v)) return(FALSE)
  tb <- table(v); tgt <- paste(sort(c(x, y)), collapse = ",")
  if (is.na(tb[tgt])) return(FALSE)
  mx <- max(tb); tb[tgt] == mx && sum(tb == mx) == 1L
}
.isStrong <- function(staters, A, allTips) {   # for-all-outgroup rule
  out <- setdiff(allTips, A)
  if (!length(out)) return(TRUE)
  pr <- utils::combn(A, 2)
  for (p in seq_len(ncol(pr))) for (Z in out)
    if (!.favoured(staters, pr[1, p], pr[2, p], Z)) return(FALSE)
  TRUE
}
.strongClusters <- function(trees) {            # brute force over 2^n subsets
  allTips <- trees[[1]][["tip.label"]]; n <- length(allTips)
  st <- lapply(trees, .tripletClose)
  out <- character(0)
  for (mask in seq_len(2L^n - 1L)) {
    idx <- which(as.integer(intToBits(mask))[seq_len(n)] == 1L)
    if (length(idx) < 2L || length(idx) > n - 1L) next
    A <- allTips[idx]
    if (.isStrong(st, A, allTips)) out <- c(out, paste(sort(A), collapse = ","))
  }
  sort(unique(out))
}
.alignTrees <- function(trees) {
  labs <- trees[[1]][["tip.label"]]
  lapply(trees, function(tr) TreeTools::RenumberTips(tr, labs))
}
.majClades <- function(trees) {
  k <- length(trees)
  tab <- table(unlist(lapply(trees, cladeSet)))
  names(tab)[tab > k / 2]
}

test_that("RStar() recovers the worked example by plurality", {
  # Degnan et al. worked example: (a,b) wins {a,b,c} 3-1-1 -> (((a,b),c),d)
  trees <- c(
    tt("(((a,b),c),d);"), tt("(((a,b),c),d);"), tt("(((a,b),c),d);"),
    tt("(((a,c),b),d);"), tt("(((b,c),a),d);")
  )
  expect_setequal(cladeSet(RStar(trees)), c("a,b", "a,b,c"))
})

test_that("RStar() uses strict plurality, not majority", {
  # 2-1-1 split: (a,b) wins {a,b,c} with only 2 of 4 votes (below strict majority)
  trees <- c(
    tt("(((a,b),c),d);"), tt("(((a,b),c),d);"),
    tt("(((a,c),b),d);"), tt("(((b,c),a),d);")
  )
  cs <- cladeSet(RStar(trees))
  expect_true("a,b" %in% cs)        # plurality keeps it
  expect_true("a,b,c" %in% cs)
})

test_that("RStar() leaves tied triples unresolved", {
  # 2-2-0 tie on {a,b,c} -> polytomy ((a,b,c),d); no resolved cherry
  trees <- c(
    tt("(((a,b),c),d);"), tt("(((a,b),c),d);"),
    tt("(((a,c),b),d);"), tt("(((a,c),b),d);")
  )
  cs <- cladeSet(RStar(trees))
  expect_setequal(cs, "a,b,c")
  expect_false("a,b" %in% cs)
  expect_false("a,c" %in% cs)
})

test_that("RStar() is idempotent on identical trees", {
  one <- tt("((((a,b),c),(d,e)),f);")
  expect_setequal(cladeSet(RStar(c(one, one, one))), cladeSet(one))
})

test_that("RStar() returns a single input unchanged", {
  one <- tt("(((a,b),c),d);")
  expect_setequal(cladeSet(RStar(list(one))), cladeSet(one))   # k = 1
  expect_identical(RStar(one), one)                            # bare phylo
})

test_that("RStar() retains clades shared by every tree (strict refinement)", {
  # {a,b} is a clade in all three trees -> must appear in R*
  trees <- c(
    tt("(((a,b),c),(d,e));"),
    tt("(((a,b),d),(c,e));"),
    tt("(((a,b),e),(c,d));")
  )
  expect_true("a,b" %in% cladeSet(RStar(trees)))
})

test_that("RStar() clades are exactly the strong clusters of R_maj", {
  # Brute-force oracle (Lemma 1.1) on random conflicting input, small n.
  set.seed(2024)
  for (trial in seq_len(12)) {
    n <- sample(4:7, 1); k <- sample(2:7, 1)
    trees <- .alignTrees(lapply(seq_len(k), function(i) {
      tr <- ape::rtree(n, rooted = TRUE); tr[["edge.length"]] <- NULL; tr
    }))
    expect_setequal(cladeSet(RStar(trees)), .strongClusters(trees))
  }
})

test_that("RStar() handles non-binary (partly resolved) input", {
  # Fans must abstain (have no impact on R_maj); the strong-cluster oracle must
  # still hold when inputs contain a mix of fans and resolved triples.
  set.seed(7)
  for (trial in seq_len(10)) {
    n <- sample(5:7, 1); k <- sample(2:6, 1)
    trees <- .alignTrees(lapply(seq_len(k), function(i) {
      tr <- ape::rtree(n, rooted = TRUE)
      ape::di2multi(tr, tol = stats::runif(1, 0.15, 0.45))  # collapse short edges
    }))
    expect_setequal(cladeSet(RStar(trees)), .strongClusters(trees))
  }
})

test_that("RStar() refines the majority-rule consensus", {
  # Lemma 1.1: every rooted majority clade is a strong cluster, hence in R*.
  set.seed(99)
  for (trial in seq_len(20)) {
    n <- sample(5:10, 1); k <- sample(3:9, 1)
    trees <- .alignTrees(lapply(seq_len(k), function(i) {
      tr <- ape::rtree(n, rooted = TRUE); tr[["edge.length"]] <- NULL; tr
    }))
    expect_true(all(.majClades(trees) %in% cladeSet(RStar(trees))))
  }
})

test_that("RStar() returns a star for star input", {
  st <- tt("(a,b,c,d,e);")
  expect_length(cladeSet(RStar(c(st, st))), 0L)
})

test_that("RStar() scales past Local's 20-leaf limit", {
  set.seed(7)
  big <- ape::rtree(30, rooted = TRUE)
  big$edge.length <- NULL
  expect_setequal(cladeSet(RStar(c(big, big))), cladeSet(big))
})

test_that("RStar() preserves the leaf set and returns a rooted phylo", {
  trees <- c(tt("(((a,b),c),(d,e));"), tt("(((a,c),b),(d,e));"))
  res <- RStar(trees)
  expect_s3_class(res, "phylo")
  expect_true(ape::is.rooted(res))
  expect_setequal(res[["tip.label"]], c("a", "b", "c", "d", "e"))
})

test_that("RStar() runs past the former 200-leaf cap", {
  # The dense n^3 triplet tensor (and its hard 200-leaf cap) is gone: memory is
  # now O(k n^2), so large leaf counts run.  Identity past the old cap.
  set.seed(1)
  big <- ape::rtree(260, rooted = TRUE)
  big$edge.length <- NULL
  expect_setequal(cladeSet(RStar(c(big, big))), cladeSet(big))
  # The C++ core no longer caps at 200 leaves either (returns Newick, not error).
  expect_type(ConsTree:::rStarConsensus(list(), 260L), "character")
})

test_that("RStar() refines strict & majority clades at larger n (cap lifted)", {
  # Lemma 1.1 refinement at sizes the former cap forbade and the brute-force
  # oracle (n <= 12) cannot reach.
  set.seed(42)
  for (trial in seq_len(6)) {
    n <- sample(40:80, 1); k <- sample(3:7, 1)
    trees <- .alignTrees(lapply(seq_len(k), function(i) {
      tr <- ape::rtree(n, rooted = TRUE); tr[["edge.length"]] <- NULL; tr
    }))
    rs <- cladeSet(RStar(trees))
    tab <- table(unlist(lapply(trees, cladeSet)))
    expect_true(all(names(tab)[tab == k] %in% rs))       # strict refinement
    expect_true(all(names(tab)[tab > k / 2] %in% rs))    # majority refinement
  }
})

test_that("RStar() rejects non-list input and trivial leaf sets", {
  expect_error(RStar(5), "list of trees")
  expect_error(RStar("not a tree"), "list of trees")
  # Fewer than three leaves: the first tree is returned unchanged.
  twoLeaf <- tt("(a, b);")
  expect_identical(RStar(c(twoLeaf, twoLeaf)), twoLeaf)
})
