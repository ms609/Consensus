# Verify Transfer() fixes against an INSTALLED build (load_all is broken for
# this path — see TC-008). Run after `R CMD INSTALL --library=<tmplib>`:
#   Rscript verify-fixes.R
.libPaths(c("C:/Users/pjjg18/cons-tmplib", .libPaths()))
suppressMessages(library(TreeTools))
suppressMessages(library(ConsTree))
ok <- function(s) { cat(s, "\n"); flush.console() }

# Canonical split-set key: members of the side that excludes tip 1.
splitset <- function(tr, tipLabels) {
  sp <- as.logical(as.Splits(tr, tipLabels))
  if (is.null(dim(sp))) sp <- matrix(sp, nrow = 1L)
  sort(apply(sp, 1L, function(r) {
    side <- if (r[1L]) which(!r) else which(r)
    paste(sort(side), collapse = ",")
  }))
}

# 1. Sanity: homogeneous input still works
trees <- as.phylo(0:9, nTip = 8)
ok(paste("sanity NSplits:", NSplits(Transfer(trees))))

# 2. TC-002: mismatched tip SETS (subset) must now ERROR
m_mismatch <- structure(list(BalancedTree(letters[1:6]),
                             BalancedTree(letters[1:5])), class = "multiPhylo")
ok(paste("TC-002 mismatch:",
         tryCatch({ Transfer(m_mismatch); "NO ERROR (BUG)" },
                  error = function(e) paste("ERROR ok:", conditionMessage(e)))))

# 2b. Same SET, different ORDER must still WORK (setequal, not identical)
t8a <- BalancedTree(letters[1:8])
t8b <- RenumberTips(t8a, rev(letters[1:8]))
m_order <- structure(list(t8a, t8b), class = "multiPhylo")
ok(paste("TC-002 reorder:",
         tryCatch(paste("OK NSplits", NSplits(Transfer(m_order))),
                  error = function(e) paste("ERROR (BAD):", conditionMessage(e)))))

# 3. TC-005/006: C++ greedy path vs pure-R reimplementation, split-set identity
tipLabels <- TipLabels(trees[[1]]); nTip <- length(tipLabels); nTree <- length(trees)
tc_cpp <- Transfer(trees, greedy = "best", scale = TRUE)
pool <- ConsTree:::.PoolSplits(trees, tipLabels)
nSplits <- nrow(pool$splits)
DIST <- ConsTree:::.TransferDistMat(pool$splits, nTip)
sentDist <- pool$lightSide - 1L
TD <- ConsTree:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                            nTree, scale = TRUE)
compat <- ConsTree:::.CompatMat(pool$splits, nTip)
sortOrd <- order(-pool$counts, seq_len(nSplits))
st <- new.env(parent = emptyenv())
st$MATCH <- rep(NA_integer_, nSplits); st$MATCH2 <- rep(NA_integer_, nSplits)
st$incl <- rep(FALSE, nSplits)
ConsTree:::.GreedyBest(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                       compat, sortOrd, scale = TRUE, nSplits, nTip)
tc_r <- ConsTree:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)
s_cpp <- splitset(tc_cpp, tipLabels); s_r <- splitset(tc_r, tipLabels)
ok(paste("TC-005 NSplits cpp/r:", NSplits(tc_cpp), "/", NSplits(tc_r)))
ok(paste("TC-005 split-sets identical:", identical(s_cpp, s_r)))
ok(paste("cpp splits:", paste(s_cpp, collapse = " | ")))
ok(paste("r   splits:", paste(s_r,   collapse = " | ")))
