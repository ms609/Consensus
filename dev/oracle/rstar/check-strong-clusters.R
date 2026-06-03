# check-strong-clusters.R  (dev-only)
# DEFINITIVE correctness gate for RStar(), straight from Lemma 1.1 (Jansson,
# Sung, Vu & Yiu 2016): the R* tree's clusters are EXACTLY the strong clusters
# of R_maj.  For small n we brute-force ALL 2^n subsets, test each for the
# strong-cluster property with an INDEPENDENT ape::mrca tally (different code
# path from src/rstar.cpp), and assert
#     CladeSet(RStar(trees)) == { A : A is a strong cluster of R_maj }.
# Exits non-zero on any mismatch.  Run:
#   Rscript.exe dev/oracle/rstar/check-strong-clusters.R

.libPaths(c(Sys.getenv("CONSTREE_LIB", "C:/Users/pjjg18/GitHub/Consensus/.agent-cons"), .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
suppressMessages(library(ape))
source("dev/oracle/oracle.R")  # CladeSet

makeStater <- function(tree) {
  tree2 <- tree; tree2$edge.length <- rep(1, nrow(tree$edge))
  dep <- ape::node.depth.edgelength(tree2)
  m <- ape::mrca(tree)
  function(x, y, z) {
    dxy <- dep[m[x, y]]; dxz <- dep[m[x, z]]; dyz <- dep[m[y, z]]
    if (dxy == dxz && dxy == dyz) return(NA_character_)
    if (dxy > dxz && dxy > dyz) return(paste(sort(c(x, y)), collapse = ","))
    if (dxz > dxy && dxz > dyz) return(paste(sort(c(x, z)), collapse = ","))
    paste(sort(c(y, z)), collapse = ",")
  }
}
uniquelyFavoured <- function(staters, x, y, Z) {
  votes <- vapply(staters, function(st) st(x, y, Z), character(1))
  votes <- votes[!is.na(votes)]
  if (!length(votes)) return(FALSE)
  tb <- table(votes); target <- paste(sort(c(x, y)), collapse = ",")
  if (is.na(tb[target]) || tb[target] == 0) return(FALSE)
  mx <- max(tb); tb[target] == mx && sum(tb == mx) == 1L
}
cladeRuleHolds <- function(staters, C, allTips) {
  outside <- setdiff(allTips, C)
  if (!length(outside)) return(TRUE)
  pairs <- utils::combn(C, 2)
  for (p in seq_len(ncol(pairs)))
    for (Z in outside)
      if (!uniquelyFavoured(staters, pairs[1, p], pairs[2, p], Z)) return(FALSE)
  TRUE
}

# Brute-force set of strong-cluster clade strings (size 2..n-1).
strongCladeStrings <- function(trees) {
  allTips <- trees[[1]]$tip.label
  n <- length(allTips)
  staters <- lapply(trees, makeStater)
  out <- character(0)
  for (mask in seq_len(2^n - 1L)) {
    idx <- which(as.integer(intToBits(mask))[seq_len(n)] == 1L)
    if (length(idx) < 2L || length(idx) > n - 1L) next
    A <- allTips[idx]
    if (cladeRuleHolds(staters, A, allTips)) out <- c(out, paste(sort(A), collapse = ","))
  }
  sort(unique(out))
}

alignTrees <- function(trees) {
  labs <- trees[[1]]$tip.label
  lapply(trees, function(tr) RenumberTips(tr, labs))
}

runBattery <- function(label, gen, nTrial) {
  fail <- 0L
  for (trial in seq_len(nTrial)) {
    trees <- gen()
    mine   <- sort(CladeSet(RStar(trees)))
    oracle <- strongCladeStrings(trees)
    if (!setequal(mine, oracle)) {
      fail <- fail + 1L
      cat(sprintf("  [%s] MISMATCH trial %d\n", label, trial))
      cat("    extra in RStar:   ", paste(setdiff(mine, oracle), collapse = " ; "), "\n")
      cat("    missing from RStar:", paste(setdiff(oracle, mine), collapse = " ; "), "\n")
      for (w in vapply(trees, write.tree, character(1))) cat("      ", w, "\n")
    }
  }
  cat(sprintf("  %-22s trials: %d   mismatches: %d\n", label, nTrial, fail))
  fail
}

cat("====  Brute-force strong-cluster oracle (Lemma 1.1)  ====\n")
set.seed(2024)
fail <- 0L

# Battery 1: binary input.
fail <- fail + runBattery("binary", function() {
  n <- sample(4:12, 1); k <- sample(2:9, 1)
  alignTrees(lapply(seq_len(k), function(i) {
    tr <- rtree(n, rooted = TRUE); tr$edge.length <- NULL; tr
  }))
}, 60L)

# Battery 2: PARTIALLY RESOLVED (non-binary) input -- exercises the mixed
# fans-abstain / resolved-count path (di2multi collapses ~short edges).  This is
# the OQ1 "fans have no impact" claim under realistic non-binary input.
fail <- fail + runBattery("partly-resolved", function() {
  n <- sample(5:12, 1); k <- sample(2:9, 1)
  alignTrees(lapply(seq_len(k), function(i) {
    tr <- ape::di2multi(rtree(n, rooted = TRUE), tol = runif(1, 0.15, 0.45))
    tr$edge.length <- NULL; tr
  }))
}, 60L)

# Battery 3: a mix of binary and non-binary trees in the same input set.
fail <- fail + runBattery("mixed-set", function() {
  n <- sample(5:11, 1); k <- sample(3:9, 1)
  alignTrees(lapply(seq_len(k), function(i) {
    tr <- rtree(n, rooted = TRUE)
    if (i %% 2L == 0L) tr <- ape::di2multi(tr, tol = 0.3)
    tr$edge.length <- NULL; tr
  }))
}, 60L)

cat("\n")
if (fail == 0L) {
  cat("  ALL batteries: RStar clusters == strong clusters of R_maj.\n")
} else {
  cat(sprintf("  TOTAL mismatches: %d\n", fail)); quit(status = 1)
}
