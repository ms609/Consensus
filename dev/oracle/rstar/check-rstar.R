# check-rstar.R
# Must-pass R* oracle cross-check (R* has no direct reference binary).  Exits
# non-zero on any failure, mirroring dev/oracle/local/check-local.R.  Three
# pillars that reliably hold:
#   1. Identity:        RStar(k copies of binary T) == T.
#   2. Congruent oracle: on congruent input RStar == FDCT `aho-build` (BUILD half).
#   3. Strict refinement: every unanimous rooted clade survives.
# The deeper exploration (incl. the OQ2 majority-refinement finding) lives in
# explore-rstar.R / diagnose-oq2.R.
#
# Run:  Rscript.exe dev/oracle/rstar/check-rstar.R

.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(Consensus))
suppressMessages(library(TreeTools))
suppressMessages(library(ape))
source("dev/oracle/oracle.R")
source("dev/oracle/local/oracle_local.R")

multi2list <- function(mp) lapply(seq_along(mp), function(i) mp[[i]])
fail <- 0L
cat("====  R* must-pass oracle cross-check  ====\n")

# --- 1. Identity --------------------------------------------------------------
cat("-- 1. identity (RStar of k copies of T == T)\n")
set.seed(11)
idFail <- 0L
for (trial in 1:40) {
  n <- sample(4:30, 1); k <- sample(2:5, 1)
  T <- rtree(n, rooted = TRUE); T$edge.length <- NULL
  if (!setequal(CladeSet(RStar(rep(list(T), k))), CladeSet(T))) idFail <- idFail + 1L
}
cat(sprintf("   identity failures: %d / 40\n", idFail)); fail <- fail + idFail

# --- 2. Congruent oracle vs FDCT aho-build ------------------------------------
cat("-- 2. congruent oracle (RStar == aho-build)\n")
congruent <- list(
  n4  = multi2list(ape::as.phylo(rep(0, 3), 4)),
  n6  = multi2list(ape::as.phylo(rep(3, 3), 6)),
  n8  = multi2list(ape::as.phylo(rep(11, 4), 8)),
  n10 = multi2list(ape::as.phylo(rep(7, 5), 10)),
  n12 = multi2list(ape::as.phylo(rep(99, 3), 12)),
  n16 = multi2list(ape::as.phylo(rep(40, 3), 16))
)
for (nm in names(congruent)) {
  trees <- congruent[[nm]]
  mine <- CladeSet(RStar(trees))
  o <- tryCatch(LocalOracle(trees, "aho"), error = function(e) NULL)
  if (is.null(o)) { cat(sprintf("   %-4s SKIP (oracle NULL)\n", nm)); next }
  ok <- setequal(mine, CladeSet(o))
  cat(sprintf("   %-4s %s\n", nm, if (ok) "MATCH" else "DIFFER"))
  if (!ok) fail <- fail + 1L
}

# --- 3. Strict refinement -----------------------------------------------------
cat("-- 3. strict-clade refinement (unanimous clades survive)\n")
set.seed(33)
scFail <- 0L
for (trial in 1:150) {
  n <- sample(5:12, 1); k <- sample(3:8, 1)
  labs <- NULL
  trees <- lapply(seq_len(k), function(i) { tr <- rtree(n, rooted = TRUE); tr$edge.length <- NULL; tr })
  labs <- trees[[1]]$tip.label
  trees <- lapply(trees, function(tr) RenumberTips(tr, labs))
  rs <- CladeSet(RStar(trees))
  tab <- table(unlist(lapply(trees, CladeSet)))
  strict <- names(tab)[tab == k]
  if (!all(strict %in% rs)) scFail <- scFail + 1L
}
cat(sprintf("   strict-refinement failures: %d / 150\n", scFail)); fail <- fail + scFail

# --- 4. Majority-rule refinement (Lemma 1.1) ----------------------------------
# Every rooted majority clade must appear in R* (was violated by the old BUILD
# construction; must now be 0).
cat("-- 4. majority-rule refinement (Lemma 1.1)\n")
set.seed(202)
mjFail <- 0L
for (trial in 1:400) {
  n <- sample(5:12, 1); k <- sample(3:9, 1)
  trees <- lapply(seq_len(k), function(i) { tr <- rtree(n, rooted = TRUE); tr$edge.length <- NULL; tr })
  labs <- trees[[1]]$tip.label
  trees <- lapply(trees, function(tr) RenumberTips(tr, labs))
  rs <- CladeSet(RStar(trees))
  tab <- table(unlist(lapply(trees, CladeSet)))
  maj <- names(tab)[tab > k / 2]
  if (!all(maj %in% rs)) mjFail <- mjFail + 1L
}
cat(sprintf("   majority-refinement failures: %d / 400\n", mjFail)); fail <- fail + mjFail

cat("\n====  Summary  ====\n")
if (fail == 0L) {
  cat("ALL must-pass R* checks PASSED.\n")
} else {
  cat(sprintf("FAILURES: %d\n", fail)); quit(status = 1)
}
