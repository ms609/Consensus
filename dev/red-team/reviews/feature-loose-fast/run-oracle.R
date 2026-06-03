# External-reviewer oracle run for feature/loose-fast.
# Uses the freshly-built temp-lib install of the BRANCH source.
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/loose-fast/.rt-lib", .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
source("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/oracle.R")

# Build-identity check: confirm looseConsensusCpp is the real ported body,
# not a stub. A stub would return "" or a star; a real port resolves splits.
cat("== build identity ==\n")
b <- RootTree(RandomTree(paste0("t", 1:12), root = TRUE), "t1")
idem <- Loose(structure(list(b, b, b), class = "multiPhylo"))
cat(sprintf("  Loose(t,t,t) nSplits=%d (binary base nSplits=%d) -> %s\n",
            NSplits(idem), NSplits(b),
            if (NSplits(idem) == NSplits(b)) "REAL PORT" else "STUB/BROKEN"))

fails <- 0L
chk <- function(ok, msg) { cat(sprintf("  %-55s %s\n", msg, if (ok) "MATCH" else "*** DIFFER ***")); if (!ok) fails <<- fails + 1L }

cmp <- function(mine, fact, labels) setequal(SplitSet(mine, labels), SplitSet(fact, labels))
cmpPol <- function(mine, fact, labels) {
  pol <- function(tr) if (NSplits(tr) == 0L) character(0L) else
    as.character(PolarizeSplits(as.Splits(tr, tipLabels = labels)))
  setequal(pol(mine), pol(fact))
}

cat("\n== binary datasets (exact) ==\n")
datasets <- list(
  "random  n9  k21" = ape::as.phylo(0:20, 9),
  "random  n10 k31" = ape::as.phylo(0:30, 10),
  "conflict n8  k7" = ape::as.phylo(c(0,0,0,1,2,53,99), 8))
for (dn in names(datasets)) {
  trees <- datasets[[dn]]; labels <- TipLabels(trees[[1]])
  for (rt in c(0L,1L))
    chk(cmp(Loose(trees), FactConsensus(trees, "loose", rooted = rt), labels),
        sprintf("%s rooted=%d", dn, rt))
}

cat("\n== POLYTOMY / op==1 path (the path no shipped test pins exactly) ==\n")
polytomySets <- list(
  "n8 trichotomies" = c(
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2),(t3,t4,t5),(t6,t7,t8));")),
  "n9 nested polytomies" = c(
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3),t4,(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),t7,(t8,t9));")))
for (dn in names(polytomySets)) {
  trees <- polytomySets[[dn]]; labels <- TipLabels(trees[[1]])
  for (rt in c(0L,1L)) {
    mine <- Loose(trees); fact <- FactConsensus(trees, "loose", rooted = rt)
    chk(cmpPol(mine, fact, labels), sprintf("%s rooted=%d nSplit=%d", dn, rt, NSplits(mine)))
    # Also report the PLAIN comparison to see whether cmp would have flagged it
    cat(sprintf("        (plain cmp: %s)\n", if (cmp(mine, fact, labels)) "match" else "differ"))
  }
}

cat("\n== Loose at scale n>60 (op==1 heavy via perturbation) ==\n")
for (n in c(80L, 137L)) {
  set.seed(n + 1000L)
  labs <- paste0("t", seq_len(n))
  base <- RootTree(RandomTree(labs, root = TRUE), labs[[1]])
  trees <- structure(lapply(1:8, function(i) {
    tr <- base
    for (s in 1:3) { ij <- sample.int(n, 2L); tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)] }
    RootTree(RenumberTips(tr, labs), labs[[1]])
  }), class = "multiPhylo")
  mine <- Loose(trees)
  chk(cmp(mine, FactConsensus(trees, "loose", rooted = 1L), labs),
      sprintf("n=%d nSplit=%d", n, NSplits(mine)))
}

cat(sprintf("\n== TOTAL FAILURES: %d ==\n", fails))
if (fails > 0L) quit(status = 1L)
