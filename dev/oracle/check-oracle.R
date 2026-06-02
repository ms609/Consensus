# Cross-validate the R consensus methods against the reference FACT binary.
# Run with: Rscript dev/oracle/check-oracle.R
.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
source("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/oracle.R")

cmp <- function(mine, fact, labels) {
  setequal(SplitSet(mine, labels), SplitSet(fact, labels))
}

datasets <- list(
  "random  n9  k21" = ape::as.phylo(0:20, 9),
  "random  n10 k31" = ape::as.phylo(0:30, 10),
  "conflict n8  k7" = ape::as.phylo(c(0, 0, 0, 1, 2, 53, 99), 8)
)
methods <- list(strict = Strict, majority = Majority, greedy = Greedy,
                loose = Loose, majorityPlus = MajorityPlus)

trees1 <- datasets[[1]]
labels1 <- TipLabels(trees1[[1]])
cat("== Determine FACT rooted flag (strict) ==\n")
for (rt in c(0L, 1L)) {
  ok <- cmp(Strict(trees1), FactConsensus(trees1, "strict", rooted = rt), labels1)
  cat(sprintf("  strict vs FACT rooted=%d : %s\n", rt, if (ok) "MATCH" else "differ"))
}

for (rt in c(0L, 1L)) {
  cat(sprintf("\n== Cross-validation (rooted=%d) ==\n", rt))
  for (dn in names(datasets)) {
    trees <- datasets[[dn]]
    labels <- TipLabels(trees[[1]])
    cat("--", dn, "--\n")
    for (mn in names(methods)) {
      mine <- methods[[mn]](trees)
      fact <- FactConsensus(trees, mn, rooted = rt)
      ok <- cmp(mine, fact, labels)
      cat(sprintf("  %-13s mine=%2d fact=%2d  %s\n",
                  mn, NSplits(mine), NSplits(fact),
                  if (ok) "MATCH" else "*** DIFFER ***"))
    }
  }
}

# Adams is a ROOTED method: validate against the classical (slow) Adams with
# rooted = 1 (each input tree's own root), comparing ROOTED CLADES, not splits.
cat("\n== Adams cross-validation (rooted=1, clade comparison) ==\n")
for (dn in names(datasets)) {
  trees <- datasets[[dn]]
  mine <- Adams(trees)
  fact <- FactConsensus(trees, "adams", rooted = 1L)
  cm <- CladeSet(mine)
  cf <- CladeSet(fact)
  cat(sprintf("  %-16s mine=%2d fact=%2d  %s\n", dn, length(cm), length(cf),
              if (setequal(cm, cf)) "MATCH" else "*** DIFFER ***"))
}

# Multi-word bitset path (n > 60: BUCKET_SIZE = 60, so LEN > 1 -- the word-index
# arithmetic (c-1)/60, %60 and the LEN-length OR-up/compare loops).  The datasets
# above all have LEN = 1, so the multi-word packing needs its own check.
#
# Greedy is tie-break sensitive, and FACT's tie-break depends on its exact
# internal clade representation (it roots at taxon 1's neighbour; the R wrapper
# roots on taxon 1's edge), so an exact match on tie-heavy random input is NOT
# expected -- greedy is "FACT-match up to tie-break" (as the previous R
# implementation also was; AGENTS.md).  Validate the multi-word path two ways
# that ARE exact: (a) idempotence Greedy(list(t,t,t)) == t exercises the pack/
# unpack with no ties; (b) with every tree identically rooted at taxon 1, mine
# and FACT see the same clades, isolating algorithm faithfulness from the
# rooting/tie-break artefact.
cat("\n== Multi-word bitset path (n > 60) ==\n")
for (n in c(80L, 137L)) {
  set.seed(n)
  labs <- paste0("t", seq_len(n))
  base <- RandomTree(labs, root = TRUE)
  idem <- setequal(SplitSet(Greedy(structure(list(base, base, base),
                                             class = "multiPhylo")), labs),
                   SplitSet(base, labs))
  trees <- structure(lapply(1:15, function(i)
                       RootTree(RandomTree(labs, root = TRUE), labs[[1]])),
                     class = "multiPhylo")
  ok <- cmp(Greedy(trees), FactConsensus(trees, "greedy", rooted = 1L), labs)
  cat(sprintf("  n=%-3d LEN=%d  idempotent: %-5s   FACT-exact (same rooting): %s\n",
              n, (n + 59L) %/% 60L, idem,
              if (ok) "MATCH" else "*** DIFFER ***"))
}
