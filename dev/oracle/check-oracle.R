# Cross-validate the R consensus methods against the reference FACT binary.
# Run with: Rscript dev/oracle/check-oracle.R
.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
source("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/oracle.R")

cmp <- function(mine, fact, labels) {
  setequal(SplitSet(mine, labels), SplitSet(fact, labels))
}

# Make the oracle an ASSERTION, not eyeball-only: count divergences and exit
# non-zero at the end, so a regression is caught by CI / a non-zero $? rather
# than scrolling past a printed "*** DIFFER ***".
failCount <- 0L
mark <- function(ok) {
  if (!isTRUE(ok)) failCount <<- failCount + 1L
  if (isTRUE(ok)) "MATCH" else "*** DIFFER ***"
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
                  mark(ok)))
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
              mark(setequal(cm, cf))))
}

# Adams at scale.  The small datasets above barely exercise the centroid-path
# iteration of the JLS2017 algorithm; these larger and more skewed inputs do.
# Yule-style random rooted trees (independent + perturbed) at n = 50, 137, and
# caterpillars (a single long spine -- the algorithm's best case and the slow
# recursion's worst case).  Adams is unique, so each must be clade-exact vs the
# slow fact.exe (rule 512, rooted = 1).
catTree <- function(labs) {
  nwk <- labs[length(labs)]
  for (i in (length(labs) - 1L):1L) nwk <- paste0("(", labs[[i]], ",", nwk, ")")
  ape::read.tree(text = paste0(nwk, ";"))
}
cat("\n== Adams at scale (rooted=1, clade comparison) ==\n")
adamsBig <- list()
set.seed(50L)
labs50 <- paste0("t", seq_len(50L))
adamsBig[["indep   n50  k20"]] <- structure(
  lapply(1:20, function(i) RandomTree(labs50, root = TRUE)), class = "multiPhylo")
base50 <- RandomTree(labs50, root = TRUE)
adamsBig[["perturb n50  k20"]] <- structure(lapply(1:20, function(i) {
  tr <- base50
  for (s in 1:3) {
    ij <- sample.int(50L, 2L)
    tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
  }
  tr
}), class = "multiPhylo")
set.seed(137L)
labs137 <- paste0("t", seq_len(137L))
adamsBig[["indep   n137 k10"]] <- structure(
  lapply(1:10, function(i) RandomTree(labs137, root = TRUE)), class = "multiPhylo")
set.seed(7L)
labs40 <- paste0("t", seq_len(40L))
cat1 <- catTree(labs40)
adamsBig[["caterpillar identical n40 k4"]] <-
  structure(list(cat1, cat1, cat1, cat1), class = "multiPhylo")
adamsBig[["caterpillar distinct  n40 k4"]] <- structure(
  list(cat1, catTree(sample(labs40)), catTree(sample(labs40)),
       catTree(sample(labs40))), class = "multiPhylo")
for (dn in names(adamsBig)) {
  trees <- adamsBig[[dn]]
  mine <- Adams(trees)
  fact <- FactConsensus(trees, "adams", rooted = 1L)
  cm <- CladeSet(mine)
  cf <- CladeSet(fact)
  cat(sprintf("  %-26s mine=%3d fact=%3d  %s\n", dn, length(cm), length(cf),
              mark(setequal(cm, cf))))
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
              mark(ok)))
}

# Loose at scale (n > 60).  Unlike greedy, the loose fast path does NO leaf-set
# bit-packing (it is purely structural: Day's labelling + consecutive-range /
# DEPTH queries), so there is no LEN > 1 word-arithmetic to exercise -- this
# block instead validates the looseMerge / contract pipeline on a large instance.
# And because the loose consensus is UNIQUE (every split compatible with all
# inputs, no frequency tie-break), it is FACT-exact at every n, so we assert an
# exact fact.exe match directly (no "same rooting" caveat -- SplitSet compares
# unrooted bipartitions and ignores the trivial root split).  Two checks:
# (a) idempotence Loose(list(t,t,t)) == t recovers a fully resolved binary tree;
# (b) a congruent (perturbed) set keeps real splits AND matches fact.exe exactly
#     (independent random trees would share no all-compatible split -> the star
#     tree, making the assertion vacuous).
cat("\n== Loose at scale (n > 60) ==\n")
for (n in c(80L, 137L)) {
  set.seed(n + 1000L)
  labs <- paste0("t", seq_len(n))
  base <- RootTree(RandomTree(labs, root = TRUE), labs[[1]])
  idem <- setequal(SplitSet(Loose(structure(list(base, base, base),
                                            class = "multiPhylo")), labs),
                   SplitSet(base, labs))
  trees <- structure(lapply(1:8, function(i) {
    tr <- base
    for (s in 1:3) {
      ij <- sample.int(n, 2L)
      tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
    }
    RootTree(RenumberTips(tr, labs), labs[[1]])
  }), class = "multiPhylo")
  mine <- Loose(trees)
  ok <- cmp(mine, FactConsensus(trees, "loose", rooted = 1L), labs)
  cat(sprintf("  n=%-3d  idempotent: %-5s   nSplit=%-3d  FACT-exact: %s\n",
              n, idem, NSplits(mine),
              mark(ok)))
}

# Loose with POLYTOMOUS inputs.  Every dataset above (and in test-loose.R) is
# binary, so each input tree B has only 2-child nodes -- but looseMerge's op == 1
# tree-construction inserts new vertices AMONG a B node's children via the
# BEFORE/AFTER + pid bookkeeping, a path that is only fully exercised when an
# input has a node with > 2 children.  Loose is deterministic (compatible-with-
# all), so this must be FACT-exact; a divergence is a real bug.  Checked under
# both rooted flags, as the binary datasets are.
#
# NB: compare CANONICAL (polarised) splits here.  A polytomous loose consensus is
# rooted differently by Loose() (.RootLikeFirst) and by fact.exe, and as.Splits()
# orients each bipartition by descendant side -- which flips with the root -- so
# the plain SplitSet comparison `cmp` (fine for the binary datasets, where mine
# and fact happen to orient alike) reports SPURIOUS differences here.
cat("\n== Loose with polytomous inputs ==\n")
cmpPol <- function(mine, fact, labels) {
  pol <- function(tr) if (NSplits(tr) == 0L) character(0L) else
    as.character(PolarizeSplits(as.Splits(tr, tipLabels = labels)))
  setequal(pol(mine), pol(fact))
}
polytomySets <- list(
  "n8 trichotomies" = c(
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2,t3),(t4,t5),(t6,t7,t8));"),
    ape::read.tree(text = "((t1,t2),(t3,t4,t5),(t6,t7,t8));")
  ),
  "n9 nested polytomies" = c(
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3),t4,(t5,t6),(t7,t8,t9));"),
    ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6),t7,(t8,t9));")
  )
)
for (dn in names(polytomySets)) {
  trees <- polytomySets[[dn]]
  labels <- TipLabels(trees[[1]])
  for (rt in c(0L, 1L)) {
    mine <- Loose(trees)
    fact <- FactConsensus(trees, "loose", rooted = rt)
    ok <- cmpPol(mine, fact, labels)
    cat(sprintf("  %-22s rooted=%d  nSplit=%-2d  %s\n", dn, rt, NSplits(mine),
                mark(ok)))
  }
}

# Fail loud: a non-zero exit status turns this script into a real gate (CI / a
# scripted `Rscript ... || stop`), instead of a wall of text a regression could
# hide in.  (The strict rooted-flag DETERMINATION block above is diagnostic and
# deliberately not counted.)
cat("\n")
if (failCount > 0L) {
  cat(sprintf("*** %d oracle comparison(s) DIFFERED -- FAILING. ***\n", failCount))
  quit(status = 1L)
}
cat("All oracle comparisons MATCH.\n")
