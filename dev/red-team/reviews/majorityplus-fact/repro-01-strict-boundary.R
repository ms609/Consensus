# repro-01: confirm MajorityPlus enforces STRICT > (displayed > contradicted),
# not >=, and that the boundary matches the FACT reference exactly.
# RAN by reviewer on 2026-06-02 against a CLEAN-recompiled build in an isolated
# library: result NSplits(twoToTwo)==0 (drop) and NSplits(twoToOne)==2 (keep),
# both FACT-exact. A >= implementation would give NSplits(twoToTwo)==1.
#
# Usage: source from the worktree root (or pass --file=<path> to Rscript).
#   R CMD INSTALL --no-multiarch --library=.agent-cons .
#   Rscript dev/red-team/reviews/majorityplus-fact/repro-01-strict-boundary.R
.rArgs <- commandArgs(FALSE)
.rFile <- sub("^--file=", "", grep("^--file=", .rArgs, value = TRUE))
.rDir  <- if (length(.rFile)) dirname(normalizePath(.rFile)) else getwd()
.rRoot <- normalizePath(file.path(.rDir, "..", "..", "..", ".."))  # worktree root
.rLib  <- file.path(.rRoot, ".agent-cons")
if (!dir.exists(.rLib)) {
  stop("Validation library not found: ", .rLib,
       "\n  Install first: R CMD INSTALL --no-multiarch --library=.agent-cons <worktree>")
}
.libPaths(c(.rLib, .libPaths()))
suppressMessages({library(ConsTree); library(TreeTools)})
source(file.path(.rRoot, "dev", "oracle", "oracle.R"))  # FactConsensus, SplitSet
cmp <- function(mine, fact, labels) setequal(SplitSet(mine, labels), SplitSet(fact, labels))

labels    <- paste0("t", 1:5)
display   <- ape::read.tree(text = "((t1, t2), (t3, t4), t5);")  # shows {t1,t2}
conflict  <- ape::read.tree(text = "((t1, t3), (t2, t4), t5);")
conflict2 <- ape::read.tree(text = "((t1, t4), (t2, t3), t5);")

# displayed 2, contradicted 2 -> 2 is NOT > 2 -> {t1,t2} dropped (star).
twoToTwo <- structure(list(display, display, conflict, conflict2), class = "multiPhylo")
stopifnot(NSplits(MajorityPlus(twoToTwo)) == 0L)
stopifnot(cmp(MajorityPlus(twoToTwo),
              FactConsensus(twoToTwo, "majorityPlus", rooted = 1L), labels))

# displayed 2, contradicted 1 -> 2 > 1 -> {t1,t2} kept.
twoToOne <- structure(list(display, display, conflict), class = "multiPhylo")
stopifnot(NSplits(MajorityPlus(twoToOne)) == 2L)
stopifnot(cmp(MajorityPlus(twoToOne),
              FactConsensus(twoToOne, "majorityPlus", rooted = 1L), labels))

cat("repro-01 PASS: strict > boundary is FACT-exact\n")
