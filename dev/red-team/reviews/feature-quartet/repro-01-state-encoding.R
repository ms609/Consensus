# repro-01-state-encoding.R
#
# Verifies that the two independent state-encoding paths in Quartet.cpp agree:
#   (a) build_quartet_profile's hand-written if-else (called per-tree, raw splits)
#   (b) quartet_state_from_sides (called per-split in add_benefit / do_add)
#
# We do this indirectly: construct 4-tip trees with KNOWN quartet state,
# call Quartet() on them, and verify both that the result has the expected
# state AND that the internal loss arithmetic gives the expected numeric values.
#
# Specifically, for n=4 there is exactly 1 quartet {t1,t2,t3,t4}.  The three
# binary topologies are:
#   state 1 (il|jk): t1,t4 | t2,t3  → tree ((t1,t4),(t2,t3))
#   state 2 (jl|ik): t2,t4 | t1,t3  → tree ((t2,t4),(t1,t3))
#   state 3 (kl|ij): t3,t4 | t1,t2  → tree ((t3,t4),(t1,t2))
#
# (Labels here are sorted alphabetically by TreeTools: t1 < t2 < t3 < t4)
#
# Usage: Rscript dev/red-team/reviews/feature-quartet/repro-01-state-encoding.R

suppressMessages({library(ConsTree); library(TreeTools); library(ape)})

tt <- function(x) read.tree(text = x)
nwk1 <- "((t1,t4),(t2,t3));"   # state 1: t1,t4 together
nwk2 <- "((t2,t4),(t1,t3));"   # state 2: t2,t4 together
nwk3 <- "((t3,t4),(t1,t2));"   # state 3: t3,t4 together

ok  <- TRUE
cat("-- State 1 wins (3x nwk1, 1x nwk2, 1x nwk3) -> Quartet returns state 1 tree\n")
trees135 <- structure(list(tt(nwk1), tt(nwk1), tt(nwk1), tt(nwk2), tt(nwk3)),
                      class = "multiPhylo")
qc <- Quartet(trees135)
sp <- sort(unname(as.character(as.Splits(qc, TipLabels(trees135[[1]])))))
# State 1 tree has t1,t4 together: split 10000001 / 01111110 → "t1,t4"
has_t1t4 <- any(grepl("t1", sp) & grepl("t4", sp) & !grepl("t2", sp) & !grepl("t3", sp))
cat("  Split contains {t1,t4}:", has_t1t4, "\n")
if (!has_t1t4) {
  cat("  FAIL: expected state-1 resolution (t1,t4 together)\n"); ok <- FALSE
} else {
  cat("  PASS\n")
}

cat("-- 2-2-1 tie (2x nwk1, 2x nwk2, 1x nwk3) -> Quartet returns a star (benefit < 0)\n")
# For k=5: state-1 count=2, state-2 count=2, state-3 count=1.
# add_benefit(state-1) = 2*2 - 5 = -1 < 0 -> no split added.
trees221 <- structure(list(tt(nwk1), tt(nwk1), tt(nwk2), tt(nwk2), tt(nwk3)),
                      class = "multiPhylo")
qc2 <- Quartet(trees221)
nsp <- NSplits(qc2)
cat("  NSplits (expect 0 = star):", nsp, "\n")
if (nsp != 0L) {
  cat("  FAIL: expected star tree for benefit < 0\n"); ok <- FALSE
} else {
  cat("  PASS\n")
}

cat("-- Tight majority (3x nwk1, 2x nwk2, 2x nwk3): k=7, count_1=3, benefit=2*3-7=-1 -> star\n")
trees_tight <- structure(
  list(tt(nwk1), tt(nwk1), tt(nwk1), tt(nwk2), tt(nwk2), tt(nwk3), tt(nwk3)),
  class = "multiPhylo")
qc3 <- Quartet(trees_tight)
nsp3 <- NSplits(qc3)
cat("  NSplits (expect 0 = star):", nsp3, "\n")
if (nsp3 != 0L) {
  cat("  FAIL: add_benefit should be -1, star expected\n"); ok <- FALSE
} else {
  cat("  PASS\n")
}

cat("-- Threshold (4x nwk1, 2x nwk2, 2x nwk3): k=8, count_1=4, benefit=2*4-8=0 -> star (tie)\n")
trees_thresh <- structure(
  list(tt(nwk1),tt(nwk1),tt(nwk1),tt(nwk1), tt(nwk2),tt(nwk2), tt(nwk3),tt(nwk3)),
  class = "multiPhylo")
qc4 <- Quartet(trees_thresh)
nsp4 <- NSplits(qc4)
cat("  NSplits (expect 0 = star, benefit=0 not >0):", nsp4, "\n")
if (nsp4 != 0L) {
  cat("  FAIL: benefit=0 should not trigger add (strict > 0 required)\n"); ok <- FALSE
} else {
  cat("  PASS\n")
}

cat("-- Benefit positive (5x nwk1, 2x nwk2, 2x nwk3): k=9, count_1=5, benefit=2*5-9=1 > 0\n")
trees_pos <- structure(
  list(tt(nwk1),tt(nwk1),tt(nwk1),tt(nwk1),tt(nwk1), tt(nwk2),tt(nwk2), tt(nwk3),tt(nwk3)),
  class = "multiPhylo")
qc5 <- Quartet(trees_pos)
nsp5 <- NSplits(qc5)
cat("  NSplits (expect 1 = resolved):", nsp5, "\n")
if (nsp5 != 1L) {
  cat("  FAIL: positive benefit should resolve the quartet\n"); ok <- FALSE
} else {
  cat("  PASS\n")
}

if (ok) {
  cat("\nAll state-encoding checks PASS.\n")
} else {
  cat("\nSome state-encoding checks FAILED.\n")
  quit(status = 1L)
}
