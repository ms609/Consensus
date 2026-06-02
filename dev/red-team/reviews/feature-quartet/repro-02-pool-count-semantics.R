# repro-02-pool-count-semantics.R
#
# Demonstrates the pool.count semantics gap (F2 in review.md):
# pool.count tracks raw split-matrix row count, not tree-membership count.
#
# The gap is normally invisible for valid binary trees (each split appears
# exactly once per tree's as.Splits() matrix).  This script:
#   1. Shows that for valid input, pool.count == tree-membership count (normal path).
#   2. Constructs a hand-crafted splits_list with a within-tree duplicate row to
#      show the overcounting.  The C++ core is called directly to avoid the R
#      wrapper's valid-tree requirement.
#
# Usage: Rscript dev/red-team/reviews/feature-quartet/repro-02-pool-count-semantics.R

suppressMessages({library(ConsTree); library(TreeTools); library(ape)})

cat("-- 1. Normal path: pool.count == tree-membership count\n")
# Two distinct 5-tip trees.  The split {t1,t2}|{t3,t4,t5} appears in tree 1
# but not tree 2.  Majority threshold half = 1.  pool.count should be 1 -> not majority.
t1 <- read.tree(text = "((t1,t2),(t3,t4,t5));")
t2 <- read.tree(text = "((t1,t3),(t2,t4,t5));")
trees <- structure(list(t1, t2), class = "multiPhylo")
labs <- TipLabels(t1)
# With init="majority" and n_tree=2, half=1, a split needs count > 1 to be majority.
# {t1,t2} appears in only 1 tree -> count 1 -> not majority -> init starts from empty.
qc_maj <- Quartet(trees, init = "majority")
qc_emp <- Quartet(trees, init = "empty")
cat("  Majority-init splits:", NSplits(qc_maj),
    "  Empty-init splits:", NSplits(qc_emp), "\n")
# Both should give the same result since no split is majority.
cat("  Same result:", isTRUE(all.equal(
  sort(unname(as.character(as.Splits(qc_maj, labs)))),
  sort(unname(as.character(as.Splits(qc_emp, labs))))
)), "\n")

cat("\n-- 2. Degenerate path: hand-crafted duplicate row demonstrates overcounting\n")
# Build a 5-tip split matrix for tree1 where the SAME canonical split appears twice.
# Canonical split: tip 0 on side 0 (bit 0 = 0), tips 2,3,4 on side 1.
# Binary: 0b11100 = 0x1C (tips 2,3,4 set; tips 0,1 not set; n_bytes=1 for 5 tips).
# After canonicalization (bit 0 already 0): canonical = 0x1C.
# Row 1: 0x1C (canonical)
# Row 2: 0x03 (complement of 0x1C in 5 bits; bit 0 = 1 -> after flip = 0x1C again)
# Two different raw rows -> same canonical -> pool.count inflated.
splits_dup <- structure(
  matrix(as.raw(c(0x1C, 0x03)), nrow = 2, ncol = 1),
  nTip = 5L, tip.label = labs, class = "Splits"
)
splits_single <- structure(
  matrix(as.raw(c(0x1C)), nrow = 1, ncol = 1),
  nTip = 5L, tip.label = labs, class = "Splits"
)
# A second tree with a different split, count=1.
splits_other <- structure(
  matrix(as.raw(c(0x12)), nrow = 1, ncol = 1),  # tips 1,4 on side 1
  nTip = 5L, tip.label = labs, class = "Splits"
)

# With n_tree=2 and half=1:
#   Normal (no dup): split {2,3,4} count=1 -> not majority -> no majority init effect.
#   Dup path: split {2,3,4} count=2 -> count > half=1 -> classified as majority!
res_dup    <- ConsTree:::cpp_quartet_consensus(
  list(unclass(splits_dup), unclass(splits_other)), 5L,
  init_majority = TRUE, init_extended = FALSE, greedy_best_flag = FALSE)
res_normal <- ConsTree:::cpp_quartet_consensus(
  list(unclass(splits_single), unclass(splits_other)), 5L,
  init_majority = TRUE, init_extended = FALSE, greedy_best_flag = FALSE)

n_incl_dup    <- sum(res_dup$included)
n_incl_normal <- sum(res_normal$included)

cat("  Included splits with duplicate row:  ", n_incl_dup, "\n")
cat("  Included splits without duplicate:   ", n_incl_normal, "\n")
if (n_incl_dup == n_incl_normal) {
  cat("  (Results match — overcounting did not change the outcome here)\n")
  cat("  Note: the pool.count semantics gap is still present; this particular\n")
  cat("  input happens not to change the final result.\n")
} else {
  cat("  OVERCOUNTING EFFECT: duplicate row inflated count and changed init phase.\n")
  cat("  This is the F2 bug documented in review.md.\n")
}

cat("\npool-count-semantics repro complete.\n")
