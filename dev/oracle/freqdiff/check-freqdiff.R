# Cross-validate Consensus::Frequency() against the reference FDCT_new freqdiff
# binary on the same fixtures used in dev/oracle/check-oracle.R.
#
# Run with:
#   Rscript.exe dev/oracle/freqdiff/check-freqdiff.R
#
# Frequency-difference is an UNROOTED split-based method.  We compare
# canonical unrooted bipartition sets (via TreeTools::as.Splits) using
# setequal(), so topology labels don't matter.
#
# If the two sides DIFFER the script reports:
#   - which splits are in R only / C++ only
#   - the count of each split across the input trees
#   - a diagnosis of whether it is a tie-break artifact (equal-frequency
#     conflicting splits) vs a genuine discrepancy

setTimeLimit(elapsed = 55, transient = FALSE)

.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(Consensus))
suppressMessages(library(TreeTools))
source("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/freqdiff/oracle_fd.R")

# ----- helper: split occurrence counts from a pool ----------------------------

.SplitCounts <- function(trees) {
  # Returns a named integer vector: canonical split string -> count
  labels <- TipLabels(trees[[1]])
  counts <- integer(0)
  for (tr in trees) {
    if (NSplits(tr) == 0L) next
    ss <- unname(as.character(as.Splits(tr, tipLabels = labels)))
    for (s in ss) {
      # Use names() lookup to avoid subscript-out-of-bounds on missing names
      prev <- counts[s]
      counts[s] <- if (is.na(prev)) 1L else prev + 1L
    }
  }
  counts
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ----- fixtures ---------------------------------------------------------------

datasets <- list(
  "random  n9  k21" = ape::as.phylo(0:20, 9),
  "random  n10 k31" = ape::as.phylo(0:30, 10),
  "conflict n8  k7" = ape::as.phylo(c(0, 0, 0, 1, 2, 53, 99), 8)
)

# ----- main loop --------------------------------------------------------------

all_pass <- TRUE

for (dn in names(datasets)) {
  trees  <- datasets[[dn]]
  labels <- TipLabels(trees[[1]])
  k      <- length(trees)

  cat("--", dn, "--\n")

  mine_tree <- tryCatch(Frequency(trees), error = function(e) {
    cat(sprintf("  Frequency() ERROR: %s\n", conditionMessage(e)))
    NULL
  })
  ref_tree  <- tryCatch(FreqDiffOracle(trees), error = function(e) {
    cat(sprintf("  FreqDiffOracle() ERROR: %s\n", conditionMessage(e)))
    NULL
  })

  if (is.null(mine_tree) || is.null(ref_tree)) {
    all_pass <- FALSE
    next
  }

  mine_splits <- SplitSetFD(mine_tree, labels)
  ref_splits  <- SplitSetFD(ref_tree,  labels)

  ok <- setequal(mine_splits, ref_splits)
  cat(sprintf("  Frequency() mine=%2d  freqdiff ref=%2d  %s\n",
              length(mine_splits), length(ref_splits),
              if (ok) "MATCH" else "*** DIFFER ***"))

  if (!ok) {
    all_pass <- FALSE
    only_mine <- setdiff(mine_splits, ref_splits)
    only_ref  <- setdiff(ref_splits,  mine_splits)

    # Split occurrence counts across input trees
    cnt <- .SplitCounts(trees)

    if (length(only_mine) > 0L) {
      cat("  Splits in R only:\n")
      for (s in only_mine) {
        v <- cnt[s]; v <- if (is.na(v)) 0L else v
        cat(sprintf("    %s  (count = %d / %d)\n", s, v, k))
      }
    }
    if (length(only_ref) > 0L) {
      cat("  Splits in C++ only:\n")
      for (s in only_ref) {
        v <- cnt[s]; v <- if (is.na(v)) 0L else v
        cat(sprintf("    %s  (count = %d / %d)\n", s, v, k))
      }
    }

    # Tie-break diagnosis: for each discrepant split, find conflicting splits
    # and check if any have the same count.
    all_splits <- union(mine_splits, ref_splits)
    pool <- tryCatch({
      sp_obj <- Consensus:::.PoolSplits(trees)
      list(splits = sp_obj$splits, counts = sp_obj$counts,
           members = sp_obj$members, labels = sp_obj$labels)
    }, error = function(e) NULL)

    if (!is.null(pool)) {
      canon <- unname(as.character(as.Splits(pool$splits,
                                             tipLabels = pool$labels)))
      cat("  Tie-break diagnosis (from R pool):\n")
      for (s in union(only_mine, only_ref)) {
        idx <- which(canon == s)
        if (length(idx) == 0L) {
          cat(sprintf("    %s  not in pool\n", s))
          next
        }
        cnt_s <- pool$counts[idx]
        # Round-trip through the logical matrix to preserve the Splits S4 class
        # (single-bracket `[` drops it; see AGENTS.md "Gotcha: Splits is S4").
        oneSplit <- as.Splits(pool$members[idx, , drop = FALSE],
                              tipLabels = pool$labels)
        compat <- as.logical(CompatibleSplits(oneSplit, pool$splits))
        rival_counts <- pool$counts[!compat]
        tie <- any(rival_counts == cnt_s)
        cat(sprintf("    %s  count=%d  max_rival=%d  tie=%s\n",
                    s, cnt_s,
                    if (length(rival_counts) > 0L) max(rival_counts) else 0L,
                    if (tie) "YES" else "no"))
      }
    }
  }
  cat("\n")
}

cat(if (all_pass) "ALL MATCH\n" else "SOME DIFFER — see above\n")
