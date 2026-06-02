# Shared helpers for the ConsTree consensus-method benchmarks.
# Sourced by baseline.R (capture pre-change timings) and compare.R (prove a
# speedup against a saved baseline).  Dev-only; not shipped (`.Rbuildignore`).
#
# Timing uses base `system.time` (median of a few reps) rather than `bench`, to
# avoid adding a package dependency for a dev script -- the expected wins are
# order-of-magnitude, so elapsed-time resolution is ample.

suppressMessages({
  library(TreeTools)
})

# The six methods under audit.  Names match `ConsTree::`.
BENCH_METHODS <- c("Greedy", "Loose", "MajorityPlus",
                   "Frequency", "Adams", "RStar")

# Hard size caps documented in the package (skip cells that would error).
BENCH_CAPS <- c(RStar = 200L)

# Generate `nTree` trees on `nTip` leaves under one of two regimes:
#   "independent" -- k independent random topologies (incongruent; large split
#                    pool, the worst case for the O(s^2) R pipeline).
#   "perturbed"   -- one base topology, each replicate perturbed by a few random
#                    tip-label swaps (mostly congruent; small split pool).
# Uses TreeTools::RandomTree (per project convention; not ape::rtree).
makeTrees <- function(nTip, nTree, regime = c("independent", "perturbed"),
                      seed = 1L, nSwap = 2L) {
  regime <- match.arg(regime)
  set.seed(seed)
  labels <- paste0("t", seq_len(nTip))
  if (regime == "independent") {
    trees <- lapply(seq_len(nTree), function(i) RandomTree(labels, root = TRUE))
  } else {
    base <- RandomTree(labels, root = TRUE)
    trees <- lapply(seq_len(nTree), function(i) {
      tr <- base
      if (nSwap > 0L) for (s in seq_len(nSwap)) {
        ij <- sample.int(nTip, 2L)
        tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
      }
      tr
    })
  }
  structure(trees, class = "multiPhylo")
}

# Time `fn(trees)`: median elapsed over `reps`, guarded by a per-call elapsed
# `timeout` (seconds).  Returns the median seconds (NA on timeout/error) and the
# output split count (an inexpensive correctness anchor across baseline/after).
timeCall <- function(fn, trees, reps = 3L, timeout = 30) {
  times <- numeric(reps)
  out <- NULL
  failed <- FALSE
  for (i in seq_len(reps)) {
    res <- tryCatch({
      setTimeLimit(elapsed = timeout, transient = TRUE)
      on.exit(setTimeLimit(), add = TRUE)
      system.time(out <- fn(trees))[["elapsed"]]
    }, error = function(e) { failed <<- TRUE; NA_real_ })
    setTimeLimit()
    if (failed) break
    times[i] <- res
  }
  nSplit <- if (!failed && !is.null(out)) {
    tryCatch(as.integer(TreeTools::NSplits(out)), error = function(e) NA_integer_)
  } else NA_integer_
  list(sec = if (failed) NA_real_ else stats::median(times), nSplit = nSplit)
}

# Run every method over a (nTip x nTree x regime) grid against the currently
# installed ConsTree, returning a long data frame.  `methods` lets callers
# restrict to e.g. just "Greedy".
benchGrid <- function(nTip, nTree,
                      regimes = c("independent", "perturbed"),
                      methods = BENCH_METHODS,
                      reps = 3L, timeout = 30) {
  rows <- list()
  for (regime in regimes) for (nt in nTip) for (k in nTree) {
    trees <- makeTrees(nt, k, regime, seed = nt * 100000L + k)
    for (m in methods) {
      cap <- BENCH_CAPS[m]
      if (!is.na(cap) && nt > cap) {
        rows[[length(rows) + 1L]] <- data.frame(
          method = m, nTip = nt, nTree = k, regime = regime,
          sec = NA_real_, nSplit = NA_integer_, note = "exceeds cap",
          stringsAsFactors = FALSE)
        next
      }
      fn <- get(m, envir = asNamespace("ConsTree"))
      r <- timeCall(fn, trees, reps = reps, timeout = timeout)
      rows[[length(rows) + 1L]] <- data.frame(
        method = m, nTip = nt, nTree = k, regime = regime,
        sec = r$sec, nSplit = r$nSplit,
        note = if (is.na(r$sec)) "timeout/error" else "",
        stringsAsFactors = FALSE)
      message(sprintf("  %-13s n=%-4d k=%-4d %-11s  %s",
                      m, nt, k, regime,
                      if (is.na(r$sec)) "NA" else sprintf("%.4fs (%d splits)",
                                                          r$sec, r$nSplit)))
    }
  }
  do.call(rbind, rows)
}
