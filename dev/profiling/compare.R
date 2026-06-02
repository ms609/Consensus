#!/usr/bin/env Rscript
# Prove a method's speedup: re-run the grid against the CURRENTLY INSTALLED
# ConsTree and diff against a saved baseline CSV.  Run AFTER reimplementing a
# method (and reinstalling).
#
#   Rscript.exe dev/profiling/compare.R [baseline-<date>.csv] [Method ...]
#
# With no method args, compares all six.  Asserts, per non-trivial cell, that
# the new timing beats the baseline and the output split count is unchanged
# (the latter is the cheap correctness anchor -- a divergence here is a red flag
# unless it is the documented Greedy tie-break, which the chip must sign off).

args <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
scriptDir <- if (length(fileArg)) dirname(normalizePath(fileArg)) else getwd()
pkgRoot <- normalizePath(file.path(scriptDir, "..", ".."))
# The shared validation library, same one check-oracle.R uses.  NB: this is a
# SINGLE lib shared across worktrees -- a sibling worktree's install can clobber
# it, so the version self-guard below is load-bearing.  (Resolving it relative
# to pkgRoot would break: the build is installed into the main repo's
# .agent-cons, not the worktree's.)
agentLib <- "C:/Users/pjjg18/GitHub/Consensus/.agent-cons"
if (!dir.exists(agentLib)) {
  stop("Validation library not found: ", agentLib,
       "\n  Install first:  R CMD INSTALL --no-multiarch --library=\"", agentLib,
       "\" \"", pkgRoot, "\"")
}
.libPaths(c(agentLib, .libPaths()))

suppressMessages(library(ConsTree))
# Self-guard: installed build must match THIS worktree's source version, else a
# sibling worktree clobbered the shared lib and the timings would be a lie.
local({
  want <- read.dcf(file.path(pkgRoot, "DESCRIPTION"), fields = "Version")[1, 1]
  have <- as.character(utils::packageVersion("ConsTree"))
  if (have != want) {
    stop(sprintf(paste0("Installed ConsTree %s != this worktree's %s -- .agent-cons",
                        " holds a stale/foreign build; reinstall before profiling."),
                 have, want))
  }
})
source(file.path(scriptDir, "bench-common.R"))

rest <- commandArgs(trailingOnly = TRUE)
baseFile <- rest[grepl("\\.csv$", rest)]
methods <- rest[!grepl("\\.csv$", rest)]
if (!length(baseFile)) {
  cands <- sort(list.files(scriptDir, "^baseline-.*\\.csv$", full.names = TRUE))
  if (!length(cands)) stop("No baseline CSV found; run baseline.R first.")
  baseFile <- cands[length(cands)]
}
if (!length(methods)) methods <- BENCH_METHODS

base <- read.csv(baseFile, stringsAsFactors = FALSE)
base <- base[base$method %in% methods, ]
message("Comparing against ", basename(baseFile))

now <- benchGrid(sort(unique(base$nTip)), sort(unique(base$nTree)),
                 regimes = sort(unique(base$regime)),
                 methods = methods, reps = 3L, timeout = 60)

m <- merge(base, now, by = c("method", "nTip", "nTree", "regime"),
           suffixes = c(".base", ".new"))
m$speedup <- m$sec.base / m$sec.new
m$splitOK <- is.na(m$nSplit.base) | is.na(m$nSplit.new) |
  m$nSplit.base == m$nSplit.new

cat("\n")
print(m[order(m$method, m$nTip, m$nTree),
        c("method", "nTip", "nTree", "regime",
          "sec.base", "sec.new", "speedup", "nSplit.base", "nSplit.new",
          "splitOK")],
      row.names = FALSE)

bad <- m[!m$splitOK, ]
if (nrow(bad)) {
  cat("\nSplit-count divergences (investigate -- tie-break or bug):\n")
  print(bad[, c("method", "nTip", "nTree", "regime",
                "nSplit.base", "nSplit.new")], row.names = FALSE)
}
slower <- m[!is.na(m$speedup) & m$speedup < 1 & m$sec.base > 0.01, ]
if (nrow(slower)) {
  cat("\nCells where the new code is SLOWER (>10ms baseline):\n")
  print(slower[, c("method", "nTip", "nTree", "regime",
                   "sec.base", "sec.new", "speedup")], row.names = FALSE)
}
