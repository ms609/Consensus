# snapshot-legacy.R  (dev-only; run ONCE, BEFORE editing src/rstar.cpp)
# Captures CladeSet(RStar(.)) of the CURRENT (definition-exact, O(kn^3)+O(n^4),
# 200-leaf-capped) build over the deterministic legacy-grid, serialising it to
# dev/oracle/rstar/legacy-clades.rds.  check-vs-legacy.R later asserts the new
# fast build reproduces every entry exactly.  This is the strongest regression
# guard for R* (which has NO reference binary): the old code IS the reference.
#
# Run from the package root:  Rscript.exe dev/oracle/rstar/snapshot-legacy.R

.libPaths(c(Sys.getenv("CONSTREE_LIB", "C:/Users/pjjg18/GitHub/Consensus/.agent-cons"), .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })
source("dev/oracle/oracle.R")             # CladeSet()
source("dev/oracle/rstar/legacy-grid.R")  # legacyTrials()

# --- gate: refuse to run unless this is the OLD capped build ------------------
isOld <- tryCatch({ ConsTree:::rStarConsensus(list(), 201L); FALSE },
                  error = function(e) grepl("200 leaves", conditionMessage(e)))
if (!isTRUE(isOld)) {
  stop("snapshot-legacy.R must run against the OLD (200-leaf-capped) build; ",
       "the installed build does not error at n = 201. Aborting to avoid a ",
       "self-referential snapshot.")
}
cat("Legacy build confirmed (version ",
    as.character(packageVersion("ConsTree")), ", cap present).\n", sep = "")

trials <- legacyTrials()
cat("Snapshotting", length(trials), "tree-sets ...\n")
snap <- vector("list", length(trials))
names(snap) <- names(trials)
for (key in names(trials)) snap[[key]] <- sort(CladeSet(RStar(trials[[key]])))

saveRDS(snap, "dev/oracle/rstar/legacy-clades.rds")
nClades <- vapply(snap, length, integer(1))
cat(sprintf("Wrote dev/oracle/rstar/legacy-clades.rds  (%d trials; clades/trial: min %d, median %g, max %d)\n",
            length(snap), min(nClades), stats::median(nClades), max(nClades)))
