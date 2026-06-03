#!/usr/bin/env Rscript
# Capture pre-change timing baselines for the six audited consensus methods,
# against the CURRENTLY INSTALLED ConsTree.  Run this BEFORE reimplementing a
# method, so compare.R can later prove the speedup.
#
#   Rscript.exe dev/profiling/baseline.R
#
# Writes dev/profiling/baseline-<date>.csv (git-tracked).  Override the date
# stamp (Date.now() is unavailable to keep runs reproducible) via the first arg.

# Locate package root from this script's directory (dev/profiling/..).
args <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
scriptDir <- if (length(fileArg)) dirname(normalizePath(fileArg)) else getwd()
pkgRoot <- normalizePath(file.path(scriptDir, "..", ".."))

# Prefer the isolated agent library if present (matches the oracle/test setup).
agentLib <- file.path(pkgRoot, ".agent-cons")
if (dir.exists(agentLib)) .libPaths(c(agentLib, .libPaths()))

suppressMessages(library(ConsTree))
source(file.path(scriptDir, "..", "oracle", "build-identity.R"))
assertConsTreeBuild(pkgRoot = pkgRoot)
source(file.path(scriptDir, "bench-common.R"))

stamp <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(stamp) || !nzchar(stamp)) stamp <- format(Sys.Date())

# Moderate grid that completes in ~10-20 min while still exposing the O(s^2) R
# pipeline's pain (n>=100).  Widen freely for a deeper sweep once the fast paths
# land.  Timing-out cells cost ~one `timeout` each (timeCall breaks on failure).
nTipGrid  <- c(10L, 50L, 100L, 200L)
nTreeGrid <- c(10L, 50L)

message("ConsTree baseline -- package version ",
        as.character(utils::packageVersion("ConsTree")))
res <- benchGrid(nTipGrid, nTreeGrid,
                 regimes = c("independent", "perturbed"),
                 methods = BENCH_METHODS,
                 reps = 2L, timeout = 20)

outFile <- file.path(scriptDir, sprintf("baseline-%s.csv", stamp))
write.csv(res, outFile, row.names = FALSE)
message("\nWrote ", outFile, " (", nrow(res), " rows)")
