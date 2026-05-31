# BLOCKER repro #1: BHVDistance() hangs (infinite loop) on a legal tree with a
# zero-length INTERIOR edge that is incompatible with the other tree's split.
# Cause: min_weight_vc() in src/bhv.cpp normalises vertex weights by sA = sum(wA)
# (line ~68-75). When a GTP group's A-side squared lengths sum to 0, cap = wA/sA
# = 0/0 = NaN. NaN capacities make BFS find no path and the cover degenerate, so
# gtp_no_common()'s split never reduces the subproblem -> non-terminating queue.
#
# RUN: Rscript repro-01-zerolen-hang.R   (will NOT return; wrap in `timeout`)
# STATUS: RAN under `timeout 20` -> process killed (EXIT 143), i.e. confirmed hang.
.libPaths(c('.agent-rev', .libPaths()))
suppressMessages(library(Consensus))
A <- ape::read.tree(text="(0:1,(1:0,2:0):0,3:1,4:1,5:1);")  # interior split {1,2}, length 0
B <- ape::read.tree(text="(0:1,(2:1,3:1):1,1:1,4:1,5:1);")  # {2,3}:1 incompatible with {1,2}
cat("calling BHVDistance (expect hang) ...\n"); flush.console()
cat("d =", BHVDistance(A, B), "\n")  # never reached
