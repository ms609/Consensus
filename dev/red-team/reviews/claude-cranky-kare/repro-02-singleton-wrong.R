# MAJOR repro #2: BHVDistance() returns a WRONG (too large) distance for a tree
# containing an internal degree-2 (singleton) node, because .TreeToBHV() does not
# collapse singletons. Such a node yields two edges inducing the SAME bipartition,
# so the tree is given a duplicated split coordinate and lands at the wrong point
# in BHV space. A tree and its collapse.singles() equivalent are mathematically
# identical yet report distance 2*sqrt(2).
#
# Reachable in practice via ape::drop.tip(..., collapse.singles = FALSE).
# STATUS: RAN -> prints 2.828427 (should be 0).
.libPaths(c('.agent-rev', .libPaths()))
suppressMessages(library(Consensus)); suppressMessages(library(ape))
B  <- read.tree(text="(0:1,((1:1,(2:1,(3:1,(4:1,5:1):1):1):1):2):1,6:1,7:1);")  # singleton node
Bc <- collapse.singles(B)                                                        # identical tree
cat("has singleton:", any(table(B$edge[,1]) == 1), "\n")
cat("d(B, collapse.singles(B)) =", BHVDistance(B, Bc), "  (correct answer: 0)\n")
# Also: tree_at endpoint corruption
A <- read.tree(text="(0:1,(1:1,((2:1,3:1):2,(4:1,5:1):2):3):1,(6:1,7:1):4);")
cat("d(.BHVTreeAt(A,B,1), B) =", BHVDistance(Consensus:::.BHVTreeAt(A,B,1), B),
    "  (should be 0)\n")
