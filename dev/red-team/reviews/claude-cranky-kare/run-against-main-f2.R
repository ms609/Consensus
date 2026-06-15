# Run cranky-kare F2 (singleton -> wrong distance) against CURRENT main source.
# Loads the package from source (compiles src/) so we test main's bhv.cpp, not a
# stale install or the deleted branch's code.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(ape))
cat("=== F2: internal singleton node -> duplicated split ===\n")
B  <- read.tree(text="(0:1,((1:1,(2:1,(3:1,(4:1,5:1):1):1):1):2):1,6:1,7:1);")
Bc <- collapse.singles(B)
cat("has singleton:", any(table(B$edge[,1]) == 1), "\n")
d1 <- BHVDistance(B, Bc)
cat("d(B, collapse.singles(B)) =", d1, "  (correct answer: 0)\n")
A <- read.tree(text="(0:1,(1:1,((2:1,3:1):2,(4:1,5:1):2):3):1,(6:1,7:1):4);")
d2 <- BHVDistance(ConsTree:::.BHVTreeAt(A, B, 1), B)
cat("d(.BHVTreeAt(A,B,1), B) =", d2, "  (should be 0)\n")
cat("=== F2 VERDICT:", if (isTRUE(all.equal(d1, 0)) && isTRUE(all.equal(d2, 0)))
      "FIXED (both ~0)" else "LIVE (nonzero distance for identical trees)", "===\n")
