# diagnose-clade-rule.R  (dev-only)
# Measures BOTH directions in which BUILD-based RStar() can deviate from Degnan
# et al.'s strict clade rule (clade C exists iff for EVERY internal pair {x,y}
# and EVERY outside Z, the triplet (x,y)|Z is uniquely favoured).  Uses an
# INDEPENDENT ape::mrca-based triplet tally (a different code path from
# src/rstar.cpp).
#
#   OQ2 (UNDER-resolution): majority clades that RStar DROPS but that DO satisfy
#        the rule -> lost only to the flat-collapse of an inconsistent component.
#   OQ3 (OVER-resolution):  clades that RStar PRODUCES but that do NOT satisfy
#        the rule -> BUILD's exists-an-outgroup edge rule resolved a clade the
#        strict for-all-outgroup definition would leave unresolved.  This can
#        happen with NO inconsistency and NO collapse (e.g. n=4, only (ab)|c
#        favoured, all other triples tied -> BUILD emits {a,b}, but (ab)|d is a
#        tie, so Degnan would not).
#
# Run:  Rscript.exe dev/oracle/rstar/diagnose-clade-rule.R

.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
suppressMessages(library(ape))
source("dev/oracle/oracle.R")  # CladeSet

makeStater <- function(tree) {
  tree2 <- tree; tree2$edge.length <- rep(1, nrow(tree$edge))
  dep <- ape::node.depth.edgelength(tree2)
  m <- ape::mrca(tree)
  function(x, y, z) {
    dxy <- dep[m[x, y]]; dxz <- dep[m[x, z]]; dyz <- dep[m[y, z]]
    if (dxy == dxz && dxy == dyz) return(NA_character_)
    if (dxy > dxz && dxy > dyz) return(paste(sort(c(x, y)), collapse = ","))
    if (dxz > dxy && dxz > dyz) return(paste(sort(c(x, z)), collapse = ","))
    paste(sort(c(y, z)), collapse = ",")
  }
}
uniquelyFavoured <- function(staters, x, y, Z) {
  votes <- vapply(staters, function(st) st(x, y, Z), character(1))
  votes <- votes[!is.na(votes)]
  if (!length(votes)) return(FALSE)
  tb <- table(votes); target <- paste(sort(c(x, y)), collapse = ",")
  if (is.na(tb[target]) || tb[target] == 0) return(FALSE)
  mx <- max(tb); tb[target] == mx && sum(tb == mx) == 1L
}
cladeRuleHolds <- function(staters, C, allTips) {
  outside <- setdiff(allTips, C)
  if (!length(outside)) return(TRUE)
  pairs <- utils::combn(C, 2)
  for (p in seq_len(ncol(pairs))) {
    for (Z in outside) if (!uniquelyFavoured(staters, pairs[1, p], pairs[2, p], Z)) return(FALSE)
  }
  TRUE
}
alignTrees <- function(trees) {
  labs <- trees[[1]]$tip.label
  lapply(trees, function(tr) RenumberTips(tr, labs))
}
majClades <- function(trees) {
  k <- length(trees); tab <- table(unlist(lapply(trees, CladeSet)))
  names(tab)[tab > k / 2]
}

set.seed(202)  # same trials as explore-rstar.R section C
nTrial <- 400L
# OQ2 accounting
oq2Dropped <- 0L; oq2RuleHolds <- 0L; oq2RuleFails <- 0L
# OQ3 accounting
oq3Clades <- 0L; oq3Trials <- 0L; oq3Standalone <- 0L; producedTotal <- 0L
oq3Examples <- list()

for (trial in seq_len(nTrial)) {
  n <- sample(5:12, 1); k <- sample(3:9, 1)
  trees <- alignTrees(lapply(seq_len(k), function(i) {
    tr <- rtree(n, rooted = TRUE); tr$edge.length <- NULL; tr
  }))
  allTips <- trees[[1]]$tip.label
  rs <- tryCatch(CladeSet(RStar(trees)), error = function(e) NULL)
  if (is.null(rs)) next
  staters <- lapply(trees, makeStater)

  # OQ2: dropped majority clades
  dropped <- setdiff(majClades(trees), rs)
  for (cl in dropped) {
    oq2Dropped <- oq2Dropped + 1L
    if (cladeRuleHolds(staters, strsplit(cl, ",")[[1]], allTips)) {
      oq2RuleHolds <- oq2RuleHolds + 1L
    } else oq2RuleFails <- oq2RuleFails + 1L
  }

  # OQ3: produced clades that violate the strict rule
  producedTotal <- producedTotal + length(rs)
  over <- rs[!vapply(rs, function(cl) cladeRuleHolds(staters, strsplit(cl, ",")[[1]], allTips), logical(1))]
  if (length(over)) {
    oq3Clades <- oq3Clades + length(over)
    oq3Trials <- oq3Trials + 1L
    if (!length(dropped)) oq3Standalone <- oq3Standalone + 1L  # OQ3 with no OQ2 drop
    if (length(oq3Examples) < 4) {
      oq3Examples[[length(oq3Examples) + 1L]] <- list(
        n = n, k = k, over = over, hadDrop = length(dropped) > 0,
        trees = vapply(trees, write.tree, character(1)))
    }
  }
}

cat("====  Clade-rule diagnosis (both directions)  ====\n")
cat(sprintf("  trials: %d\n\n", nTrial))
cat("-- OQ2 (under-resolution): dropped majority clades\n")
cat(sprintf("   dropped clades:               %d\n", oq2Dropped))
cat(sprintf("   ... satisfy clade rule:       %d  (=> lost only to flat-collapse)\n", oq2RuleHolds))
cat(sprintf("   ... violate clade rule:       %d\n\n", oq2RuleFails))
cat("-- OQ3 (over-resolution): produced clades violating the strict rule\n")
cat(sprintf("   produced clades (total):      %d\n", producedTotal))
cat(sprintf("   ... violate clade rule:       %d  (over-resolved vs strict Degnan)\n", oq3Clades))
cat(sprintf("   trials with >=1 over-res:     %d / %d\n", oq3Trials, nTrial))
cat(sprintf("   ... of which had NO OQ2 drop: %d  (=> OQ3 stands alone)\n", oq3Standalone))
if (length(oq3Examples)) {
  cat("\n  --- OQ3 examples ---\n")
  for (ex in oq3Examples) {
    cat(sprintf("  n=%d k=%d  over-resolved: %s   (OQ2 drop also? %s)\n",
                ex$n, ex$k, paste(ex$over, collapse = " ; "), ex$hadDrop))
  }
}
