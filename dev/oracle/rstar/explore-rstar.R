# explore-rstar.R
# R* validation exploration (dev-only).  Three pillars, because R* has NO direct
# reference binary:
#   A. Identity:   RStar(k identical copies of binary T) must equal T exactly.
#   B. Oracle:     on CONGRUENT input (all trees agree on every resolved triplet)
#                  the plurality set == the unanimous set, so RStar must equal the
#                  FDCT `aho-build` BUILD oracle.  Validates the assembly half.
#   C. Refinement: on CONFLICTING input, R* clades must contain (i) all strict
#                  (unanimous) rooted clades [must always hold] and (ii) all
#                  rooted MAJORITY clades [the check that bites on conflict; a
#                  failure is a FINDING about the OQ2 flat-collapse default,
#                  per the design notes — not a test to weaken].
#
# Run from the package root:
#   Rscript.exe dev/oracle/rstar/explore-rstar.R

.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))
suppressMessages(library(ape))
source("dev/oracle/oracle.R")          # CladeSet()
source("dev/oracle/local/oracle_local.R")  # LocalOracle(., "aho")

multi2list <- function(mp) lapply(seq_along(mp), function(i) mp[[i]])

# Rooted clades present in strictly more than `frac` of the trees.
clTally <- function(trees) {
  k <- length(trees)
  tab <- table(unlist(lapply(trees, CladeSet)))
  list(k = k, tab = tab)
}
majClades    <- function(trees) { t <- clTally(trees); names(t$tab)[t$tab >  t$k / 2] }
strictClades <- function(trees) { t <- clTally(trees); names(t$tab)[t$tab == t$k] }

# Align a set of rtree()s onto a shared label order (R* needs identical leaf sets)
alignTrees <- function(trees) {
  labs <- trees[[1]]$tip.label
  lapply(trees, function(tr) RenumberTips(tr, labs))
}

# =====================================================================
# A. Identity: RStar(k copies of binary T) == T
# =====================================================================
cat("====  A. Identity (RStar of k identical binary trees == T)  ====\n")
set.seed(101)
idFail <- 0L
for (trial in 1:60) {
  n <- sample(4:40, 1)
  k <- sample(2:6, 1)
  T <- rtree(n, rooted = TRUE); T$edge.length <- NULL
  r <- RStar(rep(list(T), k))
  if (!setequal(CladeSet(r), CladeSet(T))) {
    idFail <- idFail + 1L
    if (idFail <= 3) cat(sprintf("  FAIL n=%d k=%d: %s\n", n, k, write.tree(T)))
  }
}
cat(sprintf("  identity failures: %d / 60\n\n", idFail))

# =====================================================================
# B. Oracle cross-check on CONGRUENT input (vs FDCT aho-build)
# =====================================================================
cat("====  B. Congruent oracle: RStar == LocalOracle(.,'aho')  ====\n")
cat(sprintf("  %-12s  %s\n", "Fixture", "Result"))
congruent <- list(
  ident_n4  = multi2list(ape::as.phylo(rep(0, 3), 4)),
  ident_n6  = multi2list(ape::as.phylo(rep(3, 3), 6)),
  ident_n8  = multi2list(ape::as.phylo(rep(11, 4), 8)),
  ident_n10 = multi2list(ape::as.phylo(rep(7, 5), 10)),
  ident_n12 = multi2list(ape::as.phylo(rep(99, 3), 12)),
  ident_n16 = multi2list(ape::as.phylo(rep(40, 3), 16))
)
oFail <- 0L
for (nm in names(congruent)) {
  trees <- congruent[[nm]]
  mine <- tryCatch(CladeSet(RStar(trees)), error = function(e) NULL)
  orac <- tryCatch({ o <- LocalOracle(trees, "aho"); if (is.null(o)) NULL else CladeSet(o) },
                   error = function(e) NULL)
  if (is.null(orac)) { cat(sprintf("  %-12s  SKIP (oracle NULL)\n", nm)); next }
  ok <- setequal(mine, orac)
  if (!ok) oFail <- oFail + 1L
  cat(sprintf("  %-12s  %s\n", nm, if (ok) "MATCH" else "DIFFER"))
  if (!ok) {
    cat("     extra:  ", paste(setdiff(mine, orac), collapse = "; "), "\n")
    cat("     missing:", paste(setdiff(orac, mine), collapse = "; "), "\n")
  }
}
cat(sprintf("  oracle DIFFERs: %d / %d\n\n", oFail, length(congruent)))

# =====================================================================
# C. Refinement property on CONFLICTING random input
# =====================================================================
cat("====  C. Refinement on random conflicting input  ====\n")
set.seed(202)
nTrial <- 400L
strictViol <- 0L; majViol <- 0L
majViolCases <- list()
for (trial in seq_len(nTrial)) {
  n <- sample(5:12, 1)
  k <- sample(3:9, 1)
  trees <- alignTrees(lapply(seq_len(k), function(i) {
    tr <- rtree(n, rooted = TRUE); tr$edge.length <- NULL; tr
  }))
  rs <- tryCatch(CladeSet(RStar(trees)), error = function(e) NULL)
  if (is.null(rs)) next
  sc <- strictClades(trees)
  mc <- majClades(trees)
  if (!all(sc %in% rs)) strictViol <- strictViol + 1L
  if (!all(mc %in% rs)) {
    majViol <- majViol + 1L
    if (length(majViolCases) < 4) {
      majViolCases[[length(majViolCases) + 1L]] <- list(
        n = n, k = k,
        missing = setdiff(mc, rs),
        trees = vapply(trees, write.tree, character(1))
      )
    }
  }
}
cat(sprintf("  trials: %d\n", nTrial))
cat(sprintf("  STRICT-clade refinement violations:   %d  (MUST be 0)\n", strictViol))
cat(sprintf("  MAJORITY-clade refinement violations: %d  (>0 => OQ2 finding)\n", majViol))
if (length(majViolCases)) {
  cat("\n  --- example majority-refinement violations (candidate findings) ---\n")
  for (cs in majViolCases) {
    cat(sprintf("  n=%d k=%d  missing clades: %s\n", cs$n, cs$k,
                paste(cs$missing, collapse = " ; ")))
    for (w in cs$trees) cat("      ", w, "\n")
  }
}
cat("\n====  done  ====\n")
