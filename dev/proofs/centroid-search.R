# Adversarial search for a "centroid" in QuartetConsensus output
# =============================================================================
# Question (the PI's "can QC ever centroid?"): can QuartetConsensus resolve a
# split (a relationship) that is present in NO input tree -- i.e. place a
# grouping that every input tree contradicts?
#
# This script searches hard for such a counterexample across two regimes:
#   (A) single rogue on a FIXED backbone (only the rogue's position varies);
#   (B) genuine topological conflict (random-topology mixtures, large unstable
#       clades resolved randomly per tree, fully random trees).
#
# It reports, per regime, the MINIMUM split-frequency of any non-trivial split
# that QC resolves, and the number of configs in which QC resolved a split of
# frequency 0 (a true centroid).
#
# RESULT (seed 1, as run 2026-05-31, build of branch consensus-rogues):
#   (A) single-rogue : min split-freq = 0.5000 ; centroid configs = 0
#   (B) conflict      : min split-freq = 0.0500 (= 1/20, one tree) ; centroid = 0
# i.e. QC resolves minority splits (down to a single tree) but NEVER a split
# absent from all trees -- consistent with the constructive proof in
# quartet-consensus-no-centroid.md.
#
# METHODOLOGICAL NOTE (a bug that produced 74 false "centroids" before it was
# caught): split frequency must be computed on a CANONICAL representation of the
# bipartition.  Comparing only the "smaller side" fails for balanced (n/2 | n/2)
# splits, because the smaller-side tiebreak is inconsistent between trees, so a
# present split looks absent (frequency spuriously 0).  Canonicalise to the side
# containing the lexicographically-first tip.
# =============================================================================
suppressMessages({library(TreeTools); library(Consensus); library(ape)})

graft <- function(bb, tip, lab = "R") {
  ape::read.tree(text = sub(tip, paste0("(", tip, ",", lab, ")"),
                            ape::write.tree(bb), fixed = TRUE))
}

# Canonical bipartition keys of a tree: each non-trivial split as the sorted
# label-set of the side containing the lexicographically-first tip (tie-safe).
canonKeys <- function(tr) {
  tr <- ape::collapse.singles(tr)
  m <- as.logical(as.Splits(tr)); tp <- colnames(m); first <- sort(tp)[1]
  ks <- vapply(seq_len(nrow(m)), function(i) {
    s <- m[i, ]; side <- tp[s]; if (!(first %in% side)) side <- tp[!s]
    paste(sort(side), collapse = ",")
  }, character(1))
  sz <- lengths(strsplit(ks, ","))
  ks[sz >= 2 & sz <= length(tp) - 2]            # non-trivial only
}

# Minimum frequency, among the input trees, of any split QC resolves.
minSplitFreq <- function(ins) {
  qc <- try(ape::collapse.singles(Consensus::QuartetConsensus(ins)), silent = TRUE)
  if (inherits(qc, "try-error")) return(NA_real_)
  qk <- canonKeys(qc); if (!length(qk)) return(1)
  inKeys <- lapply(ins, canonKeys)
  min(vapply(qk, function(key) mean(vapply(inKeys, function(K) key %in% K, logical(1))),
             numeric(1)))
}

runRegime <- function(N, generator, label) {
  set.seed(1)
  worst <- 1; nZero <- 0
  for (trial in seq_len(N)) {
    ins <- generator()
    mn <- minSplitFreq(ins); if (is.na(mn)) next
    worst <- min(worst, mn)
    if (mn < 1e-9) nZero <- nZero + 1
  }
  cat(sprintf("%-14s trials=%4d | min QC split-freq = %.4f | centroid (freq-0) configs = %d\n",
              label, N, worst, nZero))
}

# ---- Regime A: single rogue on a fixed backbone -----------------------------
genRogue <- function() {
  n <- sample(6:10, 1); k <- sample(c(12, 24, 40), 1); labs <- paste0("t", seq_len(n))
  bb <- switch(sample(3, 1), BalancedTree(labs), PectinateTree(labs),
               ape::rtree(n, rooted = FALSE, tip.label = sample(labs), br = NULL))
  gen <- sample(4, 1)
  pos <- if (gen == 1) sample(labs, k, TRUE)
    else if (gen == 2) sample(labs, k, TRUE, prob = runif(n)^3)
    else if (gen == 3) sample(labs[1:min(3, n)], k, TRUE)
    else sample(labs[c(1, n)], k, TRUE)
  structure(lapply(pos, function(p) graft(bb, p)), class = "multiPhylo")
}

# ---- Regime B: genuine topological conflict ---------------------------------
genConflict <- function() {
  n <- sample(8:12, 1); k <- sample(c(20, 40), 1); labs <- paste0("t", seq_len(n))
  g <- sample(3, 1)
  if (g == 1) {                                  # mixture of m distinct topologies
    m <- sample(2:4, 1)
    pool <- lapply(seq_len(m), function(.) ape::rtree(n, rooted = FALSE, tip.label = labs, br = NULL))
    structure(lapply(seq_len(k), function(.) pool[[sample(m, 1)]]), class = "multiPhylo")
  } else if (g == 2) {                           # stable backbone + unstable clade
    nu <- sample(5:7, 1); unst <- labs[seq_len(nu)]; stab <- labs[(nu + 1):n]
    base <- paste0("(", paste(stab, collapse = ","), ");")
    structure(lapply(seq_len(k), function(.) {
      ut <- ape::rtree(nu, rooted = TRUE, tip.label = sample(unst), br = NULL)
      uw <- sub(";", "", ape::write.tree(ut))
      ape::read.tree(text = sub(stab[1], paste0("(", stab[1], ",", uw, ")"), base, fixed = TRUE))
    }), class = "multiPhylo")
  } else {                                       # fully random trees (max conflict)
    structure(lapply(seq_len(k), function(.) ape::rtree(n, rooted = FALSE, tip.label = labs, br = NULL)),
              class = "multiPhylo")
  }
}

runRegime(1500, genRogue,    "single-rogue")
runRegime(1200, genConflict, "conflict")
