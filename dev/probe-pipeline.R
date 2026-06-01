# Probe v3: validate the all-R selection pipeline by reproducing Majority/Strict
# from scratch and comparing to TreeTools::Consensus; sanity-check Loose/Greedy.
.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(TreeTools))

RF <- TreeDist::RobinsonFoulds

for (spec in list(0:5, 0:20, c(0, 1, 2, 53, 99, 4))) {
  trees <- ape::as.phylo(spec, 9)
  k <- length(trees)
  labels <- TipLabels(trees[[1]])
  pooled <- do.call(c, lapply(trees, as.Splits, tipLabels = labels))

  u <- unique(pooled)
  counts <- tabulate(match(as.character(pooled), as.character(u)),
                     nbins = length(u))
  uMat <- as.logical(u)
  mkTree <- function(rows) {
    as.phylo(as.Splits(uMat[rows, , drop = FALSE], tipLabels = labels))
  }

  thresh <- as.integer(k * 0.5) + 1L
  treeMaj <- mkTree(counts >= thresh)
  treeStrict <- mkTree(counts == k)

  cmat <- as.matrix(CompatibleSplits(u, u))
  compatAll <- apply(cmat, 1, all)            # loose
  treeLoose <- mkTree(compatAll)

  # greedy: count-descending, accept if compatible with all accepted so far
  ord <- order(counts, decreasing = TRUE)
  accepted <- integer(0)
  for (i in ord) {
    if (length(accepted) == 0 || all(cmat[i, accepted])) accepted <- c(accepted, i)
  }
  treeGreedy <- mkTree(seq_len(nrow(uMat)) %in% accepted)

  cat(sprintf(
    "k=%2d distinct=%2d | MAJ RF=%g  STRICT RF=%g | NSplits  strict=%d loose=%d maj=%d greedy=%d\n",
    k, length(u),
    RF(treeMaj, Consensus(trees, p = 0.5)),
    RF(treeStrict, Consensus(trees, p = 1)),
    NSplits(treeStrict), NSplits(treeLoose), NSplits(treeMaj), NSplits(treeGreedy)))
}
