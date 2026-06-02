# Adversarial battery for FreqDiff scratch-bounds (3n/2n) and radix (5n/k) caps.
# Run against an _GLIBCXX_ASSERTIONS build: an OOB vector[] aborts the process;
# a radix cap hit throws "radix sort n/k overflow" via Rcpp::stop.
.libPaths(c('C:/Users/pjjg18/GitHub/Consensus/.agent-cons', .libPaths()))
suppressMessages({library(ConsTree); library(TreeTools); library(ape)})

set.seed(1)

caterpillar <- function(n) {
  # pectinate / fully-asymmetric tree: short heavy path won't exist; instead
  # build a "broom": one big balanced clade hanging off a short backbone.
  ape::stree(n, type = "left", tip.label = paste0("t", seq_len(n)))
}

broom <- function(n) {
  # backbone of length 2, one giant star-ish subtree -> long centroid path side
  t <- ape::rtree(n, tip.label = paste0("t", seq_len(n)))
  t
}

run_case <- function(label, trees) {
  trees <- lapply(trees, function(x) {x$edge.length <- NULL; x})
  class(trees) <- "multiPhylo"
  res <- tryCatch({
    fr <- Frequency(trees)
    sprintf("OK  splits=%d", length(as.Splits(fr)))
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("[%s] %s\n", label, res))
  invisible(NULL)
}

# Sweep n near 3n boundary triggers; incongruent ensembles maximise special nodes
for (n in c(4,5,6,7,8,16,31,32,33,63,64,65,100,127,128,129,200,255,256,257)) {
  k <- 8
  trees <- replicate(k, RandomTree(n, root = FALSE), simplify = FALSE)
  run_case(sprintf("rand n=%d k=%d", n, k), trees)
}

# k >> n : exercise weight-compression branch (radix->k > radix->n)
for (n in c(5, 8, 16, 32)) {
  k <- 50
  trees <- replicate(k, RandomTree(n, root = FALSE), simplify = FALSE)
  run_case(sprintf("k>>n n=%d k=%d", n, k), trees)
}

# Pathological shapes: caterpillar vs balanced (max incongruence -> deep contraction)
for (n in c(8, 16, 32, 64, 128)) {
  t1 <- ape::stree(n, type = "left",  tip.label = paste0("t", seq_len(n)))
  t2 <- ape::stree(n, type = "right", tip.label = paste0("t", seq_len(n)))
  t3 <- as.phylo(ape::stree(n, type = "balanced", tip.label = paste0("t", seq_len(n))))
  # rotate labels to force conflict
  t4 <- ape::stree(n, type = "left", tip.label = paste0("t", c(seq(2,n),1)))
  run_case(sprintf("cat/bal n=%d", n), list(t1,t2,t3,t4))
}

# star trees (single internal node) and near-star
for (n in c(8, 32, 64)) {
  st <- ape::stree(n, tip.label = paste0("t", seq_len(n)))  # star
  rt <- RandomTree(n, root = FALSE)
  run_case(sprintf("star+rand n=%d", n), list(st, rt, st, rt))
}

cat("BATTERY COMPLETE\n")
