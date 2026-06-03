# Correctness signal in the now-UNCAPPED n > 200 regime, where the brute-force
# oracle (n<=12) and the new-vs-legacy diff (n<=200) cannot reach.  We cannot
# brute-force strong clusters at n=300, but two NECESSARY conditions are cheap:
#   (1) R* refines the strict consensus  (every unanimous clade survives), and
#   (2) R* refines the majority-rule consensus (Lemma: R* >= Majority).
# Both use an INDEPENDENT reference (TreeTools::Consensus), different code path.
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/rstar/.agent-cons", .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })

cladeVec <- function(tree) {
  tree <- TreeTools::Preorder(tree); edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]; paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sizes <- lengths(strsplit(cs, ",")); sort(unique(cs[sizes > 1 & sizes < nTip]))
}

set.seed(123)
ok <- TRUE
for (n in c(250L, 300L, 500L)) {
  k <- 12L
  base <- RandomTree(paste0("t", seq_len(n)), root = TRUE)
  trees <- lapply(seq_len(k), function(i) {
    tr <- base
    for (s in 1:4) { ij <- sample.int(n, 2); tr$tip.label[ij] <- tr$tip.label[rev(ij)] }
    RenumberTips(tr, base$tip.label)
  })
  t0 <- proc.time()[["elapsed"]]
  rs  <- cladeVec(RStar(trees))
  el  <- proc.time()[["elapsed"]] - t0
  maj <- cladeVec(TreeTools::Consensus(trees, p = 0.5))
  str <- cladeVec(TreeTools::Consensus(trees, p = 1))
  majOK <- all(maj %in% rs); strOK <- all(str %in% rs)
  ok <- ok && majOK && strOK
  cat(sprintf("n=%-4d k=%d  RStar=%4.2fs  clades: RStar=%d maj=%d strict=%d | maj<=RStar:%s strict<=RStar:%s\n",
              n, k, el, length(rs), length(maj), length(str), majOK, strOK))
  if (!majOK) cat("   MAJORITY clades missing from RStar:", paste(setdiff(maj, rs), collapse=" | "), "\n")
  if (!strOK) cat("   STRICT clades missing from RStar:",   paste(setdiff(str, rs), collapse=" | "), "\n")
}
cat(if (ok) "\nPASS: R* refines strict & majority at n=250,300,500.\n"
    else    "\nFAIL: a necessary refinement condition violated.\n")
