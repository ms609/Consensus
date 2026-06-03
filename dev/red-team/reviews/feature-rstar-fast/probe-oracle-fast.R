# Faster oracle: n=8..10 only (2^10 feasible), many trials, exercises both
# verification branches. Plus a targeted big-cluster case to force inside-branch.
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/rstar/.agent-cons", .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })
cladeSet <- function(tree) {
  tree <- TreeTools::Preorder(tree); edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]; paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1)); sizes <- lengths(strsplit(cs, ",")); sort(unique(cs[sizes > 1 & sizes < nTip]))
}
tripletClose <- function(tree) {
  tree2 <- tree; tree2[["edge.length"]] <- rep(1, nrow(tree[["edge"]]))
  dep <- ape::node.depth.edgelength(tree2); m <- ape::mrca(tree)
  function(x, y, z) { dxy <- dep[m[x,y]]; dxz <- dep[m[x,z]]; dyz <- dep[m[y,z]]
    if (dxy==dxz && dxy==dyz) return(NA_character_)
    if (dxy>dxz && dxy>dyz) return(paste(sort(c(x,y)),collapse=","))
    if (dxz>dxy && dxz>dyz) return(paste(sort(c(x,z)),collapse=","))
    paste(sort(c(y,z)),collapse=",") }
}
favoured <- function(st, x, y, Z) { v <- vapply(st, function(s) s(x,y,Z), ""); v <- v[!is.na(v)]
  if (!length(v)) return(FALSE); tb <- table(v); tgt <- paste(sort(c(x,y)),collapse=",")
  if (is.na(tb[tgt])) return(FALSE); mx <- max(tb); tb[tgt]==mx && sum(tb==mx)==1L }
isStrong <- function(st, A, all) { out <- setdiff(all, A); if (!length(out)) return(TRUE)
  pr <- utils::combn(A,2); for (p in seq_len(ncol(pr))) for (Z in out)
    if (!favoured(st, pr[1,p], pr[2,p], Z)) return(FALSE); TRUE }
strongClusters <- function(trees) { all <- trees[[1]]$tip.label; n <- length(all)
  st <- lapply(trees, tripletClose); out <- character(0)
  for (mask in seq_len(2L^n-1L)) { idx <- which(as.integer(intToBits(mask))[seq_len(n)]==1L)
    if (length(idx)<2L||length(idx)>n-1L) next; A <- all[idx]
    if (isStrong(st, A, all)) out <- c(out, paste(sort(A),collapse=",")) }; sort(unique(out)) }
alignTrees <- function(trees) { labs <- trees[[1]]$tip.label; lapply(trees, function(tr) TreeTools::RenumberTips(tr, labs)) }

set.seed(42)
fails <- 0L; ntest <- 0L
for (trial in seq_len(250)) {
  n <- sample(8:10, 1); k <- sample(2:9, 1)
  regime <- sample(c("indep","perturbed","partly"), 1)
  if (regime=="perturbed") { base <- ape::rtree(n,rooted=TRUE); base$edge.length <- NULL
    trees <- lapply(seq_len(k), function(i){ tr <- base; for (s in 1:2){ij<-sample.int(n,2);tr$tip.label[ij]<-tr$tip.label[rev(ij)]};tr})
  } else { trees <- lapply(seq_len(k), function(i){ tr <- ape::rtree(n,rooted=TRUE)
    if (regime=="partly") tr <- ape::di2multi(tr, tol=0.3); tr$edge.length <- NULL; tr}) }
  trees <- alignTrees(trees); ntest <- ntest + 1L
  got <- cladeSet(RStar(trees)); exp <- strongClusters(trees)
  if (!setequal(got, exp)) { fails <- fails+1L; cat("MISMATCH",trial,"n",n,"k",k,regime,"\n")
    cat("  extra:  ", paste(setdiff(got,exp),collapse=" | "),"\n")
    cat("  missing:", paste(setdiff(exp,got),collapse=" | "),"\n"); if (fails>=8) break }
}
cat(sprintf("\nFAST oracle: %d / %d exact (n=8..10).\n", ntest-fails, ntest))
