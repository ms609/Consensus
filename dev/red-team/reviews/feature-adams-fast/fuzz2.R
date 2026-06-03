source('C:/Users/pjjg18/GitHub/worktrees/Consensus/adams/dev/red-team/reviews/feature-adams-fast/brute_adams.R')

# Random MULTIFURCATING rooted tree: start from a random binary rooted tree and
# collapse a fraction of internal edges (di2multi-style) to create polytomies,
# including possibly a high-degree root.
randMultifurc <- function(labs, collapseProb = 0.45) {
  tr <- TreeTools::RandomTree(labs, root = TRUE)
  edge <- tr[["edge"]]
  nTip <- length(labs)
  internalEdges <- which(edge[, 2] > nTip)  # edges leading to internal nodes
  if (length(internalEdges)) {
    drop <- internalEdges[runif(length(internalEdges)) < collapseProb]
    if (length(drop)) {
      # collapse by setting branch length 0 and di2multi
      el <- rep(1, nrow(edge)); el[drop] <- 0
      tr[["edge.length"]] <- el
      tr <- ape::di2multi(tr, tol = 1e-8)
      tr[["edge.length"]] <- NULL
    }
  }
  TreeTools::Preorder(tr)
}

set.seed(7)
fails <- 0; tot <- 0; maxdeg <- 0
for (n in c(5,6,8,10,12,15,18)) {
  for (k in c(2,3,4,6,10)) {
    for (rep in 1:10) {
      labs <- paste0("t", seq_len(n))
      trees <- structure(lapply(seq_len(k), function(i)
        randMultifurc(labs, collapseProb = runif(1, 0.2, 0.7))),
        class = "multiPhylo")
      # record max root degree seen
      for (tt in trees) {
        e <- tt[["edge"]]; rt <- e[1,1]
        maxdeg <- max(maxdeg, sum(e[,1] == rt))
      }
      cpp <- tryCatch(cladeSetLab(Adams(trees)),
                      error = function(e) paste("ERR:", conditionMessage(e)))
      bru <- bruteAdamsClades(trees)
      tot <- tot + 1
      if (!identical(cpp, bru)) {
        fails <- fails + 1
        if (fails <= 8) {
          cat(sprintf("MISMATCH n=%d k=%d rep=%d\n", n, k, rep))
          cat("  cpp  :", paste(cpp, collapse=" ; "), "\n")
          cat("  brute:", paste(bru, collapse=" ; "), "\n")
          for (tt in trees) cat("    ", ape::write.tree(tt), "\n")
        }
      }
    }
  }
}
cat(sprintf("\nMULTIFURCATING ROOTED: %d/%d mismatches; max root degree seen=%d\n",
            fails, tot, maxdeg))
