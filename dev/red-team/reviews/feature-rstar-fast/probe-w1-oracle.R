# NOTE: this is the PRE-FIX diagnostic that proved W1.  Now that RStar() rejects
# unequal leaf sets, every ragged trial errors instead: post-fix this prints
# "190 errored" (not "129 diverge").  Kept as the record of the original defect.
#
# Settle W1 definitively. On ragged input (trees on differing leaf subsets, the
# union = trees[[1]]'s full set), is the C++ result == the natural "each tree
# abstains on triples touching a taxon it lacks" semantics, or does it diverge
# (= silent wrong answer)?  trees[[1]] always carries the full set so the wrapper
# uses the full n; other trees may drop 1-2 leaves (no pre-align; subset => the
# wrapper does NOT error).
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/rstar/.agent-cons", .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })

cladeVec <- function(tree) {
  if (!inherits(tree, "phylo")) return(NA_character_)
  tree <- TreeTools::Preorder(tree); edge <- tree[["edge"]]; nTip <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) { below <- TreeTools::DescendantEdges(edge[,1], edge[,2], edge=e)
    tips <- edge[below,2]; paste(sort(tl[tips[tips <= nTip]]), collapse=",") }, character(1))
  sizes <- lengths(strsplit(cs, ",")); sort(unique(cs[sizes > 1 & sizes < nTip])) }

# Stater that ABSTAINS (NA) on any triple touching a label absent from this tree.
makeStaterRagged <- function(tree) {
  tree2 <- tree; tree2$edge.length <- rep(1, nrow(tree$edge))
  dep <- ape::node.depth.edgelength(tree2); m <- ape::mrca(tree); L <- tree$tip.label
  function(x, y, z) {
    ix <- match(x, L); iy <- match(y, L); iz <- match(z, L); if (anyNA(c(ix, iy, iz))) return(NA_character_)
    dxy <- dep[m[ix,iy]]; dxz <- dep[m[ix,iz]]; dyz <- dep[m[iy,iz]]
    if (dxy==dxz && dxy==dyz) return(NA_character_)
    if (dxy>dxz && dxy>dyz) return(paste(sort(c(x,y)),collapse=","))
    if (dxz>dxy && dxz>dyz) return(paste(sort(c(x,z)),collapse=","))
    paste(sort(c(y,z)),collapse=",") } }
favoured <- function(st, x, y, Z) { v <- vapply(st, function(s) s(x,y,Z), ""); v <- v[!is.na(v)]
  if (!length(v)) return(FALSE); tb <- table(v); tgt <- paste(sort(c(x,y)),collapse=",")
  if (is.na(tb[tgt])) return(FALSE); mx <- max(tb); tb[tgt]==mx && sum(tb==mx)==1L }
isStrong <- function(st, A, all) { out <- setdiff(all, A); if (!length(out)) return(TRUE)
  pr <- utils::combn(A,2); for (p in seq_len(ncol(pr))) for (Z in out) if (!favoured(st, pr[1,p], pr[2,p], Z)) return(FALSE); TRUE }
strongRagged <- function(trees, fullLabs) { n <- length(fullLabs); st <- lapply(trees, makeStaterRagged); out <- character(0)
  for (mask in seq_len(2L^n-1L)) { idx <- which(as.integer(intToBits(mask))[seq_len(n)]==1L)
    if (length(idx)<2L||length(idx)>n-1L) next; A <- fullLabs[idx]
    if (isStrong(st, A, fullLabs)) out <- c(out, paste(sort(A),collapse=",")) }; sort(unique(out)) }

set.seed(7)
fails <- 0L; errs <- 0L; ragged <- 0L; ntest <- 0L
for (trial in seq_len(200)) {
  n <- sample(6:8, 1); k <- sample(3:7, 1); fullLabs <- paste0("t", seq_len(n))
  base <- RandomTree(fullLabs, root = TRUE)
  trees <- lapply(seq_len(k), function(i) { tr <- base
    for (s in 1:2) { ij <- sample.int(n, 2); tr$tip.label[ij] <- tr$tip.label[rev(ij)] }; tr })
  # Drop 1-2 leaves from a few trees (NOT tree 1). This makes the input ragged.
  isRagged <- FALSE
  for (i in 2:k) if (runif(1) < 0.6) {
    d <- sample(fullLabs, sample(1:2, 1)); keep <- setdiff(trees[[i]]$tip.label, d)
    if (length(keep) >= 3) { trees[[i]] <- ape::drop.tip(trees[[i]], d); isRagged <- TRUE } }
  if (!isRagged) next
  ragged <- ragged + 1L; ntest <- ntest + 1L
  got <- tryCatch(cladeVec(RStar(trees)), error = function(e) { errs <<- errs + 1L; structure("ERR", class="e") })
  if (inherits(got, "e")) next
  oracle <- strongRagged(trees, fullLabs)
  if (!setequal(got, oracle)) { fails <- fails + 1L
    if (fails <= 4) { cat(sprintf("DIVERGE trial %d (n=%d k=%d)\n", trial, n, k))
      cat("  RStar :", paste(got, collapse=" | "), "\n  abstain-oracle:", paste(oracle, collapse=" | "), "\n")
      for (tr in trees) cat("   ", write.tree(tr), "\n") } } }
cat(sprintf("\nRagged-input: %d ragged trials | %d diverge from abstain-semantics | %d errored\n",
            ragged, fails, errs))
