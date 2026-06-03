# Independent brute-force Adams, from the DEFINITION (Jansson et al. 2017, p335):
#   pi(S) = product over trees of pi(T_j); two leaves share a block iff in EVERY
#   tree they descend from the same child of lca_{T_j}(current set).
#   Recurse per block until |block| <= 2.  Adams tree is unique.
#
# This implementation derives "which child of the LCA" WITHOUT KeepTip and
# WITHOUT the C++ nub machinery -- purely from ancestor sets on the ORIGINAL
# trees -- so it shares no bug mode with either.
.libPaths(c('C:/Users/pjjg18/GitHub/worktrees/Consensus/adams/.agent-cons', .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))

# Precompute, per tree: parent[], children[], and for each node its tip-set (as a
# logical over taxa 1..n on the SHARED labelling).
prep <- function(tr, labels) {
  tr <- Preorder(RenumberTips(tr, labels))
  edge <- tr[["edge"]]
  nTip <- length(labels)
  nNode <- max(edge)
  parent <- integer(nNode); parent[] <- NA_integer_
  kids <- vector("list", nNode)
  for (i in seq_len(nrow(edge))) {
    p <- edge[i, 1]; c <- edge[i, 2]
    parent[c] <- p
    kids[[p]] <- c(kids[[p]], c)
  }
  root <- edge[1, 1]
  # tipset[node] = integer vector of taxa under node
  tipset <- vector("list", nNode)
  # fill via postorder (children before parent): process edges in reverse preorder
  for (v in seq_len(nTip)) tipset[[v]] <- v
  ord <- rev(unique(c(edge[, 1], edge[, 2])))  # children appear after parents in preorder
  for (v in ord) {
    if (v > nTip) {
      tipset[[v]] <- sort(unlist(lapply(kids[[v]], function(c) tipset[[c]])))
    }
  }
  # ancestor lookup: for a taxon leaf, the path of nodes root..leaf
  pathTo <- function(taxon) {
    p <- taxon; chain <- p
    while (!is.na(parent[p])) { p <- parent[p]; chain <- c(chain, p) }
    rev(chain)  # root first
  }
  list(parent = parent, kids = kids, root = root, tipset = tipset,
       nTip = nTip, pathTo = pathTo)
}

# LCA of a set of taxa in a prepped tree = deepest node whose tipset is a superset.
lcaNode <- function(P, taxa) {
  # walk down from root: among children, find the unique child whose tipset
  # contains ALL taxa; if none, current node is the LCA.
  v <- P$root
  repeat {
    nextv <- NA_integer_
    for (c in P$kids[[v]]) {
      if (all(taxa %in% P$tipset[[c]])) { nextv <- c; break }
    }
    if (is.na(nextv)) return(v)
    v <- nextv
  }
}

# "which child of lca(taxa) does `t` descend from" -- returns that child node id
# (or the leaf itself if the leaf IS the lca, impossible for |taxa|>=2 unless t
# not under lca).
childUnderLCA <- function(P, lca, t) {
  for (c in P$kids[[lca]]) if (t %in% P$tipset[[c]]) return(c)
  NA_integer_  # t not under lca -- shouldn't happen for taxa within the block
}

bruteAdamsClades <- function(trees) {
  labels <- TipLabels(trees[[1]])
  Ps <- lapply(trees, prep, labels = labels)
  clades <- character(0)
  rec <- function(taxa) {
    if (length(taxa) <= 2L) return(invisible())
    # signature per taxon = (child-of-lca in tree1, ..., in tree k)
    lcas <- lapply(Ps, lcaNode, taxa = taxa)
    sig <- vapply(taxa, function(t) {
      paste(vapply(seq_along(Ps), function(j)
        as.character(childUnderLCA(Ps[[j]], lcas[[j]], t)), character(1)),
        collapse = "|")
    }, character(1))
    blocks <- split(taxa, sig)
    for (B in blocks) {
      if (length(B) >= 2L)
        clades[[length(clades) + 1L]] <<- paste(sort(labels[B]), collapse = ",")
      rec(B)
    }
  }
  rec(seq_len(length(labels)))
  sort(unique(clades))
}

# cladeSet of an Adams output (label sets)
cladeSetLab <- function(tree) {
  tree <- Preorder(tree)
  edge <- tree[["edge"]]; nTip <- NTip(tree); tl <- tree[["tip.label"]]
  cl <- vapply(seq_len(nrow(edge)), function(e) {
    below <- DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]
    paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sz <- lengths(strsplit(cl, ","))
  sort(unique(cl[sz > 1 & sz < nTip]))
}
