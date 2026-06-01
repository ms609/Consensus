# Adams consensus.  Unlike the split-selection family in `selection.R`, the Adams
# consensus is a *rooted* method and can contain clusters present in no input
# tree, so it is not built by selecting pooled splits.  Instead it follows
# Adams' (1972) recursive definition: partition the current taxon set by the
# top-level (root-children) grouping shared across every input tree, then recurse
# within each block.  This mirrors FACT's classical `adamsRecursionSlow`; the
# O(kn log n) variant is deferred.

#' Adams consensus tree
#'
#' `Adams()` returns the Adams consensus \insertCite{Adams1972}{ConsTree}, a
#' summary of a set of **rooted** trees.
#'
#' The Adams consensus is defined recursively.  At the top level, the leaves are
#' partitioned into the coarsest grouping that every input tree refines: two
#' leaves fall in the same block if and only if, in every input tree, they sit
#' below the same child of the most recent common ancestor of the full leaf set.
#' Each block becomes a child of the consensus root, and the construction repeats
#' within each block.
#'
#' Unlike the split-based methods ([`Strict()`], [`Majority()`], [`Loose()`],
#' [`Greedy()`], [`MajorityPlus()`], [`Frequency()`]), the Adams consensus can
#' contain a cluster that appears in no individual input tree, and it depends on
#' how the input trees are rooted: it is a statement about nesting, not about
#' unrooted bipartitions.  Input trees are treated as rooted on their current
#' root; root the trees as you intend before calling `Adams()`.
#'
#' @inheritParams Strict
#'
#' @return `Adams()` returns the consensus tree, an object of class `phylo`,
#' rooted by construction.
#'
#' @examples
#' # Two rooted trees that disagree only on the position of one leaf
#' trees <- c(ape::read.tree(text = "(((a, b), c), d);"),
#'            ape::read.tree(text = "(((a, b), d), c);"))
#' Adams(trees) # keeps the clade (a, b); leaves c, d unresolved at the root
#'
#' @seealso Closely related: [`Strict()`], [`Majority()`], [`Loose()`],
#' [`Greedy()`], [`MajorityPlus()`], [`Frequency()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape read.tree
#' @importFrom TreeTools Preorder RenumberTips TipLabels
#' @export
Adams <- function(trees) {
  if (inherits(trees, "phylo")) {
    return(trees)
  }
  if (!is.list(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  trees <- c(trees)
  trees <- trees[!vapply(trees, is.null, logical(1))]
  nTree <- length(trees)
  if (nTree < 2L) {
    return(if (nTree) trees[[1]] else NULL)
  }
  labels <- TipLabels(trees[[1]])
  n <- length(labels)
  if (n < 3L) {
    return(trees[[1]])
  }
  # Relabel every tree with integer codes 1..n in a shared order, so that the
  # recursion can build a Newick string safely (free of label-escaping issues).
  treesInt <- lapply(trees, function(tr) {
    tr <- RenumberTips(tr, labels)
    tr[["tip.label"]] <- as.character(seq_len(n))
    Preorder(tr)
  })
  nwk <- paste0(.AdamsNewick(as.character(seq_len(n)), treesInt), ";")
  tree <- read.tree(text = nwk)
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return:
  tree
}

# Recursively build the Newick string (without trailing ";") of the Adams
# consensus on the integer-coded leaf set `taxa`, across the relabelled `trees`.
.AdamsNewick <- function(taxa, trees) {
  m <- length(taxa)
  if (m == 1L) {
    return(taxa)
  }
  if (m == 2L) {
    return(paste0("(", taxa[[1]], ",", taxa[[2]], ")"))
  }
  blocks <- .AdamsPartition(taxa, trees)
  if (length(blocks) < 2L) {
    # No tree refines the set: an unresolved (star) node.
    return(paste0("(", paste(taxa, collapse = ","), ")"))
  }
  parts <- vapply(blocks, .AdamsNewick, character(1), trees = trees)
  # Return:
  paste0("(", paste(parts, collapse = ","), ")")
}

# Partition `taxa` into the blocks of the Adams product: leaves grouped so that
# two share a block exactly when, in every tree, they descend from the same child
# of the most recent common ancestor of `taxa`.
#' @importFrom TreeTools DescendantEdges KeepTip NTip Preorder
.AdamsPartition <- function(taxa, trees) {
  sig <- vapply(trees, function(tr) {
    kept <- Preorder(KeepTip(tr, taxa))
    edge <- kept[["edge"]]
    root <- edge[1L, 1L]            # Preorder lists root-emanating edges first
    nt <- NTip(kept)
    rootKids <- which(edge[, 1L] == root)
    blockId <- integer(nt)
    for (b in seq_along(rootKids)) {
      below <- DescendantEdges(edge[, 1L], edge[, 2L], edge = rootKids[[b]])
      childTips <- edge[below, 2L]
      childTips <- childTips[childTips <= nt]
      blockId[childTips] <- b
    }
    # Block of each requested taxon, in the order of `taxa`:
    blockId[match(taxa, kept[["tip.label"]])]
  }, integer(length(taxa)))
  if (is.null(dim(sig))) {
    sig <- matrix(sig, nrow = length(taxa))
  }
  signature <- apply(sig, 1L, paste, collapse = ",")
  # Return:
  unname(split(taxa, signature))
}

