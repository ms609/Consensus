# Adams consensus.  Unlike the split-selection family in `selection.R`, the Adams
# consensus is a *rooted* method and can contain clusters present in no input
# tree, so it is not built by selecting pooled splits.  Instead it follows
# Adams' (1972) recursive definition: partition the current taxon set by the
# top-level (root-children) grouping shared across every input tree, then recurse
# within each block.  The heavy lifting is done in C++ (src/cons_adams.cpp) by the
# O(kn log n) centroid-path algorithm of Jansson, Li & Sung (2017); this wrapper
# only marshals each tree (on its OWN root) to integer-coded preorder edges and
# relabels the integer-labelled Newick it returns.

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
#' @details
#' This implementation employs the asymptotically efficient \eqn{O(kn \log n)}
#' algorithm of \insertCite{JanssonLiSung2017}{ConsTree} (their
#' `New_Adams_consensus_k`, from the FACT toolkit of Jansson and colleagues, used
#' with permission), realised in C++.  Rather than recomputing Adams' partition
#' from scratch at every recursion level, it follows a centroid path through each
#' input tree in unison, expanding the shared "spine" of the consensus
#' iteratively and recursing only on the off-spine blocks (each at most half the
#' leaves).  The Adams consensus tree is unique, so the result is identical to the
#' classical recursive definition; only the running time differs.
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
  # Marshal each tree on its OWN root (Adams is rooted; do not re-root at taxon 1
  # as the split methods do).  RenumberTips aligns ape tip i with labels[i] in
  # every tree, so the integer labels the C++ returns map straight back.
  edgeList <- lapply(trees, function(tr) {
    Preorder(RenumberTips(tr, labels))[["edge"]]
  })
  nwk <- paste0(adamsConsensusCpp(edgeList, n), ";")
  tree <- read.tree(text = nwk)
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return: rooted by construction -- do not re-root.
  tree
}
