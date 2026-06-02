# Machinery shared by the split-selection consensus methods (`Loose()`,
# `Greedy()`, `MajorityPlus()`).  Each pools the bipartitions ("splits") that
# occur across the input trees, then retains a subset according to a
# method-specific rule before rebuilding a tree.  Selection operates on the
# distinct splits and their occurrence counts, exactly the substrate produced
# by `TreeTools::Consensus()`; reconstruction and rooting reuse `TreeTools`,
# keeping behaviour consistent with `Strict()` and `Majority()`.

# Validate `trees`, drop NULLs, and short-circuit the degenerate cases (a bare
# `phylo`, fewer than two trees, or fewer than four leaves) by returning a
# `trivial` tree.  Shared by `.PoolSplits()` (the R split-selection pipeline) and
# the C++ fast paths (e.g. `Greedy()`), which need the cleaned tree list and the
# shared leaf labels.
#' @importFrom TreeTools TipLabels
.PrepareTrees <- function(trees) {
  if (inherits(trees, "phylo")) {
    return(list(trivial = trees))
  }
  if (!is.list(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  trees <- c(trees)
  trees <- trees[!vapply(trees, is.null, logical(1))]
  nTree <- length(trees)
  if (nTree < 2L) {
    return(list(trivial = if (nTree) trees[[1]] else NULL))
  }
  labels <- TipLabels(trees[[1]])
  if (length(labels) < 4L) {
    return(list(trivial = trees[[1]]))
  }
  # Return:
  list(trivial = NULL, trees = trees, labels = labels, nTree = nTree,
       firstTree = trees[[1]])
}

# Edge matrices for the FACT split-method C++ ports.  Each tree is renumbered to
# the shared `labels` and rooted at taxon 1 (`labels[[1]]`), so the rooted
# clusters the C++ extracts correspond consistently to unrooted bipartitions --
# the convention FACT's reference C++ uses (`tree.cpp`: an unrooted tree is
# rooted at the node adjacent to taxon 1).  Returned in `Preorder` (parent before
# child), as `buildTreeFromEdge()` requires.
#' @importFrom TreeTools Preorder RenumberTips RootTree
.FactEdges <- function(trees, labels) {
  lapply(trees, function(tr) {
    tr <- RenumberTips(tr, labels)
    tr <- RootTree(tr, labels[[1]])
    Preorder(tr)[["edge"]]
  })
}

# Pool the splits occurring across `trees`, returning their distinct values and
# occurrence counts -- or, when no computation is required, a `trivial` tree.
#' @importFrom TreeTools as.Splits
.PoolSplits <- function(trees) {
  prep <- .PrepareTrees(trees)
  if (!is.null(prep[["trivial"]])) {
    return(list(trivial = prep[["trivial"]]))
  }
  trees <- prep[["trees"]]
  labels <- prep[["labels"]]
  splitList <- lapply(trees, as.Splits, tipLabels = labels)
  pooled <- do.call(c, splitList)
  distinct <- unique(pooled)
  distinctKeys <- as.character(distinct)
  counts <- tabulate(match(as.character(pooled), distinctKeys),
                     nbins = length(distinct))
  # For each tree, the indices of the distinct splits it displays:
  membership <- lapply(splitList,
                       function(s) match(as.character(s), distinctKeys))
  # Return:
  list(trivial = NULL,
       splits = distinct,
       members = as.logical(distinct),
       counts = counts,
       membership = membership,
       labels = labels,
       nTree = prep[["nTree"]],
       firstTree = trees[[1]])
}

# Pairwise compatibility matrix among the distinct splits in `prep`.
#' @importFrom TreeTools CompatibleSplits
.CompatibilityMatrix <- function(prep) {
  # Return:
  as.matrix(CompatibleSplits(prep[["splits"]], prep[["splits"]]))
}

# Root `tree` to match the first input tree, as `TreeTools::Consensus()` does:
# place the root on the split that separates the first tree's root group from the
# rest.  Factored out of `.SelectedConsensus()` so the C++ fast paths (`Greedy()`
# and the Loose/MajorityPlus chips) reuse identical rooting.
#' @importFrom TreeTools DescendantEdges NTip Preorder RootTree
.RootLikeFirst <- function(tree, firstTree) {
  first <- Preorder(firstTree)
  edge <- first[["edge"]]
  rootTips <- edge[DescendantEdges(edge[, 1], edge[, 2], edge = 1), 2]
  rootTips <- rootTips[rootTips <= NTip(first)]
  # Return:
  RootTree(tree, first[["tip.label"]][rootTips])
}

# Rebuild a tree from the splits selected by `keep` (a logical vector over
# `prep$splits`), rooted to match the first input tree as `TreeTools::Consensus`
# does.
#' @importFrom ape as.phylo
#' @importFrom TreeTools as.Splits
.SelectedConsensus <- function(keep, prep) {
  selected <- as.Splits(prep[["members"]][keep, , drop = FALSE],
                        tipLabels = prep[["labels"]])
  # Return:
  .RootLikeFirst(as.phylo(selected), prep[["firstTree"]])
}

#' Loose consensus tree
#'
#' `Loose()` returns the loose consensus, also known as the semi-strict or
#' combinable-component consensus \insertCite{Bremer1990}{ConsTree}.  It
#' contains every split that is contradicted by none of the input trees;
#' equivalently, every split that is compatible with each input tree.
#'
#' The loose consensus is always at least as resolved as the strict consensus
#' ([`Strict()`]) and never includes a grouping that conflicts with any input
#' tree.  It is incomparable with the majority-rule consensus ([`Majority()`]):
#' a split present in most trees may still be contradicted by a minority, and so
#' excluded from the loose consensus, whereas a split occurring in a single tree
#' is retained if no other tree contradicts it.
#'
#' This implementation ports the asymptotically efficient `looseConsensusFast`
#' algorithm of \insertCite{JanssonShenSung2016}{ConsTree} from their FACT
#' toolkit (used with permission): the input trees are merged into a one-way
#' compatible tree by repeated linear-time consecutive-range queries, the
#' clusters that are compatible with every input are then marked, and the rest
#' contracted away -- avoiding the explicit pairwise compatibility matrix used
#' previously.
#'
#' @inheritParams Strict
#'
#' @return `Loose()` returns the consensus tree, an object of class `phylo`,
#' rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' Loose(trees)
#'
#' @seealso Closely related: [`Strict()`], [`Majority()`], [`Greedy()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape read.tree
#' @export
Loose <- function(trees) {
  prep <- .PrepareTrees(trees)
  if (!is.null(prep[["trivial"]])) {
    return(prep[["trivial"]])
  }
  labels <- prep[["labels"]]
  nwk <- looseConsensusCpp(.FactEdges(prep[["trees"]], labels), length(labels))
  tree <- read.tree(text = paste0(nwk, ";"))
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return:
  .RootLikeFirst(tree, prep[["firstTree"]])
}

#' Greedy (extended majority-rule) consensus tree
#'
#' `Greedy()` returns the greedy consensus, also termed the extended
#' majority-rule consensus \insertCite{Bryant2003}{ConsTree}.  Distinct splits
#' are considered in decreasing order of their frequency across the input trees;
#' each is added to the growing consensus if it is compatible with every split
#' already accepted.  The result is typically more resolved than the
#' majority-rule consensus ([`Majority()`]), and contains it.
#'
#' Splits that occur equally often are considered in a fixed, reproducible order,
#' so the result is deterministic.  Where several mutually incompatible splits
#' are equally frequent, a different (but equally valid) greedy resolution may be
#' returned by other software.
#'
#' This implementation ports the asymptotically efficient `greedyConsensusFast`
#' algorithm of \insertCite{JanssonShenSung2016}{ConsTree} from their FACT
#' toolkit (used with permission): the distinct clusters are extracted in a
#' single post-order sweep of each tree and added in decreasing order of
#' frequency whenever compatible with the tree built so far, avoiding the
#' explicit pairwise compatibility matrix used previously.
#'
#' @inheritParams Strict
#'
#' @return `Greedy()` returns the consensus tree, an object of class `phylo`,
#' rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' Greedy(trees)
#'
#' @seealso Closely related: [`Strict()`], [`Majority()`], [`Loose()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape read.tree
#' @export
Greedy <- function(trees) {
  prep <- .PrepareTrees(trees)
  if (!is.null(prep[["trivial"]])) {
    return(prep[["trivial"]])
  }
  labels <- prep[["labels"]]
  nwk <- greedyConsensusCpp(.FactEdges(prep[["trees"]], labels), length(labels))
  tree <- read.tree(text = paste0(nwk, ";"))
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return:
  .RootLikeFirst(tree, prep[["firstTree"]])
}

#' Majority-rule (+) consensus tree
#'
#' `MajorityPlus()` returns the majority-rule (+) consensus
#' \insertCite{JanssonShenSung2016}{ConsTree}: a clade is retained when it
#' occurs in more input trees than contradict it -- i.e. when the number of
#' trees displaying the clade exceeds the number of trees incompatible with it.
#' A tree that is compatible with a clade without displaying it counts neither
#' for nor against.
#'
#' Every majority-rule split is retained (a split in more than half the trees is
#' contradicted by fewer than half), so `MajorityPlus()` contains the
#' majority-rule consensus ([`Majority()`]) and may add further splits that are
#' supported more often than they are contradicted.  The retained splits are
#' necessarily mutually compatible, so they define a valid tree.
#'
#' \insertCite{JanssonShenSung2016;textual}{ConsTree} give an optimal
#' \eqn{O(kn)} algorithm for this consensus, implemented in their FACT toolkit;
#' here the same tree is computed directly from the pooled splits.
#'
#' @inheritParams Strict
#'
#' @return `MajorityPlus()` returns the consensus tree, an object of class
#' `phylo`, rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' MajorityPlus(trees)
#'
#' @seealso Closely related: [`Majority()`], [`Greedy()`], [`Loose()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @export
MajorityPlus <- function(trees) {
  prep <- .PoolSplits(trees)
  if (!is.null(prep[["trivial"]])) {
    return(prep[["trivial"]])
  }
  conflict <- !.CompatibilityMatrix(prep)
  nSplit <- length(prep[["counts"]])
  # incidence[i, j]: does tree i display distinct split j?
  incidence <- matrix(FALSE, prep[["nTree"]], nSplit)
  for (i in seq_len(prep[["nTree"]])) {
    incidence[i, prep[["membership"]][[i]]] <- TRUE
  }
  # conflictsPerTree[j, i]: number of splits in tree i incompatible with split j
  conflictsPerTree <- conflict %*% t(incidence)
  incompatibleTrees <- rowSums(conflictsPerTree > 0)
  keep <- prep[["counts"]] > incompatibleTrees
  # Return:
  .SelectedConsensus(keep, prep)
}

#' Frequency-difference consensus tree
#'
#' `Frequency()` returns the frequency-difference consensus: a split is retained
#' when it occurs strictly more often than every split that conflicts with it.
#' Equivalently, among each set of mutually incompatible splits, the consensus
#' keeps a split only if it is strictly more frequent than all its rivals.
#'
#' The frequency-difference consensus is at least as resolved as the
#' majority-rule consensus ([`Majority()`]) -- a split in more than half the
#' trees is necessarily more frequent than any conflicting split -- and is
#' contained within the greedy consensus ([`Greedy()`]).  The retained splits
#' are mutually compatible, so they define a valid tree.
#'
#' An efficient algorithm for constructing this consensus was given by
#' \insertCite{Jansson2024}{ConsTree}; the present implementation computes the
#' same tree directly from pooled split frequencies.
#'
#' @inheritParams Strict
#'
#' @return `Frequency()` returns the consensus tree, an object of class `phylo`,
#' rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' Frequency(trees)
#'
#' @seealso Closely related: [`Majority()`], [`MajorityPlus()`],
#' [`Greedy()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @export
Frequency <- function(trees) {
  prep <- .PoolSplits(trees)
  if (!is.null(prep[["trivial"]])) {
    return(prep[["trivial"]])
  }
  conflict <- !.CompatibilityMatrix(prep)
  counts <- prep[["counts"]]
  keep <- vapply(seq_along(counts), function(i) {
    rivals <- conflict[i, ]
    counts[[i]] > if (any(rivals)) max(counts[rivals]) else 0L
  }, logical(1))
  # Return:
  .SelectedConsensus(keep, prep)
}

