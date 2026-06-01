#' Strict consensus tree
#'
#' `Strict()` returns the strict consensus of a set of trees: the tree that
#' contains exactly those splits (clades) present in _every_ input tree
#' \insertCite{Day1985}{ConsTree}.
#'
#' This is a thin wrapper around [`TreeTools::Consensus()`] with `p = 1`,
#' provided so that every consensus method in this package is reachable through
#' a common, consistently named interface.
#'
#' @param trees A list of trees, or a `multiPhylo` object; all entries must
#' share the same leaf labels.
#'
#' @return `Strict()` returns the consensus tree, an object of class `phylo`,
#' rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' Strict(trees)
#'
#' @seealso Less conservative summaries: [`Majority()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom TreeTools Consensus
#' @export
Strict <- function(trees) {
  Consensus(trees, p = 1)
}

#' Majority-rule consensus tree
#'
#' `Majority()` returns the majority-rule consensus
#' \insertCite{MargushMcMorris1981}{ConsTree}: the tree containing each split
#' that occurs in more than half of the input trees.  Raising `p` retains only
#' splits present in a greater proportion of trees, up to the strict consensus
#' at `p = 1`.
#'
#' A thin wrapper around [`TreeTools::Consensus()`].
#'
#' @inheritParams Strict
#' @param p Numeric between 0.5 and 1: the minimum proportion of trees that must
#' contain a split for it to be retained.  `p = 0.5` (the default) gives the
#' majority-rule consensus; `p = 1` gives the strict consensus.
#'
#' @return `Majority()` returns the consensus tree, an object of class `phylo`,
#' rooted as in the first entry of `trees`.
#'
#' @examples
#' trees <- ape::as.phylo(0:5, 8)
#' Majority(trees)
#'
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom TreeTools Consensus
#' @export
Majority <- function(trees, p = 0.5) {
  Consensus(trees, p = p)
}

#' @rdname Majority
#' @export
MajorityRule <- Majority

