#' Local consensus tree
#'
#' `Local()` returns the local consensus
#' \insertCite{JanssonRajabySung2018}{ConsTree} of a set of **rooted** trees.
#'
#' The local consensus is the most conservative tree consistent with the rooted
#' triplets shared by *every* input tree.  Two variants are offered: the minimum
#' rooted local consensus (MinRLC, `type = "rooted"`) and the minimum induced
#' local consensus (MinILC, `type = "induced"`), which differ in how the
#' resolution of the result is scored.  The tree is assembled from the Aho-graph
#' decomposition of the common-triplet set by dynamic programming over the
#' subsets of leaves.
#'
#' @details
#' The MinRLC and MinILC variants, and the exact exponential-time algorithm used
#' here to construct them, are due to \insertCite{JanssonRajabySung2018}{ConsTree};
#' this is a direct port of the reference C++ from their `FDCT_new` toolkit, used
#' with permission.
#'
#' **Complexity note:** the algorithm is exponential, so `Local()` is limited to
#' `n <= 20` leaves.  Running time depends not only on `n` but on how *congruent*
#' the input trees are: when few rooted triplets are shared (highly incongruent
#' trees, particularly with `type = "induced"`), the dynamic programming can
#' become intractable even within the 20-leaf limit.  A long-running call can be
#' interrupted (e.g. with Ctrl-C).
#'
#' **No-valid-consensus:** when the entire leaf set forms a single inseparable
#' Aho-graph component (no common triplets separate any pair), the algorithm has
#' no valid consensus; `Local()` then returns a star tree (all leaves attached
#' directly to the root, fully unresolved).  The reference binary reports
#' "No valid consensus found." for this case.
#'
#' Input trees are treated as rooted on their current root; root the trees as
#' you intend before calling `Local()`.
#'
#' @inheritParams Strict
#' @param type Character specifying whether to compute the minimum rooted local
#'   consensus (`"rooted"`, the default; MinRLC) or the minimum induced local
#'   consensus (`"induced"`; MinILC).
#'
#' @return `Local()` returns the consensus tree, an object of class `phylo`,
#'   rooted by construction.
#'
#' @examples
#' # Two trees that agree on one cherry but disagree on overall topology
#' t1 <- ape::read.tree(text = "(1,((2,3),4));")
#' t2 <- ape::read.tree(text = "(1,((2,4),3));")
#' Local(list(t1, t2), "rooted")  # keeps clade {2,3,4} only
#'
#' @seealso Closely related: [`Strict()`], [`Majority()`], [`Adams()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape read.tree
#' @importFrom TreeTools NTip Preorder RenumberTips TipLabels
#' @export
Local <- function(trees, type = c("rooted", "induced")) {
  type  <- match.arg(type)
  minrs <- (type == "rooted")

  # Coerce to a plain list
  if (inherits(trees, "phylo")) return(trees)
  if (!is.list(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  trees <- c(trees)
  trees <- trees[!vapply(trees, is.null, logical(1))]
  nTree <- length(trees)

  if (nTree < 2L) return(if (nTree) trees[[1L]] else NULL)

  labels <- TipLabels(trees[[1L]])
  if (any(vapply(trees[-1L], function(tr)
    !setequal(TipLabels(tr), labels), logical(1)))) {
    stop("all trees must have the same tip labels")
  }
  n      <- length(labels)

  if (n < 3L) return(trees[[1L]])

  if (n > 20L) {
    stop(
      "Local() is exact-exponential and limited to 20 leaves (n = ", n, ")."
    )
  }

  # Relabel every tree 1..n in a shared canonical order and put in Preorder.
  edgeList <- lapply(trees, function(tr) {
    tr <- Preorder(RenumberTips(tr, labels))
    tr[["edge"]]
  })

  nwk <- localConsensus(edgeList, n, minrs)

  if (!nzchar(nwk)) {
    # nocov start
    # No-valid-consensus sentinel (empty string from the C++ core): return a
    # star tree.  Unreachable for valid input -- the common-triplet set is a
    # subset of the first tree's triplets, hence consistent, so the Aho-graph
    # decomposition always separates and never reports a single inseparable
    # full component.  Retained to mirror the reference binary's behaviour.
    star <- ape::read.tree(
      text = paste0("(", paste(seq_len(n), collapse = ","), ");")
    )
    star[["tip.label"]] <- labels
    return(star)
    # nocov end
  }

  tree <- ape::read.tree(text = paste0(nwk, ";"))
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return:
  tree
}

