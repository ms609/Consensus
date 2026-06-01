#' R* consensus tree
#'
#' `RStar()` returns the R* consensus \insertCite{Degnan2009}{ConsTree} of a set
#' of **rooted** trees.
#'
#' The R* consensus is a rooted-triplet method.  For every set of three leaves it
#' tallies, across the input trees, the three possible resolved rooted triplets
#' (`ab|c`, `ac|b`, `bc|a`) and keeps the one that is *uniquely favoured*: the
#' resolution appearing in strictly more input trees than **each** of the other
#' two, considered separately.  This is a strict **plurality** rule, not a
#' majority rule: a triplet can be uniquely favoured with far fewer than half the
#' votes (e.g. a 2--1--1 split among four trees), and any tie leaves those taxa
#' unresolved.  The kept triplets form the set of *majority resolved triplets*,
#' \eqn{R_{maj}}.  The R* tree is then the unique tree whose clades are exactly
#' the **strong clusters** of \eqn{R_{maj}}
#' \insertCite{Degnan2009,Jansson2016a}{ConsTree}: a leaf set `A` is a clade if
#' and only if, for *every* pair of leaves in `A` and *every* leaf `x` outside
#' `A`, the triplet grouping that pair against `x` is uniquely favoured.
#' Equivalently, R* is the most resolved tree that displays no resolved triplet
#' outside \eqn{R_{maj}}.
#'
#' R* is always a *refinement* of the majority-rule consensus
#' (every majority clade also appears in `RStar()`) and is a statistically
#' consistent estimator of a species tree from gene trees.  Unlike [`Local()`]
#' it is **polynomial**, so it is not restricted to small leaf counts; `RStar()`
#' caps `n` at 200 purely as a memory safeguard on its dense triplet tensor (a
#' limit quite different in nature from `Local()`'s exact-exponential 20-leaf
#' bound).
#'
#' Like [`Adams()`], R* is a rooted method: triplet states depend on the
#' rooting, so input trees are treated as rooted on their current root.  Root the
#' trees as you intend before calling `RStar()`.
#'
#' @section Construction and conventions:
#' \itemize{
#'   \item *Fans.*  When a non-binary input tree leaves three leaves unresolved
#'     (a fan), that tree does not count toward any resolution of that triplet;
#'     fans have no impact on \eqn{R_{maj}}
#'     \insertCite{Jansson2016a}{ConsTree}.
#'   \item *Ties.*  If no resolution of a triple is uniquely favoured, that
#'     triple contributes nothing, leaving the affected taxa unresolved.
#'   \item *Existence and uniqueness.*  The strong-cluster construction always
#'     yields a single, well-defined tree (Lemma 1.1 of
#'     \insertCite{Jansson2016a}{ConsTree}); there is no incompatibility or
#'     "build-failure" case to resolve.
#'   \item *Algorithm.*  This implementation is correctness-first: an
#'     \eqn{O(kn^3)} triplet tally followed by an \eqn{O(n^4)} strong-cluster
#'     assembly.  The sub-cubic and near-quadratic algorithms of
#'     \insertCite{Jansson2013a,Jansson2016a}{ConsTree} are a deferred speed
#'     optimisation.
#' }
#'
#' @inheritParams Strict
#'
#' @return `RStar()` returns the consensus tree, an object of class `phylo`,
#' rooted by construction.
#'
#' @examples
#' # Five trees whose majority signal recovers the species tree (((a,b),c),d):
#' trees <- c(
#'   ape::read.tree(text = "(((a, b), c), d);"),
#'   ape::read.tree(text = "(((a, b), c), d);"),
#'   ape::read.tree(text = "(((a, b), c), d);"),
#'   ape::read.tree(text = "(((a, c), b), d);"),
#'   ape::read.tree(text = "(((b, c), a), d);")
#' )
#' RStar(trees) # (a, b) wins {a,b,c} by plurality (3 vs 1 vs 1)
#'
#' @seealso Closely related: [`Strict()`], [`Majority()`], [`Adams()`],
#' [`Local()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape read.tree
#' @importFrom TreeTools Preorder RenumberTips TipLabels
#' @export
RStar <- function(trees) {
  # Coerce to a plain list.
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
    return(if (nTree) trees[[1L]] else NULL)
  }

  labels <- TipLabels(trees[[1L]])
  n <- length(labels)

  if (n < 3L) {
    return(trees[[1L]])
  }

  if (n > 200L) {
    stop(
      "RStar() caps the dense triplet tensor at 200 leaves (n = ", n, ")."
    )
  }

  # Relabel every tree 1..n in a shared canonical order and put in Preorder so
  # the C++ LCA pass (which assumes parent precedes child) is valid.
  edgeList <- lapply(trees, function(tr) {
    tr <- Preorder(RenumberTips(tr, labels))
    tr[["edge"]]
  })

  nwk <- rStarConsensus(edgeList, n)

  tree <- read.tree(text = paste0(nwk, ";"))
  tree[["tip.label"]] <- labels[as.integer(tree[["tip.label"]])]
  # Return:
  tree
}

