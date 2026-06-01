#' Information-maximizing quartet consensus tree
#'
#' Construct a consensus tree that maximizes the net concordant quartet
#' information shared with a set of input trees, using a greedy add-and-prune
#' heuristic that can also drop rogue taxa.
#'
#' Each quartet (four-leaf statement) carries one trit of information.  For
#' every quartet that the consensus resolves, `QuartetConsensus()` scores
#' \eqn{\mathrm{agree} - \mathrm{penalty} \times \mathrm{disagree}}, where,
#' among the input trees that resolve the quartet, _agree_ is the number that
#' resolve it as the consensus does and _disagree_ the number that resolve it
#' differently (input polytomies neither agree nor disagree).  Quartets that
#' the consensus leaves unresolved, or that involve a dropped taxon, score
#' zero.  The objective is the **absolute sum** of these signed contributions
#' (no normalization by quartet count), so the star tree scores zero and the
#' criterion is calibrated for both resolution and leaf count: a clade or a
#' taxon is retained only when it adds net positive information, so the score
#' neither rises by adding uninformative taxa nor collapses to a trivially
#' resolved few-leaf tree \insertCite{Smith2019,Smith2022}{Consensus}.
#'
#' With the default `penalty = 1` and binary input trees this recovers the
#' approximate median tree under the symmetric quartet distance
#' \insertCite{Takazawa2026}{Consensus}.
#'
#' Because the quartet objective gives greater weight to deep branches (which
#' resolve more quartets), quartet consensus trees tend to be more resolved
#' than majority-rule trees, especially when phylogenetic signal is low.
#'
#' @param trees An object of class `multiPhylo`: the input trees.
#'   All trees must share the same tip labels.
#'   Trees may be non-binary (polytomies are handled correctly).
#' @param init Character string specifying the initial tree:
#'   - `"majority"` (default): start from the majority-rule consensus.
#'   - `"empty"`: start from a star tree (purely additive).
#'   - `"extended"`: start from the extended (greedy) majority-rule consensus.
#' @param greedy Character string specifying the greedy strategy:
#'   - `"best"` (default): evaluate all candidates and pick the single
#'     highest-benefit action at each step.
#'   - `"first"`: pick the first improving action encountered (faster but
#'     may give a slightly worse result).
#' @param neverDrop Controls rogue taxon dropping:
#'   - `TRUE` (default): never drop taxa.
#'   - `FALSE`: any taxon may be dropped.
#'   - A character vector of tip labels: those tips are protected from
#'     dropping; all other tips are candidates.
#' @param penalty Numeric: the misinformation penalty, equal to the ratio
#'   \eqn{b / a} of the cost of a misleading quartet to the reward for a
#'   correct one.  A quartet is resolved only when the proportion of resolving
#'   input trees that support one resolution exceeds
#'   \eqn{\mathrm{penalty} / (1 + \mathrm{penalty})}.  The default `1` requires
#'   a simple majority (> 1/2, "penalize misinformation as much as you reward
#'   information"); `penalty = 0.5` requires only that a resolution beat the
#'   1-in-3 random baseline (> 1/3, the chance-corrected threshold); larger
#'   values are progressively more conservative.
#'
#' @details
#' The algorithm pools all splits observed across input trees and maintains
#' a quartet profile: for each of the \eqn{\binom{n}{4}}{C(n,4)} quartets,
#' a count of how many input trees resolve it as each of the three possible
#' topologies.  Splits are greedily added to (or removed from) the consensus
#' when doing so increases the objective; candidate splits must be compatible
#' with all currently included splits.  When `neverDrop` enables taxon
#' dropping, the greedy loop also considers removing rogue taxa, taking the
#' single best-improving action (add split, remove split, or drop taxon) at
#' each step.  A rogue is dropped only when the quartets it forces the
#' consensus to resolve are, on balance, misleading; it can equally be left in
#' place within a polytomy when that conveys more information.
#'
#' The function supports trees with up to 100 tips.  For larger trees,
#' the explicit quartet enumeration becomes prohibitively expensive.
#'
#' @return A tree of class `phylo`.  When taxon dropping is enabled,
#'   the tree may have fewer tips than the input trees, and attributes
#'   `"dropped"` (character vector of dropped tip labels, in drop order)
#'   and `"drop_scores"` (the objective score after each drop; higher is
#'   better) are attached.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library("TreeTools")
#'
#' # Generate bootstrap-like trees
#' trees <- as.phylo(1:30, nTip = 8)
#'
#' # Quartet consensus
#' qc <- QuartetConsensus(trees)
#' plot(qc)
#'
#' # Compare resolution with majority-rule
#' mr <- UnrootTree(Consensus(trees, p = 0.5))
#' cat("Majority-rule splits:", NSplits(mr), "\n")
#' cat("Quartet consensus splits:", NSplits(qc), "\n")
#'
#' @importFrom ape as.phylo
#' @importFrom TreeTools as.Splits TipLabels NSplits Consensus StarTree
#' @export
QuartetConsensus <- function(trees,
                             init = c("majority", "empty", "extended"),
                             greedy = c("best", "first"),
                             neverDrop = TRUE,
                             penalty = 1) {
  init <- match.arg(init)
  greedy <- match.arg(greedy)

  if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be an object of class 'multiPhylo'.")
  }
  nTree <- length(trees)
  if (nTree < 2L) stop("Need at least 2 trees.")

  if (!is.numeric(penalty) || length(penalty) != 1L ||
      !is.finite(penalty) || penalty < 0) {
    stop("`penalty` must be a single non-negative number.")
  }

  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)

  if (nTip < 4L) stop("Need at least 4 tips for quartet consensus.")
  if (nTip > 100L) {
    stop("QuartetConsensus supports at most 100 tips. ",
         "The explicit quartet enumeration is O(n^4).")
  }

  for (i in seq_len(nTree)[-1]) {
    labs_i <- TipLabels(trees[[i]])
    if (!setequal(labs_i, tipLabels)) {
      extra   <- setdiff(labs_i, tipLabels)
      missing <- setdiff(tipLabels, labs_i)
      stop("Tree ", i, " has different tip labels from tree 1.",
           if (length(missing)) paste0("\n  Missing in tree ", i, ": ",
                                       paste(missing, collapse = ", ")),
           if (length(extra))   paste0("\n  Unexpected in tree ", i, ": ",
                                       paste(extra, collapse = ", ")))
    }
  }

  # Resolve neverDrop to an integer vector (1-based) or NULL
  if (isTRUE(neverDrop)) {
    neverDropR <- NULL
  } else if (isFALSE(neverDrop)) {
    neverDropR <- integer(0)
  } else {
    # Character vector of protected labels
    neverDrop <- as.character(neverDrop)
    bad <- setdiff(neverDrop, tipLabels)
    if (length(bad)) {
      stop("neverDrop labels not found in trees: ",
           paste(bad, collapse = ", "))
    }
    neverDropR <- match(neverDrop, tipLabels)
  }

  # Convert each tree to a raw split matrix
  splitsList <- lapply(trees, function(tr) {
    sp <- as.Splits(tr, tipLabels)
    unclass(sp)
  })

  res <- cpp_quartet_consensus(
    splitsList, nTip,
    init_majority = (init == "majority"),
    init_extended = (init == "extended"),
    greedy_best_flag = (greedy == "best"),
    never_drop_r = neverDropR,
    penalty_r = penalty
  )

  # Build tree from pre-filtered splits (already remapped to active tips)
  activeTipLabels <- tipLabels[res$active_tips]
  nActiveTip <- res$n_active
  splitsMat <- res$splits

  if (nrow(splitsMat) == 0L) {
    result <- StarTree(activeTipLabels)
  } else {
    sp <- structure(splitsMat, nTip = nActiveTip,
                    tip.label = activeTipLabels, class = "Splits")
    result <- as.phylo(sp)
  }

  # Attach drop metadata
  if (!is.null(neverDropR)) {
    droppedTipIdx <- res$dropped_tips
    attr(result, "dropped") <- tipLabels[droppedTipIdx]
    attr(result, "drop_scores") <- res$drop_scores
  }

  result
}
