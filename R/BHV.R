#' @useDynLib Consensus, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Validate that all trees share one leaf set and return a canonical ordering.
.BHVTipLabels <- function(trees) {
  labels <- lapply(trees, function(t) t[["tip.label"]])
  ref <- labels[[1]]
  for (l in labels) {
    if (length(l) != length(ref) || !setequal(l, ref)) {
      stop("All trees must share the same leaf labels.")
    }
  }
  sort(ref)
}

# Extract a tree into the BHV representation consumed by the C++ core:
# a 0/1 split-membership matrix (rows = interior splits, columns = `tipLabels`
# order), the matching interior-edge lengths, and a per-tip pendant length.
# BHV treespace is built on bipartitions (Brown & Owen treat trees as rooted by
# fixing one leaf), so degree-2 nodes are suppressed first: an internal singleton
# or a degree-2 root has two edges inducing the same split, which must be combined
# into a single coordinate (otherwise the tree lands at the wrong point).
#' @importFrom ape collapse.singles is.rooted unroot
#' @importFrom TreeTools as.Splits
.TreeToBHV <- function(tree, tipLabels) {
  if (is.null(tree[["edge.length"]])) {
    stop("Trees must have `edge.length` to compute BHV distances.")
  }
  tree <- collapse.singles(tree)
  if (is.rooted(tree)) {
    tree <- unroot(tree)
  }
  sp <- as.Splits(tree, tipLabels = tipLabels)
  nTip <- length(tipLabels)
  if (length(sp) == 0L) {
    membership <- matrix(0L, 0L, nTip)
    lengths <- numeric(0)
  } else {
    M <- as.logical(sp)
    if (is.null(dim(M))) {
      M <- matrix(as.integer(M), nrow = 1L)
    } else {
      M <- matrix(as.integer(M), nrow = nrow(M))
    }
    membership <- M
    nodeOfSplit <- as.integer(rownames(sp))
    lengths <- tree[["edge.length"]][match(nodeOfSplit, tree[["edge"]][, 2])]
  }
  leaf <- numeric(nTip)
  for (i in seq_len(nTip)) {
    tipNum <- match(tipLabels[i], tree[["tip.label"]])
    edgeRow <- match(tipNum, tree[["edge"]][, 2])
    leaf[i] <- if (is.na(edgeRow)) 0 else tree[["edge.length"]][edgeRow]
  }
  list(membership = membership, lengths = lengths, leaf = leaf)
}

# Rebuild a `phylo` from a BHV representation (C++ `to_r` output): clade
# membership (canonical, excluding the reference tip), interior lengths and
# pendant lengths.  Edges shorter than the tolerance have already been dropped.
#' @importFrom ape as.phylo
#' @importFrom TreeTools as.Splits StarTree
.BHVToTree <- function(rep, tipLabels) {
  nTip <- length(tipLabels)
  membership <- rep[["membership"]]
  if (is.null(nrow(membership)) || nrow(membership) == 0L) {
    tree <- StarTree(tipLabels)
    tree[["edge.length"]] <- rep[["leaf"]][match(tree[["tip.label"]], tipLabels)]
    return(tree)
  }
  spl <- as.Splits(membership > 0, tipLabels = tipLabels)
  tree <- as.phylo(spl)
  edge <- tree[["edge"]]
  nEdge <- nrow(edge)
  edgeLen <- numeric(nEdge)
  treeSplits <- as.Splits(tree, tipLabels = tipLabels)
  splitNode <- as.integer(rownames(treeSplits))
  # match rebuilt splits to the input clades to recover interior lengths
  hit <- match(treeSplits, spl)
  for (s in seq_along(splitNode)) {
    if (!is.na(hit[s])) {
      row <- which(edge[, 2] == splitNode[s])
      if (length(row)) edgeLen[row] <- rep[["lengths"]][hit[s]]
    }
  }
  for (tip in seq_len(nTip)) {
    label <- tree[["tip.label"]][tip]
    row <- which(edge[, 2] == tip)
    if (length(row)) edgeLen[row] <- rep[["leaf"]][match(label, tipLabels)]
  }
  tree[["edge.length"]] <- edgeLen
  # Return:
  tree
}

#' Geodesic (BHV) distance between two trees
#'
#' `BHVDistance()` returns the geodesic distance between two phylogenetic trees
#' with edge lengths in the Billera-Holmes-Vogtmann (BHV) treespace
#' \insertCite{BilleraHolmesVogtmann2001}{Consensus}, computed with the
#' polynomial-time GTP algorithm of \insertCite{OwenProvan2011;textual}{Consensus}.
#'
#' @param tree1,tree2 Trees of class `phylo`, sharing the same leaf labels and
#' carrying `edge.length`.
#'
#' @return `BHVDistance()` returns a single non-negative number: the length of
#' the geodesic between `tree1` and `tree2` in BHV treespace.
#'
#' @examples
#' set.seed(0)
#' edged <- function(tree) {
#'   tree$edge.length <- runif(nrow(tree$edge))
#'   tree
#' }
#' t1 <- edged(TreeTools::RandomTree(8, root = FALSE))
#' t2 <- edged(TreeTools::RandomTree(t1$tip.label, root = FALSE))
#' BHVDistance(t1, t2)
#'
#' @references \insertAllCited{}
#' @family BHV summaries
#' @export
BHVDistance <- function(tree1, tree2) {
  tl <- .BHVTipLabels(list(tree1, tree2))
  a <- .TreeToBHV(tree1, tl)
  b <- .TreeToBHV(tree2, tl)
  cpp_bhv_distance(a[["membership"]], a[["lengths"]], a[["leaf"]],
                   b[["membership"]], b[["lengths"]], b[["leaf"]])
}

# Tree at parameter `lambda` along the geodesic (chiefly for testing/teaching).
.BHVTreeAt <- function(tree1, tree2, lambda) {
  tl <- .BHVTipLabels(list(tree1, tree2))
  a <- .TreeToBHV(tree1, tl)
  b <- .TreeToBHV(tree2, tl)
  rep <- cpp_bhv_tree_at(a[["membership"]], a[["lengths"]], a[["leaf"]],
                         b[["membership"]], b[["lengths"]], b[["leaf"]], lambda)
  .BHVToTree(rep, tl)
}

# Coerce `trees` (a `multiPhylo`, a list, or a single `phylo`) to a plain list.
.BHVTreeList <- function(trees) {
  if (inherits(trees, "phylo")) {
    return(list(trees))
  }
  trees <- unclass(trees)
  if (!all(vapply(trees, inherits, logical(1), "phylo"))) {
    stop("`trees` must be a `phylo`, a `multiPhylo`, or a list of `phylo`.")
  }
  # Return:
  trees
}

#' Distances between every pair of trees
#'
#' `BHVPairwiseDistances()` returns the geodesic ([BHVDistance()]) between each
#' pair of trees in a set, as a `dist` object.
#'
#' @param trees A list of trees, or a `multiPhylo` object; all entries must
#' share the same leaf labels and carry `edge.length`.
#'
#' @return `BHVPairwiseDistances()` returns a [`stats::dist`] object of geodesic
#' distances.
#'
#' @examples
#' set.seed(0)
#' trees <- lapply(1:4, function(i) {
#'   tree <- TreeTools::RandomTree(7, root = FALSE)
#'   tree$edge.length <- runif(nrow(tree$edge))
#'   tree
#' })
#' BHVPairwiseDistances(trees)
#'
#' @references \insertAllCited{}
#' @family BHV summaries
#' @importFrom stats as.dist
#' @export
BHVPairwiseDistances <- function(trees) {
  trees <- .BHVTreeList(trees)
  tl <- .BHVTipLabels(trees)
  reps <- lapply(trees, .TreeToBHV, tipLabels = tl)
  n <- length(reps)
  d <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    a <- reps[[i]]
    for (j in (i + 1L):n) {
      b <- reps[[j]]
      d[i, j] <- d[j, i] <- cpp_bhv_distance(
        a[["membership"]], a[["lengths"]], a[["leaf"]],
        b[["membership"]], b[["lengths"]], b[["leaf"]])
    }
  }
  as.dist(d)
}

#' Fréchet mean and variance in BHV treespace
#'
#' `BHVMean()` returns the Fréchet (Karcher) mean of a set of trees in BHV
#' treespace: the tree that minimizes the sum of squared geodesic distances to
#' the sample \insertCite{BrownOwen2020}{Consensus}.  As there is no known
#' closed form, it is approximated by the iterative law-of-large-numbers
#' algorithm of \insertCite{Sturm2003;textual}{Consensus} and
#' \insertCite{MillerOwenProvan2015;textual}{Consensus}: starting from a sample
#' tree, each step walks a fraction \eqn{1/(k+1)} of the way along the geodesic
#' towards a randomly chosen sample tree.
#'
#' `BHVVariance()` returns the Fréchet variance: by default the mean squared
#' geodesic distance from the sample to its mean, \eqn{(1/r)\sum_i d(\bar
#' T, T_i)^2}; with `type = "sum"`, the total \eqn{\sum_i d(\bar T, T_i)^2}.
#'
#' The mean is "sticky": perturbing one sample tree need not move it, and it is
#' pulled towards lower-dimensional (less resolved) orthants, so it may be
#' unresolved even when the sample trees are binary
#' \insertCite{BrownOwen2020}{Consensus}.
#'
#' @inheritParams BHVPairwiseDistances
#' @param tolerance Numeric convergence threshold, _relative_ to the sample
#' standard deviation: iteration stops once `cauchyLength` consecutive steps
#' each move the estimate less than `tolerance` times the sample standard
#' deviation.  Smaller values give a more precise mean at the cost of more
#' iterations.
#' @param maxIter Integer specifying the maximum number of iterations.
#' @param cauchyLength Integer specifying the number of consecutive small steps
#' required to declare convergence.
#'
#' @return `BHVMean()` returns the mean tree, an object of class `phylo`, with
#' attributes `iterations` (number of steps taken) and `converged`.  Because the
#' step length shrinks as \eqn{1/(k+1)}, `converged = TRUE` indicates that
#' successive estimates have stopped moving appreciably (the stopping rule was
#' met before `maxIter`), not a guaranteed bound on the distance to the exact
#' mean; tighten `tolerance` for greater precision.
#'
#' @examples
#' set.seed(0)
#' trees <- lapply(1:25, function(i) {
#'   tree <- TreeTools::RandomTree(6, root = FALSE)
#'   tree$edge.length <- runif(nrow(tree$edge))
#'   tree
#' })
#' mean <- BHVMean(trees)
#' BHVVariance(trees, mean = mean)
#'
#' @references \insertAllCited{}
#' @family BHV summaries
#' @export
BHVMean <- function(trees, tolerance = 1e-4, maxIter = 100000L,
                    cauchyLength = 10L) {
  trees <- .BHVTreeList(trees)
  tl <- .BHVTipLabels(trees)
  reps <- lapply(trees, .TreeToBHV, tipLabels = tl)
  res <- cpp_bhv_mean(
    lapply(reps, `[[`, "membership"),
    lapply(reps, `[[`, "lengths"),
    lapply(reps, `[[`, "leaf"),
    length(tl), as.integer(maxIter), tolerance, as.integer(cauchyLength))
  tree <- .BHVToTree(res, tl)
  attr(tree, "iterations") <- res[["iterations"]]
  attr(tree, "converged") <- res[["converged"]]
  # Return:
  tree
}

#' @rdname BHVMean
#' @param mean Object of class `phylo` specifying a pre-computed mean tree;
#' computed via [BHVMean()] if `NULL` (the default).
#' @param type Character specifying whether to return the mean squared distance
#' (`"average"`, the default) or the total squared distance (`"sum"`).
#' @return `BHVVariance()` returns a single non-negative number.
#' @export
BHVVariance <- function(trees, mean = NULL, type = c("average", "sum")) {
  type <- match.arg(type)
  trees <- .BHVTreeList(trees)
  if (is.null(mean)) {
    mean <- BHVMean(trees)
  }
  ss <- sum(vapply(trees, function(t) BHVDistance(mean, t)^2, numeric(1)))
  # Return:
  if (type == "average") ss / length(trees) else ss
}
