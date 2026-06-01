#' Average consensus tree
#'
#' `Average()` returns the *average consensus*
#' \insertCite{LapointeCucumel1997}{Consensus}: the tree whose path-length
#' (patristic) distances most closely match the average of the path-length
#' distances of the input trees.  Informally, it places each leaf at its mean
#' position across the input trees, making it a natural distance-based summary of
#' a posterior sample -- complementing the split-based [`Strict()`] and
#' [`Majority()`] methods.
#'
#' The procedure has two steps \insertCite{LapointeCucumel1997}{Consensus}:
#'
#' 1. Compute the path-length distance matrix of each input tree (using branch
#'    lengths where present, otherwise counting edges), optionally rescaling each
#'    matrix, and average the matrices.
#' 2. Find the tree whose own path-length distances best fit this average matrix,
#'    in the least-squares sense.
#'
#' Because the average of several path-length matrices is usually not itself
#' realisable by any tree (it violates the four-point condition), step 2 is a
#' fit, not an inversion.  By default `Average()` approximates it with the fast
#' balanced minimum-evolution tree; `method = "ls"` instead performs the exact
#' least-squares search, which -- being NP-hard \insertCite{Day1987}{Consensus}
#' -- uses tree rearrangements, as did the original \acronym{FITCH}
#' implementation.
#'
#' @inheritParams Strict
#' @param method Character specifying how to build the tree from the average
#' distance matrix:
#' - `"fastme.bal"` (the default) returns the balanced minimum-evolution tree
#'   \insertCite{DesperGascuel2002}{Consensus}: a fast, accurate approximation of
#'   the least-squares tree;
#' - `"ls"` searches for the least-squares tree itself -- the criterion under
#'   which Lapointe & Cucumel's averaging guarantee holds -- using
#'   \pkg{TreeSearch} and \pkg{phangorn}; slower, and currently a single-start
#'   rearrangement search;
#' - `"nj"`, `"bionj"` and `"fastme.ols"` return the corresponding distance tree
#'   \insertCite{SaitouNei1987,Gascuel1997}{Consensus}.
#' @param weights Optional numeric vector, one entry per tree, giving the weight
#' of each tree in the average (e.g. posterior probabilities).  The default,
#' `NULL`, weights every tree equally &ndash; appropriate for a posterior sample,
#' in which a tree's frequency already encodes its probability.
#' @param scale Character specifying whether to rescale each tree's distance
#' matrix before averaging.  `"none"` (the default) leaves matrices unscaled,
#' appropriate when the trees are already commensurable (e.g. a single posterior
#' sample).  `"max"` divides each matrix by its largest entry, the
#' standardization recommended by Lapointe & Cucumel when combining trees from
#' heterogeneous sources whose absolute distances are not comparable.
#' @param edgeLengths Logical specifying whether to use branch lengths when
#' computing path-length distances.  The default, `NA`, uses branch lengths when
#' *every* tree has them and otherwise counts edges; `TRUE` requires branch
#' lengths; `FALSE` always counts edges (a topology-only summary).
#' @param outgroup Optional tip label(s) on which to root the result.  The
#' default, `NULL`, returns an unrooted tree: path-length distances are
#' unaffected by rooting, so the method is intrinsically unrooted.
#' @param check.labels Logical specifying whether to confirm that every tree
#' describes the same leaves.  The default, `TRUE`, is safer; `FALSE` is faster
#' when the trees are known to share an identical leaf set.
#' @param \dots Further arguments (e.g. `maxIter`, `maxHits`, `EdgeSwapper`)
#' passed to [`TreeSearch::TreeSearch()`] when `method = "ls"`.
#'
#' @return `Average()` returns the average consensus tree, an object of class
#' `phylo` with fitted branch lengths, unrooted unless `outgroup` is given.
#'
#' @examples
#' trees <- ape::rmtree(5, 8)         # five random eight-leaf trees
#' Average(trees)                     # fast (balanced minimum evolution) default
#' \donttest{
#' Average(trees, method = "ls")      # faithful least-squares fit (slower)
#' }
#'
#' @seealso Split-based summaries: [`Strict()`], [`Majority()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape bionj cophenetic.phylo fastme.bal fastme.ols nj root unroot
#' @importFrom stats as.dist setNames
#' @importFrom TreeTools RenumberTips TipLabels
#' @export
Average <- function(trees,
                    method = c("fastme.bal", "ls", "nj", "bionj", "fastme.ols"),
                    weights = NULL,
                    scale = c("none", "max"),
                    edgeLengths = NA,
                    outgroup = NULL,
                    check.labels = TRUE,
                    ...) {
  method <- match.arg(method)
  scale <- match.arg(scale)

  trees <- .AverageTrees(trees, check.labels)
  if (length(trees) == 1L) {
    return(.RootResult(trees[[1]], outgroup))
  }

  labs <- TipLabels(trees[[1]])
  averageDist <- .AverageDistance(trees, labs, weights, scale, edgeLengths)
  tree <- .FitTree(averageDist, method, labs, ...)

  .RootResult(tree, outgroup)
}

# Coerce input to a `multiPhylo` of trees on a shared, consistently ordered leaf
# set.  Unlike `TreeTools::Consensus()`, branch lengths are *retained*, as the
# average consensus depends on them.
.AverageTrees <- function(trees, check.labels) {
  if (inherits(trees, "phylo")) {
    return(structure(list(trees), class = "multiPhylo"))
  }
  if (!is.list(trees) || is.data.frame(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  trees <- structure(unclass(trees), class = "multiPhylo")
  if (length(trees) == 0L) {
    stop("`trees` contains no trees.")
  }

  labs1 <- TipLabels(trees[[1]])
  if (length(trees) > 1L && isTRUE(check.labels)) {
    for (i in seq_along(trees)[-1]) {
      if (!setequal(TipLabels(trees[[i]]), labs1)) {
        stop("All trees must share the same leaf labels; `Average()` does not ",
             "yet support overlapping or partial taxon sets.")
      }
    }
    trees <- RenumberTips(trees, labs1)
  }

  trees
}

# Path-length distance matrix of each tree, in a fixed label order, optionally
# rescaled; then a (weighted) mean.
.AverageDistance <- function(trees, labs, weights, scale, edgeLengths) {
  haveLengths <- vapply(trees, function(tr) !is.null(tr[["edge.length"]]),
                        logical(1))
  useLengths <- if (is.na(edgeLengths)) all(haveLengths) else isTRUE(edgeLengths)
  if (isTRUE(edgeLengths) && !all(haveLengths)) {
    stop("`edgeLengths = TRUE`, but not every tree has branch lengths.")
  }

  matrices <- lapply(trees, function(tr) {
    if (!useLengths) {
      tr[["edge.length"]] <- rep_len(1, nrow(tr[["edge"]]))
    }
    m <- cophenetic.phylo(tr)[labs, labs]
    if (scale == "max") {
      maxVal <- max(m)
      if (maxVal > 0) {
        m <- m / maxVal
      }
    }
    m
  })

  if (is.null(weights)) {
    weights <- rep_len(1, length(matrices))
  } else {
    if (length(weights) != length(matrices)) {
      stop("`weights` must have one entry per tree.")
    }
    if (anyNA(weights) || any(weights < 0) || sum(weights) <= 0) {
      stop("`weights` must be non-negative and not all zero.")
    }
  }
  weights <- weights / sum(weights)

  Reduce(`+`, Map(function(m, w) m * w, matrices, weights))
}

# Build a tree from an average distance matrix.
.FitTree <- function(averageDist, method, labs, ...) {
  d <- as.dist(averageDist)
  switch(method,
    nj = nj(d),
    bionj = bionj(d),
    fastme.bal = fastme.bal(d),
    fastme.ols = fastme.ols(d),
    ls = .LeastSquaresTree(averageDist, labs, ...)
  )
}

# Heuristic least-squares tree: rearrangement search scored by the residual sum
# of squares of NNLS-fitted branch lengths.  This R-level search is correct but
# slow; it is the interim engine pending a compiled least-squares scorer in
# 'TreeSearch'.
.LeastSquaresTree <- function(averageDist, labs, start = NULL,
                              EdgeSwapper = NULL, maxIter = 100L, maxHits = 20L,
                              verbosity = 0L, ...) {
  if (!requireNamespace("TreeSearch", quietly = TRUE)) {
    stop("`method = \"ls\"` requires the 'TreeSearch' package.")
  }
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    stop("`method = \"ls\"` requires the 'phangorn' package.")
  }
  d <- as.dist(averageDist)
  n <- length(labs)

  # Fewer than four leaves: a single (unrooted) topology, so just fit lengths.
  if (n < 4L) {
    return(.NnlsRefit(nj(d), averageDist, labs))
  }

  if (is.null(start)) {
    start <- fastme.bal(d)
  }
  start <- root(start, outgroup = labs[[1]], resolve.root = TRUE)
  start[["edge.length"]] <- NULL
  if (is.null(EdgeSwapper)) {
    EdgeSwapper <- TreeSearch::NNISwap
  }

  dataset <- setNames(vector("list", n), labs)
  attr(dataset, "distances") <- as.matrix(averageDist)[labs, labs]

  best <- TreeSearch::TreeSearch(
    start, dataset,
    InitializeData = function(dataset) dataset,
    CleanUpData = function(x) NULL,
    TreeScorer = .LeastSquaresScore,
    EdgeSwapper = EdgeSwapper,
    maxIter = maxIter, maxHits = maxHits, verbosity = verbosity, ...
  )

  .NnlsRefit(best, averageDist, labs)
}

# Residual sum of squares of the NNLS fit of a candidate topology (supplied as
# `parent`/`child` edge vectors) to the target distance matrix.
.LeastSquaresScore <- function(parent, child, dataset, ...) {
  distances <- attr(dataset, "distances")
  labs <- names(dataset)
  tree <- unroot(structure(
    list(edge = cbind(parent, child),
         Nnode = length(unique(parent)),
         tip.label = labs),
    class = "phylo"))
  fit <- phangorn::nnls.tree(as.dist(distances), tree, method = "unrooted",
                             trace = 0L)
  fitted <- cophenetic.phylo(fit)[labs, labs]
  sum((fitted[lower.tri(fitted)] - distances[lower.tri(distances)]) ^ 2)
}

# Fit non-negative least-squares branch lengths to a fixed (unrooted) topology.
.NnlsRefit <- function(tree, averageDist, labs) {
  phangorn::nnls.tree(as.dist(as.matrix(averageDist)[labs, labs]),
                      unroot(tree), method = "unrooted", trace = 0L)
}

# Return unrooted (the default) or rooted on the requested outgroup.
.RootResult <- function(tree, outgroup) {
  if (is.null(outgroup)) {
    return(tree)
  }
  root(tree, outgroup = outgroup, resolve.root = TRUE)
}
