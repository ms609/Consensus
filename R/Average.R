#' Average consensus tree
#'
#' `Average()` returns the *average consensus*
#' \insertCite{LapointeCucumel1997}{ConsTree}: the tree whose path-length
#' (patristic) distances most closely match the average of the path-length
#' distances of the input trees.  Informally, it places each leaf at its mean
#' position across the input trees, making it a natural distance-based summary of
#' a posterior sample -- complementing the split-based [`Strict()`] and
#' [`Majority()`] methods.
#'
#' The procedure has two steps \insertCite{LapointeCucumel1997}{ConsTree}:
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
#' least-squares search, which -- being NP-hard \insertCite{Day1987}{ConsTree}
#' -- uses tree rearrangements, as did the original \acronym{FITCH}
#' implementation.  Branch lengths are fitted by non-negative least squares, so
#' that the fitted distances are realisable by a tree, as the criterion requires.
#'
#' A lone input tree is its own average: it is returned (unrooted unless
#' `outgroup` is given) without refitting, and `method`, `scale`, `weights` and
#' `edgeLengths` then have no effect.
#'
#' @inheritParams Strict
#' @param method Character specifying how to build the tree from the average
#' distance matrix:
#' - `"fastme.bal"` (the default) returns the balanced minimum-evolution tree
#'   \insertCite{DesperGascuel2002}{ConsTree}: a fast, accurate approximation of
#'   the least-squares tree;
#' - `"ls"` searches for the least-squares tree itself -- the criterion under
#'   which Lapointe & Cucumel's averaging guarantee holds -- using
#'   `TreeSearch::LeastSquaresTree()`, a compiled non-negative least-squares
#'   \acronym{NNI}/\acronym{SPR} search;
#' - `"nj"`, `"bionj"` and `"fastme.ols"` return the corresponding distance tree
#'   \insertCite{SaitouNei1987,Gascuel1997}{ConsTree}.
#' @param weights Numeric vector specifying the weight of each tree in the
#' average (e.g. posterior probabilities), with one entry per tree.  Defaults
#' to `NULL`, which weights every tree equally -- appropriate for a posterior
#' sample, in which a tree's frequency already encodes its probability.
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
#' @param outgroup Character vector specifying tip label(s) on which to root
#' the result.  Defaults to `NULL`, which returns an unrooted tree: path-length
#' distances are unaffected by rooting, so the method is intrinsically unrooted.
#' @param check.labels Logical specifying whether to confirm that every tree
#' describes the same leaves.  The default, `TRUE`, is safer; `FALSE` is faster
#' when the trees are known to share an identical leaf set.
#' @param lsControl Named list of further arguments for the least-squares
#' search (`method = "ls"`), passed to `TreeSearch::LeastSquaresTree()`; for
#' example `list(spr = FALSE, maxHits = 5L, weight = "fm")` to use
#' Fitch-Margoliash weighting.  Defaults to `list()`; ignored by other methods.
#'
#' @return `Average()` returns the average consensus tree, an object of class
#' `phylo` with fitted branch lengths, unrooted unless `outgroup` is given.
#'
#' @examples
#' trees <- ape::rmtree(5, 8)         # five random eight-leaf trees
#' Average(trees)                     # fast (balanced minimum evolution) default
#' \donttest{
#' if (requireNamespace("TreeSearch", quietly = TRUE) &&
#'     exists("LeastSquaresTree", where = asNamespace("TreeSearch"),
#'            mode = "function")) {
#'   Average(trees, method = "ls")    # faithful least-squares fit (slower)
#' }
#' }
#'
#' @seealso Split-based summaries: [`Strict()`], [`Majority()`].
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape bionj cophenetic.phylo fastme.bal fastme.ols is.rooted nj root unroot
#' @importFrom stats as.dist
#' @importFrom TreeTools RenumberTips TipLabels
#' @importFrom utils modifyList
#' @export
Average <- function(trees,
                    method = c("fastme.bal", "ls", "nj", "bionj", "fastme.ols"),
                    weights = NULL,
                    scale = c("none", "max"),
                    edgeLengths = NA,
                    outgroup = NULL,
                    check.labels = TRUE,
                    lsControl = list()) {
  method <- match.arg(method)
  scale <- match.arg(scale)

  trees <- .AverageTrees(trees, check.labels)
  if (length(trees) == 1L) {
    return(.RootResult(trees[[1]], outgroup))
  }

  labs <- TipLabels(trees[[1]])
  averageDist <- .AverageDistance(trees, labs, weights, scale, edgeLengths)
  tree <- .FitTree(averageDist, method, lsControl)

  # Return:
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

  # Return:
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
.FitTree <- function(averageDist, method, lsControl) {
  d <- as.dist(averageDist)
  switch(method,
    nj = nj(d),
    bionj = bionj(d),
    fastme.bal = fastme.bal(d),
    fastme.ols = fastme.ols(d),
    ls = .LeastSquaresTree(averageDist, lsControl)
  )
}

# Least squares: delegate the topology search to TreeSearch's compiled
# non-negative least-squares kernel.  `lsControl` is a named list of further
# arguments for `TreeSearch::LeastSquaresTree()` (e.g. `spr`, `maxHits`,
# `weight`); passing them through a list rather than `...` avoids R's partial
# matching silently binding them to `Average()`'s own formals.
.LeastSquaresTree <- function(averageDist, lsControl) {
  if (!requireNamespace("TreeSearch", quietly = TRUE) ||
      !exists("LeastSquaresTree", where = asNamespace("TreeSearch"),
              mode = "function")) {
    # nocov start
    # Defensive guard for installations with a 'TreeSearch' older than the
    # required (>= 2.0.0) build that provides `LeastSquaresTree()`; unreachable
    # whenever the declared dependency is satisfied (as in the coverage run).
    stop("`method = \"ls\"` requires a version of 'TreeSearch' that provides ",
         "`LeastSquaresTree()`.")
    # nocov end
  }
  if (!is.list(lsControl)) {
    stop("`lsControl` must be a named list.")
  }
  d <- as.dist(averageDist)
  if (length(lsControl) == 0L) {
    TreeSearch::LeastSquaresTree(d, method = "nnls")
  } else {
    do.call(TreeSearch::LeastSquaresTree,
            modifyList(list(d, method = "nnls"), lsControl))
  }
}

# Return unrooted (the default, honouring the rooting-invariant contract) or
# rooted on the requested outgroup.  Only re-root when necessary, to preserve
# attributes (e.g. the `"RSS"` set by the least-squares search).
.RootResult <- function(tree, outgroup) {
  if (is.null(outgroup)) {
    if (is.rooted(tree)) {
      tree <- unroot(tree)
    }
    return(tree)
  }
  # Return:
  root(tree, outgroup = outgroup, resolve.root = TRUE)
}

