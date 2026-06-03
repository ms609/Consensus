## probe-edge.R — empirical red-team probes for Average()
## Run from repo root:
##   Rscript.exe dev/red-team/reviews/average-consensus/probe-edge.R

.libPaths(c(".agent-cons", .libPaths()))
library(ConsTree)
library(ape)

PASS <- function(label) cat(sprintf("PASS  %s\n", label))
FAIL <- function(label, detail = "") cat(sprintf("FAIL  %s  %s\n", label, detail))
INFO <- function(label, detail = "") cat(sprintf("INFO  %s  %s\n", label, detail))

## ----------------------------------------------------------------
## P1: n=3 tips — does every method complete without error?
## ----------------------------------------------------------------
t3a <- read.tree(text = "(A,(B,C));")
t3b <- read.tree(text = "((A,B),C);")
for (m in c("nj", "bionj", "fastme.bal", "fastme.ols")) {
  tryCatch({
    r <- Average(list(t3a, t3b), method = m)
    PASS(sprintf("n=3 method=%s: %d tips, class=%s", m, length(r$tip.label), class(r)))
  }, error = function(e) {
    FAIL(sprintf("n=3 method=%s", m), conditionMessage(e))
  })
}

## ----------------------------------------------------------------
## P2: n=2 tips — below triangle inequality floor for distance methods
## ----------------------------------------------------------------
t2a <- read.tree(text = "(A,B);")
t2b <- read.tree(text = "(A,B);")
for (m in c("nj", "bionj", "fastme.bal", "fastme.ols")) {
  tryCatch({
    r <- Average(list(t2a, t2b), method = m)
    PASS(sprintf("n=2 method=%s: %d tips", m, length(r$tip.label)))
  }, error = function(e) {
    FAIL(sprintf("n=2 method=%s", m), conditionMessage(e))
  })
}

## ----------------------------------------------------------------
## P3: star tree with zero branch lengths, scale="max" div-by-zero guard
## ----------------------------------------------------------------
star <- read.tree(text = "((A,B,C,D,E));")
star$edge.length <- rep(0, nrow(star$edge))
tryCatch({
  r <- Average(list(star, star), method = "nj", scale = "max")
  PASS(sprintf("star+zero-bl scale=max: class=%s", class(r)))
}, error = function(e) {
  FAIL("star+zero-bl scale=max", conditionMessage(e))
})

## ----------------------------------------------------------------
## P4: check.labels=FALSE with MISMATCHED tip sets
## P4a: different label (subset / superset) — silent wrong result?
## ----------------------------------------------------------------
good <- read.tree(text = "((A,B),(C,D));")
bad  <- read.tree(text = "((A,B),(C,E));")
tryCatch({
  r <- Average(list(good, bad), check.labels = FALSE)
  FAIL("P4a mismatched-no-check: should have errored or produced NA",
       sprintf("got %d tips: %s", length(r$tip.label),
               paste(sort(r$tip.label), collapse = ",")))
}, error = function(e) {
  INFO("P4a mismatched-no-check errors", conditionMessage(e))
})

## ----------------------------------------------------------------
## P5: label permutation — does reordering tips change the result?
##     (cophenetic indexes by name so it should NOT)
## ----------------------------------------------------------------
set.seed(42)
trees  <- rmtree(5, 8)
perm   <- sample(trees[[1]]$tip.label)
# reorder tip.label in every tree (same tree, different internal ordering)
trees2 <- lapply(trees, function(tr) {
  TreeTools::RenumberTips(tr, perm)
})
r1 <- Average(trees,  method = "nj")
r2 <- Average(trees2, method = "nj")
labs <- sort(r1$tip.label)
d1 <- ape::cophenetic.phylo(r1)[labs, labs]
d2 <- ape::cophenetic.phylo(r2)[labs, labs]
if (isTRUE(all.equal(d1, d2, tolerance = 1e-10))) {
  PASS("P5 permutation-invariance: identical distances")
} else {
  FAIL("P5 permutation-invariance: distances differ after tip relabelling")
}

## ----------------------------------------------------------------
## P6: NA branch lengths in a tree — cophenetic.phylo behaviour
## ----------------------------------------------------------------
tna <- rtree(6)
tna$edge.length[3] <- NA
tryCatch({
  r <- Average(list(tna, rtree(6)), method = "nj")
  dmat <- ape::cophenetic.phylo(ape::unroot(r))
  if (any(is.na(dmat))) {
    FAIL("P6 NA edge.length: NA propagated into result distances")
  } else {
    INFO("P6 NA edge.length: completed without error (NA silently treated as 0 or similar)")
  }
}, error = function(e) {
  INFO("P6 NA edge.length: errored", conditionMessage(e))
})

## ----------------------------------------------------------------
## P7: negative branch lengths in input
## ----------------------------------------------------------------
tneg <- rtree(6)
tneg$edge.length[1] <- -0.5
tryCatch({
  r <- Average(list(tneg, rtree(6)), method = "nj")
  dmat <- ape::cophenetic.phylo(ape::unroot(r))
  if (any(dmat < -1e-10)) {
    FAIL("P7 negative-bl: negative output distances")
  } else {
    INFO("P7 negative-bl: completed, no negative output distances")
  }
}, error = function(e) {
  INFO("P7 negative-bl: errored", conditionMessage(e))
})

## ----------------------------------------------------------------
## P8: non-logical edgeLengths — e.g. "yes", c(TRUE,FALSE)
## ----------------------------------------------------------------
trees <- rmtree(3, 6)
for (bad_arg in list("yes", c(TRUE, FALSE), 1L)) {
  tryCatch({
    r <- Average(trees, edgeLengths = bad_arg)
    INFO(sprintf("P8 edgeLengths=%s: silently treated as %s",
                 deparse(bad_arg), deparse(r$Nnode)))
  }, error   = function(e) INFO(sprintf("P8 edgeLengths=%s errors", deparse(bad_arg)),
                                conditionMessage(e)),
     warning = function(w) INFO(sprintf("P8 edgeLengths=%s warns", deparse(bad_arg)),
                                conditionMessage(w)))
}

## ----------------------------------------------------------------
## P9: weights = c(1,0) with zero-weight tree — additive matrix?
## ----------------------------------------------------------------
t1 <- read.tree(text = "((((A,X),B),C),D);")
t2 <- read.tree(text = "(((A,B),C),(D,X));")
for (m in c("nj", "fastme.bal")) {
  r <- Average(list(t1, t2), method = m, weights = c(1, 0))
  # distance matrix should equal t1's exactly
  labs <- sort(t1$tip.label)
  d_ref  <- ape::cophenetic.phylo(t1)[labs, labs]
  d_got  <- ape::cophenetic.phylo(r)[labs, labs]
  rss    <- sum((d_ref[lower.tri(d_ref)] - d_got[lower.tri(d_got)])^2)
  if (rss < 1e-8) {
    PASS(sprintf("P9 weight=(1,0) method=%s: distances match t1 (rss=%.2e)", m, rss))
  } else {
    INFO(sprintf("P9 weight=(1,0) method=%s: RSS=%.4f (approx fit, not exact)", m, rss))
  }
}

## ----------------------------------------------------------------
## P10: all-identical trees (additive average) — LS should recover exactly
## ----------------------------------------------------------------
has_ls <- requireNamespace("TreeSearch", quietly = TRUE) &&
  exists("LeastSquaresTree", where = asNamespace("TreeSearch"), mode = "function")
if (has_ls) {
  set.seed(7)
  tr <- rtree(7)
  avg <- Average(list(tr, tr, tr), method = "ls")
  labs <- sort(tr$tip.label)
  d_ref <- ape::cophenetic.phylo(tr)[labs, labs]
  d_got <- ape::cophenetic.phylo(avg)[labs, labs]
  rss <- sum((d_ref[lower.tri(d_ref)] - d_got[lower.tri(d_got)])^2)
  if (rss < 1e-8) PASS(sprintf("P10 identical-trees ls: RSS=%.2e", rss))
  else FAIL("P10 identical-trees ls: RSS too large", sprintf("%.6f", rss))
} else {
  INFO("P10 ls skipped (LeastSquaresTree not available)")
}

## ----------------------------------------------------------------
## P11: single tree returned unrooted — verify no refit distortion
## ----------------------------------------------------------------
set.seed(1)
tr <- rtree(8)
r  <- Average(list(tr))
if (!ape::is.rooted(r)) {
  PASS("P11 single-tree unrooted")
} else {
  FAIL("P11 single-tree still rooted")
}
labs <- sort(tr$tip.label)
d_orig <- ape::cophenetic.phylo(ape::unroot(tr))[labs, labs]
d_ret  <- ape::cophenetic.phylo(r)[labs, labs]
rss <- sum((d_orig[lower.tri(d_orig)] - d_ret[lower.tri(d_ret)])^2)
if (rss < 1e-10) {
  PASS(sprintf("P11 single-tree distances preserved (rss=%.2e)", rss))
} else {
  FAIL("P11 single-tree distances distorted", sprintf("rss=%.6f", rss))
}

## ----------------------------------------------------------------
## P12: data.frame input — should error
## ----------------------------------------------------------------
tryCatch({
  r <- Average(data.frame(a=1))
  FAIL("P12 data.frame input: should error")
}, error = function(e) {
  PASS(sprintf("P12 data.frame errors: %s", conditionMessage(e)))
})

## ----------------------------------------------------------------
## P13: weights all zero — should error
## ----------------------------------------------------------------
tryCatch({
  r <- Average(rmtree(2, 6), weights = c(0, 0))
  FAIL("P13 zero-weights: should error")
}, error = function(e) {
  PASS(sprintf("P13 zero-weights errors: %s", conditionMessage(e)))
})

## ----------------------------------------------------------------
## P14: single-entry multiPhylo (length-1 list)
## ----------------------------------------------------------------
tr <- rtree(5)
r  <- Average(list(tr))
if (inherits(r, "phylo") && !ape::is.rooted(r)) {
  PASS("P14 length-1 list returns unrooted phylo")
} else {
  FAIL("P14 length-1 list: unexpected result", class(r))
}

## ----------------------------------------------------------------
## P15: edgeLengths=FALSE ignores branch lengths in avg — topology-only
## ----------------------------------------------------------------
set.seed(3)
t1bl <- rtree(6)     # has branch lengths
t1np <- t1bl
t1np$edge.length <- NULL   # no branch lengths
# with edgeLengths=FALSE both should give same result as all-unit lengths
r_bl <- Average(list(t1bl, t1bl), method = "nj", edgeLengths = FALSE)
r_no <- Average(list(t1np, t1np), method = "nj", edgeLengths = FALSE)
labs <- sort(t1bl$tip.label)
d1 <- ape::cophenetic.phylo(r_bl)[labs, labs]
d2 <- ape::cophenetic.phylo(r_no)[labs, labs]
if (isTRUE(all.equal(d1, d2, tolerance = 1e-10))) {
  PASS("P15 edgeLengths=FALSE: same as no-lengths version")
} else {
  FAIL("P15 edgeLengths=FALSE differs between with/without lengths input")
}

## ----------------------------------------------------------------
## P16: check.labels=FALSE with permuted-order tips — is RenumberTips skipped?
## If RenumberTips is skipped when check.labels=FALSE, the cophenetic index
## [labs, labs] still works because cophenetic returns a named matrix.
## Verify the result is the same as with check.labels=TRUE.
## ----------------------------------------------------------------
set.seed(99)
base_trees <- rmtree(4, 8)
perm_trees  <- lapply(base_trees, function(tr) TreeTools::RenumberTips(tr, sample(tr$tip.label)))
r_check    <- Average(base_trees,  method = "nj", check.labels = TRUE)
r_nocheck  <- Average(perm_trees,  method = "nj", check.labels = FALSE)
labs <- sort(r_check$tip.label)
d_c  <- ape::cophenetic.phylo(r_check)[labs,   labs]
d_nc <- ape::cophenetic.phylo(r_nocheck)[labs, labs]
if (isTRUE(all.equal(d_c, d_nc, tolerance = 1e-8))) {
  PASS("P16 check.labels=FALSE with permuted tips: same result as TRUE")
} else {
  FAIL("P16 check.labels=FALSE with permuted tips: different result (label-order bug?)")
}

cat("\n--- done ---\n")
