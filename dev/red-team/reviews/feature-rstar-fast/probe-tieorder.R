# Stress the tie-order: std::sort is NOT stable, so equal-weight MST edges may be
# processed in ANY relative order. Does ANY tie-order cause the MST single-linkage
# family to miss a threshold component? We randomly permute within equal-weight
# groups (and also randomise Prim's choice among equal-best frontier nodes, which
# changes WHICH MST is built when s has ties -- the MST is not unique under ties).

ufFind <- function(f, x) { while (f[x] != x) { f[x] <- f[f[x]]; x <- f[x] }; x }

oldCandidates <- function(s, n) {
  cand <- list()
  for (theta in 1:(n - 2)) {
    f <- seq_len(n)
    for (a in 1:(n - 1)) for (b in (a + 1):n)
      if (s[a, b] >= theta) { ra <- ufFind(f, a); rb <- ufFind(f, b); if (ra != rb) f[ra] <- rb }
    comp <- split(seq_len(n), vapply(seq_len(n), function(x) ufFind(f, x), 1L))
    for (cc in comp) if (length(cc) >= 2 && length(cc) <= n - 1) cand[[length(cand) + 1L]] <- sort(cc)
  }
  unique(cand)
}

# Prim with RANDOM tie-break on frontier choice, then RANDOM within-weight edge order.
newCandidatesRandTie <- function(s, n) {
  bestW <- rep(-1L, n); bestTo <- rep(-1L, n); inTree <- rep(FALSE, n)
  mst <- list(); inTree[1] <- TRUE
  for (v in 2:n) { bestW[v] <- s[1, v]; bestTo[v] <- 1 }
  for (iter in 2:n) {
    cand <- which(!inTree)
    mw <- max(bestW[cand]); ties <- cand[bestW[cand] == mw]; u <- if (length(ties) == 1) ties else sample(ties, 1)
    inTree[u] <- TRUE
    mst[[length(mst) + 1L]] <- c(bestW[u], bestTo[u], u)
    for (v in 1:n) if (!inTree[v]) { sw <- s[u, v]; if (sw > bestW[v]) { bestW[v] <- sw; bestTo[v] <- u } }
  }
  w <- vapply(mst, `[`, 1, 1)
  # random permutation that still sorts by decreasing weight (shuffle ties)
  ord <- order(w, runif(length(w)), decreasing = TRUE)
  mst <- mst[ord]
  dpar <- seq_len(n); members <- lapply(seq_len(n), function(i) i); cand <- list()
  for (e in mst) {
    ru <- ufFind(dpar, e[2]); rv <- ufFind(dpar, e[3]); if (ru == rv) next
    A <- sort(c(members[[ru]], members[[rv]]))
    if (length(A) >= 2 && length(A) <= n - 1) cand[[length(cand) + 1L]] <- A
    if (length(members[[ru]]) >= length(members[[rv]])) { dpar[rv] <- ru; members[[ru]] <- A; members[[rv]] <- integer(0) }
    else { dpar[ru] <- rv; members[[rv]] <- A; members[[ru]] <- integer(0) }
  }
  unique(cand)
}
key <- function(lst) sort(vapply(lst, paste, "", collapse = ","))

set.seed(11)
missingTotal <- 0L; wc <- NULL
for (trial in seq_len(40000)) {
  n <- sample(4:9, 1); maxv <- sample(1:3, 1)
  s <- matrix(0L, n, n)
  for (a in 1:(n - 1)) for (b in (a + 1):n) { v <- sample(0:maxv, 1); s[a, b] <- v; s[b, a] <- v }
  old <- key(oldCandidates(s, n)); new <- key(newCandidatesRandTie(s, n))
  miss <- setdiff(old, new)
  if (length(miss)) { missingTotal <- missingTotal + 1L; if (is.null(wc)) wc <- list(s=s,n=n,miss=miss,old=old,new=new) }
}
cat("RANDOM-tie trials missing a threshold component:", missingTotal, "/ 40000\n")
if (!is.null(wc)) { cat("n=",wc$n,"\n"); print(wc$s); cat("MISSING:",paste(wc$miss,collapse=" | "),"\n") }
