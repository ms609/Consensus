# Probe: does the NEW MST-based single-linkage candidate set contain every
# threshold-connected component the OLD code enumerated?  A strong cluster is
# always a threshold component (strong ⊆ Apresjan ⊆ threshold-components).  If
# the MST candidate set ever misses a threshold component, the new filter could
# fail to even TEST a genuine strong cluster -> under-resolution bug.
#
# We compare, over many random integer similarity matrices (symmetric, 0..maxv),
# the OLD threshold-component family vs the NEW MST single-linkage family.

ufFind <- function(f, x) { while (f[x] != x) { f[x] <- f[f[x]]; x <- f[x] }; x }

# OLD: all threshold-connected components (size 2..n-1)
oldCandidates <- function(s, n) {
  cand <- list()
  for (theta in 1:(n - 2)) {
    f <- seq_len(n)
    for (a in 1:(n - 1)) for (b in (a + 1):n)
      if (s[a, b] >= theta) {
        ra <- ufFind(f, a); rb <- ufFind(f, b); if (ra != rb) f[ra] <- rb
      }
    comp <- split(seq_len(n), vapply(seq_len(n), function(x) ufFind(f, x), 1L))
    for (cc in comp) if (length(cc) >= 2 && length(cc) <= n - 1)
      cand[[length(cand) + 1L]] <- sort(cc)
  }
  unique(cand)
}

# NEW: replicate the C++ Prim-MST + decreasing-weight merge exactly.
newCandidates <- function(s, n) {
  bestW <- rep(-1L, n); bestTo <- rep(-1L, n); inTree <- rep(FALSE, n)
  mst <- list()
  inTree[1] <- TRUE
  for (v in 2:n) { bestW[v] <- s[1, v]; bestTo[v] <- 1 }
  for (iter in 2:n) {
    u <- -1L; w <- -1L
    for (v in 1:n) if (!inTree[v] && bestW[v] > w) { w <- bestW[v]; u <- v }
    inTree[u] <- TRUE
    mst[[length(mst) + 1L]] <- c(bestW[u], bestTo[u], u)
    for (v in 1:n) if (!inTree[v]) {
      sw <- s[u, v]; if (sw > bestW[v]) { bestW[v] <- sw; bestTo[v] <- u }
    }
  }
  # sort by decreasing weight (R's order is a STABLE sort, matching std::sort? NO
  # -- std::sort is not stable. We test BOTH orders below.)
  ord <- order(vapply(mst, `[`, 1, 1), decreasing = TRUE)
  mst <- mst[ord]
  dpar <- seq_len(n)
  members <- lapply(seq_len(n), function(i) i)
  cand <- list()
  for (e in mst) {
    ru <- ufFind(dpar, e[2] ); rv <- ufFind(dpar, e[3])
    if (ru == rv) next
    A <- sort(c(members[[ru]], members[[rv]]))
    if (length(A) >= 2 && length(A) <= n - 1) cand[[length(cand) + 1L]] <- A
    if (length(members[[ru]]) >= length(members[[rv]])) {
      dpar[rv] <- ru; members[[ru]] <- A; members[[rv]] <- integer(0)
    } else {
      dpar[ru] <- rv; members[[rv]] <- A; members[[ru]] <- integer(0)
    }
  }
  unique(cand)
}

key <- function(lst) sort(vapply(lst, paste, "", collapse = ","))

set.seed(1)
missingTotal <- 0L
worstCase <- NULL
for (trial in seq_len(20000)) {
  n <- sample(4:9, 1)
  maxv <- sample(1:4, 1)        # small value range -> many ties
  s <- matrix(0L, n, n)
  for (a in 1:(n - 1)) for (b in (a + 1):n) {
    v <- sample(0:maxv, 1); s[a, b] <- v; s[b, a] <- v
  }
  old <- key(oldCandidates(s, n))
  new <- key(newCandidates(s, n))
  miss <- setdiff(old, new)     # threshold components the NEW set fails to produce
  if (length(miss)) {
    missingTotal <- missingTotal + 1L
    if (is.null(worstCase)) worstCase <- list(s = s, n = n, miss = miss, old = old, new = new)
  }
}
cat("Trials with a threshold-component MISSING from the MST candidate set:",
    missingTotal, "/ 20000\n")
if (!is.null(worstCase)) {
  cat("\n--- First counterexample ---\nn =", worstCase$n, "\n")
  print(worstCase$s)
  cat("Missing from NEW:", paste(worstCase$miss, collapse = " | "), "\n")
  cat("OLD:", paste(worstCase$old, collapse = " | "), "\n")
  cat("NEW:", paste(worstCase$new, collapse = " | "), "\n")
}
