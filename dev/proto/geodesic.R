# Prototype of the Owen-Provan (2011) GTP geodesic distance, to validate the
# algorithm against the hand-derived 15*sqrt(2) oracle BEFORE porting to C++.
# Representation: a tree -> list of clades (integer tip-sets, canonicalised to
# EXCLUDE the reference tip = column 1) with lengths, plus per-tip leaf lengths.

suppressMessages(library(TreeTools))

treeToBHV <- function(tree, tipLabels) {
  sp <- as.Splits(tree, tipLabels = tipLabels)
  M <- as.logical(sp)
  if (is.null(dim(M))) M <- matrix(M, nrow = 1,
                                   dimnames = list(rownames(sp), NULL))
  nodeOfSplit <- as.integer(rownames(sp))
  ei <- match(nodeOfSplit, tree$edge[, 2])
  lens <- tree$edge.length[ei]
  ref <- 1L
  clades <- lapply(seq_len(nrow(M)), function(i) {
    row <- as.logical(M[i, ])
    if (row[ref]) row <- !row
    which(row)
  })
  leaf <- numeric(length(tipLabels))
  for (i in seq_along(tipLabels)) {
    er <- match(i, tree$edge[, 2])
    leaf[i] <- if (is.na(er)) 0 else tree$edge.length[er]
  }
  list(clades = clades, lens = lens, leaf = leaf, nTip = length(tipLabels))
}

.subset <- function(a, b) all(a %in% b)                  # a subseteq b
.disjoint <- function(a, b) length(intersect(a, b)) == 0L
compatible <- function(a, b) .subset(a, b) || .subset(b, a) || .disjoint(a, b)
clEqual <- function(a, b) length(a) == length(b) && all(a %in% b)

# Min-weight vertex cover of bipartite incompatibility graph via max-flow.
# incid[i, j] TRUE iff a_i incompatible with b_j (an edge to cover).
# Returns logical aCover (length na), bCover (length nb), and cover weight
# (with weights normalised so sum(wA) = sum(wB) = 1).
minWeightVC <- function(incid, wA, wB) {
  na <- length(wA); nb <- length(wB)
  sA <- sum(wA); sB <- sum(wB)
  nwA <- wA / sA; nwB <- wB / sB
  # Node ids: source = 1, a_i = 1 + i, b_j = 1 + na + j, sink = 2 + na + nb
  S <- 1L; T <- 2L + na + nb
  N <- T
  INF <- 1e18
  cap <- matrix(0, N, N)
  for (i in seq_len(na)) cap[S, 1L + i] <- nwA[i]
  for (j in seq_len(nb)) cap[1L + na + j, T] <- nwB[j]
  for (i in seq_len(na)) for (j in seq_len(nb)) {
    if (isTRUE(incid[i, j])) cap[1L + i, 1L + na + j] <- INF
  }
  flow <- matrix(0, N, N)
  eps <- 1e-12
  repeat {
    # BFS for augmenting path in residual graph
    pred <- integer(N); pred[] <- NA_integer_; pred[S] <- S
    q <- c(S); found <- FALSE
    while (length(q)) {
      u <- q[1]; q <- q[-1]
      if (u == T) { found <- TRUE; break }
      for (v in seq_len(N)) {
        if (is.na(pred[v]) && cap[u, v] - flow[u, v] > eps) {
          pred[v] <- u; q <- c(q, v)
        }
      }
    }
    if (!found) break
    # bottleneck
    v <- T; b <- INF
    while (v != S) { u <- pred[v]; b <- min(b, cap[u, v] - flow[u, v]); v <- u }
    v <- T
    while (v != S) { u <- pred[v]; flow[u, v] <- flow[u, v] + b
                     flow[v, u] <- flow[v, u] - b; v <- u }
  }
  # reachable set from S in residual graph
  reach <- logical(N); reach[S] <- TRUE; q <- c(S)
  while (length(q)) {
    u <- q[1]; q <- q[-1]
    for (v in seq_len(N)) {
      if (!reach[v] && cap[u, v] - flow[u, v] > eps) { reach[v] <- TRUE; q <- c(q, v) }
    }
  }
  aCover <- !reach[1L + seq_len(na)]          # A in cover iff NOT reachable
  bCover <-  reach[1L + na + seq_len(nb)]      # B in cover iff reachable
  weight <- sum(nwA[aCover]) + sum(nwB[bCover])
  list(aCover = aCover, bCover = bCover, weight = weight)
}

vnorm <- function(x) sqrt(sum(x^2))
TOL <- 1e-10

# GTP on a group of mutually non-common edges. Returns list of ratios; each
# ratio = list(eC, eL, fC, fL) (clades + lengths on each side).
gtpNoCommon <- function(Ac, Al, Bc, Bl) {
  rs <- list()
  queue <- list(list(a = seq_along(Ac), b = seq_along(Bc)))
  while (length(queue)) {
    r <- queue[[1]]; queue[[1]] <- NULL
    na <- length(r$a); nb <- length(r$b)
    if (na == 0 || nb == 0) { rs[[length(rs) + 1]] <- r; next }
    incid <- matrix(FALSE, na, nb)
    for (ii in seq_len(na)) for (jj in seq_len(nb))
      incid[ii, jj] <- !compatible(Ac[[r$a[ii]]], Bc[[r$b[jj]]])
    cov <- minWeightVC(incid, Al[r$a]^2, Bl[r$b]^2)
    # split iff the min-weight vertex cover has weight < 1 (Owen & Provan 2011,
    # Lemma 3.2); a weight-0 cover correctly separates compatible cross-edges.
    if (cov$weight >= 1 - TOL) {
      rs[[length(rs) + 1]] <- r; next
    }
    a1 <- r$a[cov$aCover];  b1 <- r$b[!cov$bCover]
    a2 <- r$a[!cov$aCover]; b2 <- r$b[cov$bCover]
    queue <- c(list(list(a = a1, b = b1), list(a = a2, b = b2)), queue)
  }
  lapply(rs, function(r) list(eC = Ac[r$a], eL = Al[r$a], fC = Bc[r$b], fL = Bl[r$b]))
}

# Recursively partition non-common edges into groups separated by common edges.
splitOnCommon <- function(Ac, Al, Bc, Bl) {
  ci <- cj <- NA_integer_
  for (i in seq_along(Ac)) { for (j in seq_along(Bc))
    if (clEqual(Ac[[i]], Bc[[j]])) { ci <- i; cj <- j; break }
    if (!is.na(ci)) break }
  if (is.na(ci)) {
    if (length(Ac) || length(Bc))
      return(list(list(Ac = Ac, Al = Al, Bc = Bc, Bl = Bl)))
    return(list())
  }
  cl <- Ac[[ci]]
  belowA <- vapply(Ac, function(s) length(s) < length(cl) && .subset(s, cl), logical(1))
  eqA    <- vapply(Ac, function(s) clEqual(s, cl), logical(1))
  belowB <- vapply(Bc, function(s) length(s) < length(cl) && .subset(s, cl), logical(1))
  eqB    <- vapply(Bc, function(s) clEqual(s, cl), logical(1))
  aboveA <- !belowA & !eqA
  aboveB <- !belowB & !eqB
  c(splitOnCommon(Ac[belowA], Al[belowA], Bc[belowB], Bl[belowB]),
    splitOnCommon(Ac[aboveA], Al[aboveA], Bc[aboveB], Bl[aboveB]))
}

geoDist <- function(A, B) {
  leafSq <- sum((A$leaf - B$leaf)^2)
  ac <- A$clades; al <- A$lens; bc <- B$clades; bl <- B$lens
  commonSq <- 0; matchedB <- rep(FALSE, length(bc)); keepA <- rep(TRUE, length(ac))
  for (i in seq_along(ac)) for (j in seq_along(bc))
    if (!matchedB[j] && clEqual(ac[[i]], bc[[j]])) {
      commonSq <- commonSq + (al[i] - bl[j])^2
      keepA[i] <- FALSE; matchedB[j] <- TRUE; break
    }
  groups <- splitOnCommon(ac[keepA], al[keepA], bc[!matchedB], bl[!matchedB])
  ratios <- list()
  for (g in groups) ratios <- c(ratios, gtpNoCommon(g$Ac, g$Al, g$Bc, g$Bl))
  # ratio value e/f; combine to non-descending then sum (e+f)^2
  rv <- vapply(ratios, function(r) {
    e <- vnorm(r$eL); f <- vnorm(r$fL)
    if (f == 0) Inf else if (e == 0) 0 else e / f
  }, numeric(1))
  ratios <- ratios[order(rv)]
  rsSq <- sum(vapply(ratios, function(r) (vnorm(r$eL) + vnorm(r$fL))^2, numeric(1)))
  sqrt(rsSq + commonSq + leafSq)
}

BHVDist <- function(t1, t2) {
  tl <- t1$tip.label
  geoDist(treeToBHV(t1, tl), treeToBHV(t2, tl))
}

# ---- Tests --------------------------------------------------------------
rt <- function(txt) ape::read.tree(text = txt)
T  <- rt("(0:0,(((1:0,2:0):4,(3:0,4:0):10):3,5:0):0);")   # OP Fig 1 T
Tp <- rt("(0:0,(1:0,((2:0,3:0):4,(4:0,5:0):3):10):0);")   # OP Fig 1 T'
cat(sprintf("ORACLE  d = %.10f   15*sqrt2 = %.10f   %s\n",
            BHVDist(T, Tp), 15 * sqrt(2),
            ifelse(abs(BHVDist(T, Tp) - 15 * sqrt(2)) < 1e-9, "PASS", "FAIL")))
cat(sprintf("SYMMETRY d(Tp,T) = %.10f\n", BHVDist(Tp, T)))

# Euclidean (same topology)
C1 <- rt("(0:0,(((1:0,2:0):4,(3:0,4:0):10):3,5:0):0);")
C2 <- rt("(0:0,(((1:0,2:0):1,(3:0,4:0):2):3,5:0):0);")
cat(sprintf("EUCLID  d = %.10f   sqrt(9+64) = %.10f   %s\n",
            BHVDist(C1, C2), sqrt(73),
            ifelse(abs(BHVDist(C1, C2) - sqrt(73)) < 1e-9, "PASS", "FAIL")))

# Identical trees -> 0
cat(sprintf("IDENT   d = %.10f   %s\n", BHVDist(T, T),
            ifelse(abs(BHVDist(T, T)) < 1e-9, "PASS", "FAIL")))

# Scaling: d(cT, cT') = c d(T,T')
Ts <- T; Ts$edge.length <- T$edge.length * 3
Tps <- Tp; Tps$edge.length <- Tp$edge.length * 3
cat(sprintf("SCALING d = %.10f   3*15sqrt2 = %.10f   %s\n",
            BHVDist(Ts, Tps), 3 * 15 * sqrt(2),
            ifelse(abs(BHVDist(Ts, Tps) - 3 * 15 * sqrt(2)) < 1e-9, "PASS", "FAIL")))

# Cone (single incompatible pair, otherwise common): leaves 0..5,
# both share clade {1,2,3,4} length 5; differ {1,2}:3 (T) vs {2,3}:4 (Tp), incompatible
K1 <- rt("(0:0,((((1:0,2:0):3,3:0):0,4:0):5,5:0):0);")  # clades {1,2}:3,{1,2,3}:0,{1,2,3,4}:5
K2 <- rt("(0:0,(((1:0,(2:0,3:0):4):0,4:0):5,5:0):0);")  # clades {2,3}:4,{1,2,3}:0,{1,2,3,4}:5
cat(sprintf("CONE-ish d = %.10f\n", BHVDist(K1, K2)))
