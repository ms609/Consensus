#include <Rcpp.h>
#include <vector>
#include <string>
#include <numeric>    // std::iota
#include <algorithm>  // std::max, std::sort
#include <set>

// =============================================================================
// R* consensus tree.
//
// Definition (Jansson, Sung, Vu & Yiu 2016, Algorithmica 76:1224-1244, after
// Degnan et al. 2009): given k rooted trees on the same leaf set L (n = |L|),
//   * #ab|c = number of input trees in which the resolved triplet ab|c is
//     consistent (fan triplets count for NOTHING -- they have no impact);
//   * the majority resolved-triplet set is
//        R_maj = { ab|c : #ab|c > max(#ac|b, #bc|a) }   (strict PLURALITY);
//   * the R* consensus tree is the unique tree tau with r(tau) subset of R_maj
//     maximising the number of internal nodes.
//
// Lemma 1.1 (Jansson et al. 2016) gives the construction used here:
//   the R* tree always exists, is unique, and its clusters are EXACTLY the
//   STRONG CLUSTERS of R_maj (plus L and the singletons), where A subset of L is
//   a strong cluster iff  aa'|x in R_maj  for every pair a,a' in A and every
//   x not in A (the "for-all-outgroup" rule).
//
// This file therefore:
//   1. Tallies resolved triplet states across all trees           -- O(k n^3).
//      (LCA-depth logic copied verbatim from local_consensus.cpp, the
//      38/38-oracle-validated MinRLC/MinILC port, which is left untouched.)
//   2. Selects R_maj by strict plurality per 3-taxon subset.
//   3. ASSEMBLES via strong clusters (NOT Aho-BUILD): generate candidate
//      clusters as the single-linkage threshold components of the similarity
//      s(a,b) = #{w : ab|w in R_maj} -- a provable superset of the strong
//      clusters (strong subset Apresjan subset threshold-components) -- test
//      each candidate against the strong-cluster rule, then nest the (laminar)
//      survivors into the tree.
//
// The strong clusters are guaranteed laminar (two overlapping strong clusters
// would force two different resolutions of one triple into R_maj); we assert
// this as a defensive bug-catch.
//
// Implementation is correctness-first and polynomial (O(k n^3) tally + O(n^4)
// assembly).  The O(n^2)/O(n^2 polylog) Apresjan-hierarchy algorithm of Jansson
// et al. (2016) is a deferred speed optimisation.
//
// Taxa are 0-indexed internally; the ape edge matrix is 1-indexed (tip node v
// <-> taxon v-1).
//
// AUTHOR-CONFIRMED behaviour (Jansson et al. 2016, Sect. 1.1 / Lemma 1.1):
//   * fans abstain (they "have no impact" on R_maj);
//   * ties leave a triple out of R_maj;
//   * the construction is the strong-cluster rule -- there is no BUILD-failure /
//     collapse ambiguity (the tree always exists and is unique).
// =============================================================================

// Memory guard for the dense n^3 tally tensor.  n = 200 -> 200^3 ints ~= 32 MB.
// This is a MEMORY safeguard; R* is polynomial (no exponential blow-up).
static const int RSTAR_MAX_TIP = 200;

// ---------------------------------------------------------------------------
// parent[] and depth[] for a single tree (verbatim from local_consensus.cpp).
// Preorder is guaranteed by the R wrapper (Preorder() is called before
// extracting edge matrices); a single forward pass suffices.
// ---------------------------------------------------------------------------
static void buildParentDepth(const Rcpp::IntegerMatrix& edge, int nNode,
                             std::vector<int>& parent, std::vector<int>& depth) {
  parent.assign(nNode + 1, -1);
  depth.assign(nNode + 1, 0);
  int nRow = edge.nrow();
  for (int r = 0; r < nRow; r++) {
    int p = edge(r, 0);
    int c = edge(r, 1);
    parent[c] = p;
    depth[c] = depth[p] + 1;
  }
}

static int lcaDepth(int u, int v,
                    const std::vector<int>& parent,
                    const std::vector<int>& depth) {
  while (depth[u] > depth[v]) u = parent[u];
  while (depth[v] > depth[u]) v = parent[v];
  while (u != v) { u = parent[u]; v = parent[v]; }
  return depth[u];
}

// Flat index into the dense tensor.  Convention (matching local_consensus.cpp):
// tri[a][b][c] with a < b counts the resolved triplet whose CLOSE pair is {a,b}
// and whose OUTGROUP is c.
static inline size_t triIdx(int a, int b, int c, int n) {
  return ((size_t)a * n + b) * n + c;
}

// Is the resolved triplet with close pair {p,q} and outgroup x in R_maj?
// (R_maj is encoded as the non-zero cells of `tri` after plurality selection.)
static inline bool inRmaj(const std::vector<int>& tri, int p, int q, int x, int n) {
  int lo = p < q ? p : q;
  int hi = p < q ? q : p;
  return tri[triIdx(lo, hi, x, n)] > 0;
}

// Disjoint-set find with path halving.
static int dsuFind(std::vector<int>& f, int x) {
  while (f[x] != x) { f[x] = f[f[x]]; x = f[x]; }
  return x;
}

// ---------------------------------------------------------------------------
// Assemble the R* tree from R_maj (encoded in `tri`) via strong clusters.
// Returns a Newick string (1-indexed integer labels, no trailing ';').
// ---------------------------------------------------------------------------
static std::string assembleRStar(const std::vector<int>& tri, int n) {
  // --- similarity s(a,b) = #{w : ab|w in R_maj} -------------------------------
  std::vector<int> s((size_t)n * n, 0);
  for (int a = 0; a < n; a++) {
    for (int b = a + 1; b < n; b++) {
      int cnt = 0;
      for (int w = 0; w < n; w++) {
        if (w != a && w != b && inRmaj(tri, a, b, w, n)) cnt++;
      }
      s[(size_t)a * n + b] = cnt;
      s[(size_t)b * n + a] = cnt;
    }
  }

  // --- candidate clusters: single-linkage threshold components of s -----------
  // For each integer threshold theta >= 1, the connected components of
  // {(a,b): s(a,b) >= theta} of size >= 2.  This is a provable superset of the
  // strong clusters (strong subset Apresjan subset threshold-components).
  std::set<std::vector<int>> candidates;
  std::vector<int> f(n);
  for (int theta = 1; theta <= n - 2; theta++) {
    Rcpp::checkUserInterrupt();
    std::iota(f.begin(), f.end(), 0);
    for (int a = 0; a < n; a++)
      for (int b = a + 1; b < n; b++)
        if (s[(size_t)a * n + b] >= theta) {
          int ra = dsuFind(f, a), rb = dsuFind(f, b);
          if (ra != rb) f[ra] = rb;
        }
    // Collect components.
    std::vector<std::vector<int>> comp(n);
    for (int x = 0; x < n; x++) comp[dsuFind(f, x)].push_back(x);
    for (int x = 0; x < n; x++) {
      if ((int)comp[x].size() >= 2 && (int)comp[x].size() <= n - 1) {
        candidates.insert(comp[x]);  // already ascending
      }
    }
  }

  // --- keep candidates that are strong clusters -------------------------------
  // strong(A): for all pairs a,a' in A and all x not in A, aa'|x in R_maj.
  std::vector<std::vector<int>> kept;
  std::vector<char> inA(n);
  for (std::set<std::vector<int>>::const_iterator it = candidates.begin();
       it != candidates.end(); ++it) {
    Rcpp::checkUserInterrupt();
    const std::vector<int>& A = *it;
    std::fill(inA.begin(), inA.end(), (char)0);
    for (int v : A) inA[v] = 1;
    bool strong = true;
    for (size_t i = 0; i < A.size() && strong; i++) {
      for (size_t j = i + 1; j < A.size() && strong; j++) {
        for (int x = 0; x < n && strong; x++) {
          if (!inA[x] && !inRmaj(tri, A[i], A[j], x, n)) strong = false;
        }
      }
    }
    if (strong) kept.push_back(A);
  }

  // --- defensive laminarity check ---------------------------------------------
  // Strong clusters are provably laminar; a violation indicates a bug.
  for (size_t i = 0; i < kept.size(); i++) {
    std::vector<char> mi(n, 0);
    for (int v : kept[i]) mi[v] = 1;
    for (size_t j = i + 1; j < kept.size(); j++) {
      int inter = 0, onlyJ = 0;
      for (int v : kept[j]) { if (mi[v]) inter++; else onlyJ++; }
      int onlyI = (int)kept[i].size() - inter;
      // Overlap with each having a private element => not nested, not disjoint.
      // # nocov start
      // Unreachable: the strong clusters of R_maj form a laminar family
      // (Jansson et al. 2016, Lemma 1.1), so no two overlap partially.
      // Defensive internal-consistency guard.
      if (inter > 0 && onlyI > 0 && onlyJ > 0) {
        Rcpp::stop("rStarConsensus: non-laminar strong clusters (internal bug).");
      }
      // # nocov end
    }
  }

  // --- assemble the laminar family into a tree --------------------------------
  // Nodes: every leaf singleton, every kept strong cluster, and L (the root).
  std::vector<std::vector<int>> nodes;          // member lists (ascending)
  for (int i = 0; i < n; i++) nodes.push_back(std::vector<int>(1, i));  // singletons
  for (const std::vector<int>& A : kept) nodes.push_back(A);
  std::vector<int> all(n); std::iota(all.begin(), all.end(), 0);
  nodes.push_back(all);                         // L = root (last node)
  int nNode = (int)nodes.size();
  int rootIdx = nNode - 1;

  // masks for subset tests
  std::vector<std::vector<char>> mask(nNode, std::vector<char>(n, 0));
  for (int i = 0; i < nNode; i++) for (int v : nodes[i]) mask[i][v] = 1;

  // parent[i] = smallest node strictly containing node i (unique by laminarity)
  std::vector<int> parent(nNode, -1);
  for (int i = 0; i < nNode; i++) {
    if (i == rootIdx) continue;
    int best = -1;
    for (int j = 0; j < nNode; j++) {
      if (j == i) continue;
      if (nodes[j].size() <= nodes[i].size()) continue;
      bool subset = true;
      for (int v : nodes[i]) if (!mask[j][v]) { subset = false; break; }
      if (!subset) continue;
      if (best == -1 || nodes[j].size() < nodes[best].size()) best = j;
    }
    parent[i] = best;  // root's children point at rootIdx
  }

  std::vector<std::vector<int>> children(nNode);
  for (int i = 0; i < nNode; i++) if (parent[i] != -1) children[parent[i]].push_back(i);

  // recursive Newick
  std::vector<std::string> memo(nNode);
  // iterative post-order to avoid deep recursion worries (n <= 200 is fine, but
  // this keeps it simple and stack-safe)
  std::vector<int> order;
  {
    std::vector<int> stack(1, rootIdx);
    while (!stack.empty()) {
      int u = stack.back(); stack.pop_back();
      order.push_back(u);
      for (int c : children[u]) stack.push_back(c);
    }
  }
  for (int idx = (int)order.size() - 1; idx >= 0; idx--) {
    int u = order[idx];
    if (children[u].empty()) {
      // a leaf (singleton); label is member+1
      memo[u] = std::to_string(nodes[u][0] + 1);
    } else {
      std::string out = "(";
      for (size_t c = 0; c < children[u].size(); c++) {
        if (c) out += ",";
        out += memo[children[u][c]];
      }
      out += ")";
      memo[u] = out;
    }
  }
  return memo[rootIdx];
}

// ---------------------------------------------------------------------------
// Main Rcpp entry point.  Returns the Newick string (no trailing ';').
// [[Rcpp::export]]
std::string rStarConsensus(Rcpp::List edgeList, int nTip) {
  int n = nTip;
  if (n > RSTAR_MAX_TIP) {
    Rcpp::stop("rStarConsensus: n = %d exceeds the memory guard of %d leaves "
               "(dense n^3 triplet tensor).", n, RSTAR_MAX_TIP);
  }
  int nTree = edgeList.size();

  std::vector<int> tri((size_t)n * n * n, 0);

  // ---- Step 1: tally resolved triplet states across every tree ----
  for (int t = 0; t < nTree; t++) {
    Rcpp::checkUserInterrupt();
    Rcpp::IntegerMatrix edge = edgeList[t];
    int nNode = 0;
    for (int r = 0; r < edge.nrow(); r++) {
      if (edge(r, 0) > nNode) nNode = edge(r, 0);
      if (edge(r, 1) > nNode) nNode = edge(r, 1);
    }
    std::vector<int> par, dep;
    buildParentDepth(edge, nNode, par, dep);

    // Per-tree pairwise LCA-depth matrix (taxon i <-> node i+1) -> O(1) lookups.
    std::vector<int> D((size_t)n * n, 0);
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        int d = lcaDepth(i + 1, j + 1, par, dep);
        D[(size_t)i * n + j] = d;
        D[(size_t)j * n + i] = d;
      }
    }

    // Triplet identification (verbatim from local_consensus.cpp).
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        int dij = D[(size_t)i * n + j];
        for (int k = j + 1; k < n; k++) {
          int dik = D[(size_t)i * n + k];
          int djk = D[(size_t)j * n + k];
          if (dij == dik && dij == djk) continue;          // fan -> abstain
          if (dik == djk && dij > dik)       tri[triIdx(i, j, k, n)]++;  // ij|k
          else if (dij == djk && dik > dij)  tri[triIdx(i, k, j, n)]++;  // ik|j
          else if (dij == dik && djk > dij)  tri[triIdx(j, k, i, n)]++;  // jk|i
        }
      }
    }
  }

  // ---- Step 2: plurality selection (build R_maj) ----
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int k = j + 1; k < n; k++) {
        size_t a = triIdx(i, j, k, n);   // ij|k
        size_t b = triIdx(i, k, j, n);   // ik|j
        size_t c = triIdx(j, k, i, n);   // jk|i
        int c1 = tri[a], c2 = tri[b], c3 = tri[c];
        int mx = std::max(c1, std::max(c2, c3));
        int nMax = (int)(c1 == mx) + (int)(c2 == mx) + (int)(c3 == mx);
        if (mx > 0 && nMax == 1) {
          tri[a] = (c1 == mx) ? c1 : 0;
          tri[b] = (c2 == mx) ? c2 : 0;
          tri[c] = (c3 == mx) ? c3 : 0;
        } else {
          tri[a] = tri[b] = tri[c] = 0;   // tie / all-fan -> not in R_maj
        }
      }
    }
  }

  // ---- Step 3: assemble via strong clusters ----
  return assembleRStar(tri, n);
}
