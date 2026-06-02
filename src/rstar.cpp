#include <Rcpp.h>
#include <vector>
#include <string>
#include <numeric>    // std::iota
#include <algorithm>  // std::max, std::sort, std::merge, std::swap
#include <array>
#include <iterator>   // std::back_inserter
#include <utility>    // std::pair, std::move
#include "fact_tree.h"

// =============================================================================
// R* consensus tree  (Jansson, Sung, Vu & Yiu 2016, Algorithmica 76:1224-1244,
// after Degnan et al. 2009).
//
// Definition.  Given k rooted trees on the same leaf set L (n = |L|):
//   * #ab|c = number of input trees in which the resolved triplet ab|c is
//     consistent (fan triplets count for NOTHING -- they have no impact);
//   * the majority resolved-triplet set is
//        R_maj = { ab|c : #ab|c > max(#ac|b, #bc|a) }   (strict PLURALITY);
//   * the R* tree is the unique tree tau with r(tau) subset of R_maj maximising
//     internal nodes.  Lemma 1.1: its clusters are EXACTLY the STRONG CLUSTERS
//     of R_maj, where A subset of L is strong iff  aa'|x in R_maj  for every pair
//     a,a' in A and every x not in A (the "for-all-outgroup" rule).
//
// CONSTRUCTION (Apresjan-style, as in Jansson et al. 2016): similarity ->
// laminar (single-linkage) candidate clusters -> filter to strong clusters ->
// build the tree.  The earlier implementation realised exactly this pipeline but
// (a) materialised a dense n^3 int triplet tensor -- an O(n^3) MEMORY wall that
// hard-capped n at 200 -- and (b) assembled in O(n^4).  This version keeps the
// (validated) pipeline and the similarity s(a,b) = #{x : ab|x in R_maj}, but:
//
//   STAGE 0  Per-tree O(1) LCA-depth.  For each tree build (Euler tour +
//            sparse-table RMQ, O(n log n)) and materialise its all-pairs
//            LCA-depth matrix D_t (O(n^2)).  Store the k matrices flat
//            (O(k n^2) ints) -- NO n^3 tensor.  Both the tally and the on-demand
//            R_maj test below read D_t in O(1).
//   STAGE 1  Tally.  Stream over the C(n,3) triples; per triple loop the k trees
//            (three O(1) LCA-depth reads each), take the strict-plurality winner,
//            and accumulate it into the n x n similarity matrix s.  O(k n^3) time,
//            O(n^2) memory for s.
//   STAGE 2  Candidates + filter + build.  Maximum spanning tree of s (Prim,
//            O(n^2)); its edges processed in decreasing weight give the laminar
//            single-linkage (Apresjan) dendrogram of candidate clusters -- a
//            provable superset of the strong clusters (strong subset Apresjan
//            subset single-linkage).  Filter bottom-up while building the R*
//            forest: a candidate's within-block pairs are already certified (a
//            strong child holds for every outside x), so only NEW cross-block
//            pairs are tested -- against every outside leaf via an on-demand
//            R_maj query (recomputed from D in O(k)) -- with the necessary
//            condition s(a,a') >= n-|A| as a cheap prefilter and early-exit.
//
// COMPLEXITY (honest).  Time is O(k n^3) (the tally; unchanged asymptotically --
// the sub-cubic bounds of Jansson et al. are "galactic", relying on
// Boolean-matrix-multiplication / dynamic-connectivity machinery that does MORE
// work than this for every feasible n).  Memory is O(k n^2) (the per-tree
// LCA-depth matrices) + O(n^2) (s) -- NO n^3 structure -- so the hard 200-leaf
// cap is gone; the practical limit is now runtime, not a memory wall.  Assembly
// drops from O(n^4) to ~O(n^3).  In practice: large constant-factor + memory
// wins, and the same R* tree (regression-checked exactly against the former
// build to n = 200 via dev/oracle/rstar/check-vs-legacy.R).
//
// Taxa are 0-indexed internally; tip ape-node v <-> taxon v-1 <-> fact node v-1.
// Output is integer-label Newick (1-indexed) with NO trailing ';'.
// =============================================================================

// ---------------------------------------------------------------------------
// closePair: from the three pairwise LCA depths of a triple {A,B,C} given as
// (dAB, dAC, dBC), return the close pair (the one whose LCA is strictly deepest):
//   0 -> AB|C ,  1 -> AC|B ,  2 -> BC|A ,  -1 -> fan / unresolved.
// (Logic identical to the validated former implementation: in any rooted tree
// two of the three LCAs coincide -- the triple's apex -- and the third is the
// strictly deeper close pair; all three equal => fan.)
// ---------------------------------------------------------------------------
static inline int closePair(int dAB, int dAC, int dBC) {
  if (dAB == dAC && dAB == dBC) return -1;   // fan -> abstain
  if (dAC == dBC && dAB > dAC)  return 0;    // AB|C
  if (dAB == dBC && dAC > dAB)  return 1;    // AC|B
  if (dAB == dAC && dBC > dAB)  return 2;    // BC|A
  return -1;                                 // defensive (unreachable on a valid tree)
}

// Disjoint-set find with path halving.
static int ufFind(std::vector<int>& f, int x) {
  while (f[x] != x) { f[x] = f[f[x]]; x = f[x]; }
  return x;
}

// ---------------------------------------------------------------------------
// Fill D[(a*n + b)*k + t] = depth(LCA(tip a, tip b)) in tree t, for all tips
// 0 <= a,b < n, via an Euler tour + sparse-table RMQ (O(n log n) build, O(1)
// query, O(n^2) fill).  The layout is PAIR-major / tree-minor: the k trees'
// depths for one leaf-pair are contiguous, so the O(kn^3) tally's inner loop
// over trees reads sequentially (cache-friendly).  Tip taxon a (0-indexed) is
// fact node a (buildTreeFromEdge maps ape node v -> fact node v-1); root
// depth = 0, child depth = parent depth + 1.
// ---------------------------------------------------------------------------
static void fillLcaDepths(const fact::Tree& tr, int n, int k, int t, int* D) {
  const int nNode = tr.cnt;
  std::vector<int> euler;     euler.reserve((size_t)2 * nNode);
  std::vector<int> firstPos(nNode, -1);
  std::vector<int> dep(nNode, 0);

  // Iterative Euler-tour DFS: record a node on first visit and again after each
  // returning child, so adjacent Euler entries are parent/child.
  {
    std::vector<std::pair<int, int> > st;
    st.reserve(nNode);
    dep[tr.root] = 0;
    st.push_back(std::make_pair(tr.root, 0));
    while (!st.empty()) {
      int node = st.back().first;
      int ci = st.back().second;
      if (ci == 0) { firstPos[node] = (int)euler.size(); euler.push_back(node); }
      if (ci < (int)tr.G[node].size()) {
        int c = tr.G[node][ci];
        st.back().second = ci + 1;
        dep[c] = dep[node] + 1;
        st.push_back(std::make_pair(c, 0));
      } else {
        st.pop_back();
        if (!st.empty()) euler.push_back(st.back().first);
      }
    }
  }

  const int m = (int)euler.size();
  std::vector<int> LOG(m + 1, 0);
  for (int i = 2; i <= m; ++i) LOG[i] = LOG[i / 2] + 1;
  const int K = LOG[m] + 1;
  // sp[j][i] = Euler index in [i, i + 2^j) whose node has minimum depth.
  std::vector<std::vector<int> > sp(K, std::vector<int>(m));
  for (int i = 0; i < m; ++i) sp[0][i] = i;
  for (int j = 1; j < K; ++j) {
    const int half = 1 << (j - 1);
    for (int i = 0; i + (1 << j) <= m; ++i) {
      int l = sp[j - 1][i], r = sp[j - 1][i + half];
      sp[j][i] = (dep[euler[l]] <= dep[euler[r]]) ? l : r;
    }
  }

  for (int a = 0; a < n; ++a) {
    const int fa = firstPos[a];
    D[((size_t)a * n + a) * k + t] = dep[a];
    for (int b = a + 1; b < n; ++b) {
      int l = fa, r = firstPos[b];
      if (l > r) std::swap(l, r);
      const int j = LOG[r - l + 1];
      const int i1 = sp[j][l], i2 = sp[j][r - (1 << j) + 1];
      const int lca = (dep[euler[i1]] <= dep[euler[i2]]) ? euler[i1] : euler[i2];
      const int d = dep[lca];
      D[((size_t)a * n + b) * k + t] = d;
      D[((size_t)b * n + a) * k + t] = d;
    }
  }
}

// --- direct ancestor-walk LCA depth (small-n self-check oracle only) ---------
static void buildParentDepth(const Rcpp::IntegerMatrix& edge, int nNode,
                             std::vector<int>& parent, std::vector<int>& depth) {
  parent.assign(nNode + 1, -1);
  depth.assign(nNode + 1, 0);
  for (int r = 0; r < edge.nrow(); ++r) {
    int p = edge(r, 0), c = edge(r, 1);
    parent[c] = p;
    depth[c] = depth[p] + 1;
  }
}
static int lcaDepthWalk(int u, int v, const std::vector<int>& parent,
                        const std::vector<int>& depth) {
  while (depth[u] > depth[v]) u = parent[u];
  while (depth[v] > depth[u]) v = parent[v];
  while (u != v) { u = parent[u]; v = parent[v]; }
  return depth[u];
}

// ---------------------------------------------------------------------------
// Main Rcpp entry point.  Returns the Newick string (no trailing ';').
// [[Rcpp::export]]
std::string rStarConsensus(Rcpp::List edgeList, int nTip) {
  const int n = nTip;
  const int k = edgeList.size();

  // Trivial leaf sets (the R wrapper handles n < 3, but stay self-contained).
  if (n < 3) {
    std::string out = "(";
    for (int i = 0; i < n; ++i) { if (i) out += ","; out += std::to_string(i + 1); }
    out += ")";
    return out;
  }

  // Memory guard on the k per-tree LCA-depth matrices (O(k n^2) ints).  This
  // replaces the former dense n^3 tensor; the ceiling is far higher (e.g.
  // n = 2000, k = 10 ~ 0.16 GB) and scales linearly in k, quadratically in n.
  const double bytes = (double)k * (double)n * (double)n * (double)sizeof(int);
  if (bytes > 2.4e9) {
    Rcpp::stop("rStarConsensus: the k per-tree LCA-depth matrices would need "
               "~%.1f GB (k = %d, n = %d). Reduce the number of trees or leaves.",
               bytes / 1e9, k, n);
  }

  // ---- Stage 0: per-tree O(1) LCA-depth matrices -----------------------------
  std::vector<int> D((size_t)k * n * n, 0);
  for (int t = 0; t < k; ++t) {
    Rcpp::checkUserInterrupt();
    Rcpp::IntegerMatrix edge = edgeList[t];
    fact::Tree tr = fact::buildTreeFromEdge(edge, n);
    fillLcaDepths(tr, n, k, t, D.data());
  }
  // O(1) LCA-depth accessor: depth of LCA(tip a, tip b) in tree t.  Pair-major
  // layout (see fillLcaDepths) -> the inner tree loop below reads contiguously.
  #define RSTAR_DEP(t, a, b) (D[(((size_t)(a) * n + (b)) * k + (t))])

  // Defensive: on small inputs, cross-check the O(1) LCA against a direct
  // ancestor walk (the new primitive's only failure mode).  Cheap; exercised by
  // every small test/oracle case.
  if (k >= 1 && n <= 50) {
    Rcpp::IntegerMatrix edge0 = edgeList[0];
    int nNode0 = 0;
    for (int r = 0; r < edge0.nrow(); ++r) {
      if (edge0(r, 0) > nNode0) nNode0 = edge0(r, 0);
      if (edge0(r, 1) > nNode0) nNode0 = edge0(r, 1);
    }
    std::vector<int> par, dpv;
    buildParentDepth(edge0, nNode0, par, dpv);
    for (int a = 0; a < n; ++a)
      for (int b = a + 1; b < n; ++b)
        if (lcaDepthWalk(a + 1, b + 1, par, dpv) != RSTAR_DEP(0, a, b))
          Rcpp::stop("rStarConsensus: internal LCA self-check failed "
                     "(a = %d, b = %d).", a, b);
  }

  // ---- Stage 1: tally -> similarity s(a,b) = #{x : ab|x in R_maj} ------------
  std::vector<int> s((size_t)n * n, 0);
  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    for (int j = i + 1; j < n; ++j) {
      for (int l = j + 1; l < n; ++l) {
        int cIJ = 0, cIL = 0, cJL = 0;
        for (int t = 0; t < k; ++t) {
          int cp = closePair(RSTAR_DEP(t, i, j), RSTAR_DEP(t, i, l), RSTAR_DEP(t, j, l));
          if      (cp == 0) ++cIJ;   // ij|l
          else if (cp == 1) ++cIL;   // il|j
          else if (cp == 2) ++cJL;   // jl|i
        }
        int mx = cIJ; if (cIL > mx) mx = cIL; if (cJL > mx) mx = cJL;
        if (mx == 0) continue;                                       // all fans
        if ((cIJ == mx) + (cIL == mx) + (cJL == mx) != 1) continue;  // tie
        if      (cIJ == mx) s[(size_t)i * n + j] += 1;               // ij|l in R_maj
        else if (cIL == mx) s[(size_t)i * n + l] += 1;               // il|j in R_maj
        else                s[(size_t)j * n + l] += 1;               // jl|i in R_maj
      }
    }
  }
  for (int a = 0; a < n; ++a)
    for (int b = a + 1; b < n; ++b)
      s[(size_t)b * n + a] = s[(size_t)a * n + b];

  // On-demand R_maj membership: is ab|x in R_maj?  (Recomputed from D in O(k);
  // strict plurality: ab|x wins iff #ab|x > #ax|b and #ab|x > #bx|a.)
  auto inRmaj = [&](int a, int b, int x) -> bool {
    int cab = 0, cax = 0, cbx = 0;
    for (int t = 0; t < k; ++t) {
      int cp = closePair(RSTAR_DEP(t, a, b), RSTAR_DEP(t, a, x), RSTAR_DEP(t, b, x));
      if      (cp == 0) ++cab;   // ab|x
      else if (cp == 1) ++cax;   // ax|b
      else if (cp == 2) ++cbx;   // bx|a
    }
    return cab > 0 && cab > cax && cab > cbx;
  };

  // ---- Stage 2a: maximum spanning tree of s (Prim, O(n^2)) -------------------
  std::vector<int> bestW(n, -1), bestTo(n, -1);
  std::vector<char> inTree(n, 0);
  std::vector<std::array<int, 3> > mst;        // (weight, u, v)
  mst.reserve(n - 1);
  inTree[0] = 1;
  for (int v = 1; v < n; ++v) { bestW[v] = s[(size_t)0 * n + v]; bestTo[v] = 0; }
  for (int iter = 1; iter < n; ++iter) {
    int u = -1, w = -1;
    for (int v = 0; v < n; ++v)
      if (!inTree[v] && bestW[v] > w) { w = bestW[v]; u = v; }
    inTree[u] = 1;
    std::array<int, 3> e = { bestW[u], bestTo[u], u };
    mst.push_back(e);
    for (int v = 0; v < n; ++v)
      if (!inTree[v]) {
        int sw = s[(size_t)u * n + v];
        if (sw > bestW[v]) { bestW[v] = sw; bestTo[v] = u; }
      }
  }
  std::sort(mst.begin(), mst.end(),
            [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
              return a[0] > b[0];   // decreasing weight -> single-linkage merges
            });

  // ---- Stage 2b: dendrogram candidates -> filter to strong -> build ----------
  std::vector<int> dpar(n);   std::iota(dpar.begin(), dpar.end(), 0);   // dendrogram UF
  std::vector<std::vector<int> > members(n);
  for (int i = 0; i < n; ++i) members[i].assign(1, i);
  std::vector<int> apar(n);   std::iota(apar.begin(), apar.end(), 0);   // accepted-block UF
  std::vector<int> repNode(n);                                          // block root -> out node
  std::iota(repNode.begin(), repNode.end(), 0);                         // leaf i -> out leaf i

  std::vector<std::vector<int> > outChildren((size_t)2 * n);            // out node -> children
  int nextId = n;                                                       // internal ids >= n
  std::vector<char> inA(n, 0);

  for (size_t e = 0; e < mst.size(); ++e) {
    Rcpp::checkUserInterrupt();
    int ru = ufFind(dpar, mst[e][1]);
    int rv = ufFind(dpar, mst[e][2]);
    if (ru == rv) continue;                          // (cannot happen: MST is acyclic)

    // candidate A = members[ru] U members[rv] (both sorted ascending)
    std::vector<int> A;
    A.reserve(members[ru].size() + members[rv].size());
    std::merge(members[ru].begin(), members[ru].end(),
               members[rv].begin(), members[rv].end(), std::back_inserter(A));
    const int sz = (int)A.size();

    if (sz >= 2 && sz <= n - 1) {
      for (int x : A) inA[x] = 1;

      // partition A into current accepted blocks (distinct accUF roots)
      std::vector<int> blockRoot;
      std::vector<std::vector<int> > blockMem;
      for (int a : A) {
        int r = ufFind(apar, a), bi = -1;
        for (size_t z = 0; z < blockRoot.size(); ++z)
          if (blockRoot[z] == r) { bi = (int)z; break; }
        if (bi < 0) { bi = (int)blockRoot.size(); blockRoot.push_back(r); blockMem.push_back(std::vector<int>()); }
        blockMem[bi].push_back(a);
      }
      const int p = (int)blockRoot.size();

      // strong iff every cross-block pair (a,a') has aa'|x in R_maj for all x not
      // in A.  Verify by counting the SMALLER side: if A is small, count the
      // inside x with aa'|x in R_maj and derive the outside count via s (the
      // pair's total over all x); else scan the outside directly with early-exit.
      bool strong = (p >= 2);
      const int need = n - sz;                       // = |L \ A| (#outside)
      const bool scanOutside = (need <= sz);
      for (int bi = 0; bi < p && strong; ++bi)
        for (int bj = bi + 1; bj < p && strong; ++bj)
          for (size_t ia = 0; ia < blockMem[bi].size() && strong; ++ia)
            for (size_t ja = 0; ja < blockMem[bj].size() && strong; ++ja) {
              int a = blockMem[bi][ia], ap = blockMem[bj][ja];
              const int sab = s[(size_t)a * n + ap];
              if (sab < need) { strong = false; break; }   // necessary condition
              if (scanOutside) {
                for (int x = 0; x < n; ++x)
                  if (!inA[x] && !inRmaj(a, ap, x)) { strong = false; break; }
              } else {
                int inside = 0;                            // x in A with aa'|x in R_maj
                for (int x2 : A)
                  if (x2 != a && x2 != ap && inRmaj(a, ap, x2)) ++inside;
                if (sab - inside != need) strong = false;  // some outside x failed
              }
            }

      if (strong) {                                  // accept: create an R* node
        int id = nextId++;
        for (int bi = 0; bi < p; ++bi) outChildren[id].push_back(repNode[blockRoot[bi]]);
        int base = blockRoot[0];
        for (int bi = 1; bi < p; ++bi) apar[ufFind(apar, blockRoot[bi])] = ufFind(apar, base);
        repNode[ufFind(apar, base)] = id;
      }

      for (int x : A) inA[x] = 0;
    }

    // unite dendrogram UF (attach the smaller member list under the larger)
    if (members[ru].size() >= members[rv].size()) {
      dpar[rv] = ru; members[ru] = std::move(A); members[rv].clear();
    } else {
      dpar[ru] = rv; members[rv] = std::move(A); members[ru].clear();
    }
  }

  // ---- root: the full leaf set is the root; join all top-level blocks --------
  const int rootId = nextId++;
  {
    std::vector<char> seen(n, 0);
    for (int leaf = 0; leaf < n; ++leaf) {
      int r = ufFind(apar, leaf);
      if (!seen[r]) { seen[r] = 1; outChildren[rootId].push_back(repNode[r]); }
    }
  }

  // ---- emit Newick (iterative post-order) ------------------------------------
  std::vector<std::string> memo(nextId);
  std::vector<int> order; order.reserve(nextId);
  {
    std::vector<int> stk(1, rootId);
    while (!stk.empty()) {
      int u = stk.back(); stk.pop_back();
      order.push_back(u);
      for (int c : outChildren[u]) stk.push_back(c);
    }
  }
  for (int idx = (int)order.size() - 1; idx >= 0; --idx) {
    int u = order[idx];
    if (outChildren[u].empty()) {                    // leaf: 1-indexed integer label
      memo[u] = std::to_string(u + 1);
    } else {
      std::string out = "(";
      for (size_t c = 0; c < outChildren[u].size(); ++c) {
        if (c) out += ",";
        out += memo[outChildren[u][c]];
      }
      out += ")";
      memo[u] = out;
    }
  }

  #undef RSTAR_DEP
  return memo[rootId];
}
