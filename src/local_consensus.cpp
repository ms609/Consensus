#include <Rcpp.h>
#include <vector>
#include <climits>
#include <bitset>
#include <string>
#include <cstdlib>  // abs

// =============================================================================
// Rcpp port of local_consensus.h from FDCT_new / FACT2 by Jansson et al.
// MinRLC (minrs=true)  = minRLC_exact (rooted local consensus)
// MinILC (minrs=false) = minILC_exact (induced local consensus)
//
// Taxa are 0-indexed bit positions throughout C++.  Edge matrix from R is
// 1-indexed; tip node v (1-indexed) corresponds to bit v-1 (0-indexed).
// =============================================================================

typedef unsigned int uint;

// comb2(x) = x*(x-1)/2  (used in the MinILC cost term)
static inline int comb2(int x) { return x * (x - 1) / 2; }

// ---------------------------------------------------------------------------
// Build parent[] and depth[] for a single tree from its ape edge matrix.
// parent[v] = parent of node v (1-indexed); -1 for root.
// depth[v]  = depth of node v (root = 0).
// Since ape::Preorder() is guaranteed, a single forward pass suffices
// (parent always appears before child in edge list).
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

// ---------------------------------------------------------------------------
// Depth of LCA of two 1-indexed nodes (leaves are 1..nTip).
// ---------------------------------------------------------------------------
static int lcaDepth(int u, int v,
                    const std::vector<int>& parent,
                    const std::vector<int>& depth) {
  while (depth[u] > depth[v]) u = parent[u];
  while (depth[v] > depth[u]) v = parent[v];
  while (u != v) { u = parent[u]; v = parent[v]; }
  return depth[u];
}

// ---------------------------------------------------------------------------
// Build triplet counts across all trees.
// triplets[i][j][k] (0-indexed) = number of trees in which ij|k is resolved.
// After all trees, values < nTree are zeroed (keep only universal triplets).
// ---------------------------------------------------------------------------
static std::vector<std::vector<std::vector<int>>>
getCommonTriplets(const Rcpp::List& edgeList, int nTip) {
  int nTree = edgeList.size();
  // Zero-initialise
  std::vector<std::vector<std::vector<int>>> tri(
      nTip, std::vector<std::vector<int>>(nTip, std::vector<int>(nTip, 0)));

  for (int t = 0; t < nTree; t++) {
    Rcpp::IntegerMatrix edge = edgeList[t];
    // Determine max node index
    int nNode = 0;
    for (int r = 0; r < edge.nrow(); r++) {
      if (edge(r, 0) > nNode) nNode = edge(r, 0);
      if (edge(r, 1) > nNode) nNode = edge(r, 1);
    }
    std::vector<int> par, dep;
    buildParentDepth(edge, nNode, par, dep);

    // Taxon i (0-indexed) = leaf node i+1 (1-indexed)
    for (int i = 0; i < nTip; i++) {
      for (int j = i + 1; j < nTip; j++) {
        int dij = lcaDepth(i + 1, j + 1, par, dep);
        for (int k = j + 1; k < nTip; k++) {
          int dik = lcaDepth(i + 1, k + 1, par, dep);
          int djk = lcaDepth(j + 1, k + 1, par, dep);

          // Fan: all three LCA depths equal → no resolved triplet
          if (dij == dik && dij == djk) continue;

          if (dik == djk && dij > dik) {        // ij|k
            tri[i][j][k]++;
          } else if (dij == djk && dik > dij) { // ik|j
            tri[i][k][j]++;
          } else if (dij == dik && djk > dij) { // jk|i
            tri[j][k][i]++;
          }
        }
      }
    }
  }

  // Zero out any triplet not universal
  for (int i = 0; i < nTip; i++)
    for (int j = 0; j < nTip; j++)
      for (int k = 0; k < nTip; k++)
        if (tri[i][j][k] < nTree) tri[i][j][k] = 0;

  return tri;
}

// ---------------------------------------------------------------------------
// DFS for Aho component (recursive, matches reference exactly).
// indices[v] = global (0-indexed) taxon id of local vertex v.
// ---------------------------------------------------------------------------
static uint dfs(const std::vector<std::vector<int>>& adjl,
                int v, std::vector<bool>& visited,
                const std::vector<int>& indices) {
  uint comp = (1u << indices[v]);
  visited[v] = true;
  for (int nb : adjl[v]) {
    if (!visited[nb]) comp |= dfs(adjl, nb, visited, indices);
  }
  return comp;
}

// ---------------------------------------------------------------------------
// Build Aho graph for a leaf subset and return connected-component bitmasks.
// Singletons (popcount == 1) are dropped when minrs == true.
// ---------------------------------------------------------------------------
static std::vector<uint>
buildAho(const std::vector<int>& indices,
         const std::vector<std::vector<std::vector<int>>>& tri,
         bool minrs) {
  int n = (int)indices.size();
  std::vector<std::vector<int>> adjl(n);

  for (int i = 0; i < n; i++) {
    int ii = indices[i];
    for (int j = i + 1; j < n; j++) {
      int ij = indices[j];
      // Edge i–j iff some k in subset has tri[ii][ij][ik] != 0
      for (int k = 0; k < n; k++) {
        int ik = indices[k];
        if (tri[ii][ij][ik]) {
          adjl[i].push_back(j);
          adjl[j].push_back(i);
          break;
        }
      }
    }
  }

  std::vector<bool> visited(n, false);
  std::vector<uint> components;
  for (int i = 0; i < n; i++) {
    if (!visited[i]) {
      uint comp = dfs(adjl, i, visited, indices);
      if ((uint)std::bitset<20>(comp).count() > (uint)minrs) {
        components.push_back(comp);
      }
    }
  }
  return components;
}

// ---------------------------------------------------------------------------
// Newick-string construction — replaces print_tree's Tree-building.
// Recursively builds Newick for leaf-set bitmask Lbitmask.
// Returns a string that, when wrapped in "()", gives the sub-tree.
// If there are multiple children, the returned string already has outer parens.
// ---------------------------------------------------------------------------
static std::string buildNewick(
    uint Lbitmask,
    const std::vector<std::vector<uint>>& comps,
    const std::vector<std::vector<int>>&  dpBT,
    bool minrs, int n) {

  std::vector<std::string> children;
  int m = (int)comps[Lbitmask].size();

  // Rooted: singletons in Lbitmask not covered by any component
  if (minrs) {
    uint covered = 0;
    for (int i = 0; i < m; i++) covered |= comps[Lbitmask][i];
    uint singletons = Lbitmask & ~covered;
    for (int bit = 0; bit < n; bit++) {
      if (singletons & (1u << bit)) children.push_back(std::to_string(bit + 1));
    }
  }

  if (m == 0) {
    // nothing more — handled above (singletons only)
  } else if (m == 1) {
    // Recurse into single component
    children.push_back(buildNewick(comps[Lbitmask][0], comps, dpBT, minrs, n));
  } else {
    // Walk dp_backtrack to peel off groups (exactly as print_tree)
    uint Dbitmask = (1u << m) - 1u;
    int  Xbitmask  = dpBT[Lbitmask][(1 << m) - 1];
    int  origXbitmask;

    do {
      origXbitmask = Xbitmask;
      Xbitmask = std::abs(Xbitmask);

      uint DmXbitmask = Dbitmask & ~(uint)Xbitmask;
      uint lambdaX = 0;
      for (int i = 0; i < m; i++) {
        if ((uint)Xbitmask & (1u << i)) lambdaX |= comps[Lbitmask][i];
      }

      if (!minrs && (int)std::bitset<20>(lambdaX).count() == 1) {
        for (int bit = 0; bit < n; bit++) {
          if (lambdaX & (1u << bit)) children.push_back(std::to_string(bit + 1));
        }
      } else {
        children.push_back(buildNewick(lambdaX, comps, dpBT, minrs, n));
      }

      Dbitmask = DmXbitmask;
      Xbitmask = dpBT[Lbitmask][DmXbitmask];
    } while (origXbitmask >= 0);

    // Tail: remaining bits
    Xbitmask = std::abs(Xbitmask);  // ensure positive (not used further here)
    uint lambdaX = 0;
    for (int i = 0; i < m; i++) {
      if (Dbitmask & (1u << i)) lambdaX |= comps[Lbitmask][i];
    }

    if (!minrs && (int)std::bitset<20>(lambdaX).count() == 1) {
      for (int bit = 0; bit < n; bit++) {
        if (lambdaX & (1u << bit)) children.push_back(std::to_string(bit + 1));
      }
    } else {
      children.push_back(buildNewick(lambdaX, comps, dpBT, minrs, n));
    }
  }

  // Collapse: single child → pass through; multiple → wrap in parens.
  if (children.empty()) return "";
  if ((int)children.size() == 1) return children[0];
  std::string result = "(";
  for (int i = 0; i < (int)children.size(); i++) {
    if (i > 0) result += ",";
    result += children[i];
  }
  result += ")";
  return result;
}


// ---------------------------------------------------------------------------
// Gosper's hack: next integer with the same number of set bits.
// ---------------------------------------------------------------------------
static inline uint nextBitPerm(uint v) {
  uint t = (v | (v - 1)) + 1;
  return t | (((t & (uint)-(int)t) / (v & (uint)-(int)v) >> 1) - 1);
}


// ---------------------------------------------------------------------------
// Main Rcpp entry point.
// [[Rcpp::export]]
std::string localConsensus(Rcpp::List edgeList, int nTip, bool minrs) {
  int n = nTip;
  if (n > 20) Rcpp::stop("localConsensus: n > 20 not supported (exact-exponential)");

  auto tri = getCommonTriplets(edgeList, n);

  uint states = (1u << n);

  std::vector<int>              opt(states, 0);
  std::vector<int>              dp(states, 0);
  std::vector<std::vector<uint>> comps(states);
  std::vector<std::vector<int>>  dpBT(states);  // dpBT[L] sized 2^m after allocation

  for (int i = 2; i <= n; i++) {
    uint Lbitmask = (1u << i) - 1u;  // lowest i-bit mask

    while ((Lbitmask & (1u << n)) == 0) {
      // Build indices
      std::vector<int> indices;
      for (int j = 0; j < n; j++) {
        if (Lbitmask & (1u << j)) indices.push_back(j);
      }

      comps[Lbitmask] = buildAho(indices, tri, minrs);
      int m = (int)comps[Lbitmask].size();

      // NULL sentinel: inseparable full component
      if (m == 1 &&
          (int)std::bitset<20>(comps[Lbitmask][0]).count() ==
          (int)std::bitset<20>(Lbitmask).count()) {
        return "";
      }

      dpBT[Lbitmask].assign(1 << m, 0);

      // Init DP for single components (j = 0..m-1 → bitmask 1<<j)
      for (int j = 0; j < m; j++) {
        dp[1 << j] = opt[comps[Lbitmask][j]];
        if (!minrs) {
          int cc = (int)std::bitset<20>(comps[Lbitmask][j]).count();
          dp[1 << j] += comb2(cc) * (i - cc);
        }
      }

      // DP over D ⊆ {0..m-1}, |D| >= 2
      for (int j = 2; j <= m; j++) {
        uint Dbitmask = (1u << j) - 1u;
        while ((Dbitmask & (1u << m)) == 0) {
          // The inner DP is O(3^m); on highly incongruent input m can approach
          // n, so let the user abort a runaway with Ctrl-C.
          Rcpp::checkUserInterrupt();
          // Collect indices of set bits in Dbitmask
          std::vector<int> idx2;
          for (int k = 0; k < m; k++) {
            if (Dbitmask & (1u << k)) idx2.push_back(k);
          }
          dp[Dbitmask] = INT_MAX;

          int q = (int)idx2.size();
          // Iterate over proper non-empty subsets X of D
          // (subX is a bitmask over idx2, i.e., subX ∈ [1, 2^q - 2])
          uint upperBm = (1u << q) - 1u;
          uint subX = 1u;
          while (subX < upperBm) {
            uint lambdaX = 0, lambdaDmX = 0;
            uint Xbitmask = 0, DmXbitmask = 0;
            for (int k = 0; k < q; k++) {
              if (subX & (1u << k)) {
                lambdaX  |= comps[Lbitmask][idx2[k]];
                Xbitmask |= (1u << idx2[k]);
              } else {
                lambdaDmX  |= comps[Lbitmask][idx2[k]];
                DmXbitmask |= (1u << idx2[k]);
              }
            }

            int optX   = opt[lambdaX];
            int optDmX = opt[lambdaDmX];
            if (!minrs) {
              int cx   = (int)std::bitset<20>(lambdaX).count();
              int cdmx = (int)std::bitset<20>(lambdaDmX).count();
              optX   += comb2(cx)   * (i - cx);
              optDmX += comb2(cdmx) * (i - cdmx);
            }

            int min2 = std::min(dp[DmXbitmask], optDmX);
            if (dp[Dbitmask] > optX + min2) {
              dp[Dbitmask] = optX + min2;
              dpBT[Lbitmask][Dbitmask] = (min2 == optDmX)
                                           ? -(int)Xbitmask
                                           :  (int)Xbitmask;
            }

            subX = nextBitPerm(subX);
          }

          Dbitmask = nextBitPerm(Dbitmask);
        }
      }

      opt[Lbitmask] = dp[(1 << m) - 1] + (int)minrs;

      Lbitmask = nextBitPerm(Lbitmask);
    }
  }

  // Build Newick for the full leaf set
  uint fullMask = states - 1u;
  std::string nwk = buildNewick(fullMask, comps, dpBT, minrs, n);
  return nwk;
}


// Scaffold self-check (kept for backward compatibility)
// [[Rcpp::export]]
int consensus_rcpp_selfcheck() {
  return 42;
}
