#include <Rcpp.h>
#include "fact_tree.h"

#include <algorithm>
#include <climits>
#include <string>
#include <utility>
#include <vector>

// Adams consensus -- the O(kn log n) algorithm of Jansson, Li & Sung (2017), "On
// finding the Adams consensus tree", Information and Computation 256:334-347
// (New_Adams_consensus_k, Section 3), implemented from the paper.  There is no
// portable reference C++ (FACT's adams.cpp is the classical O(kn^2) recursion and
// its fast binary mis-prints), so this is validated clade-for-clade against the
// classical slow Adams (the definition) -- dev/oracle: FACT rule 512, rooted=1.
//
// Definition (the clade contract): given k ROOTED trees on the same n leaves, at
// a node responsible for leaf-set L' the children are the blocks of the MEET of
// the per-tree partitions "which child of lca_{T_j}(L') does each leaf descend
// from"; recurse per block until |L'| <= 2.  The Adams tree is UNIQUE, so the
// centroid-path device below does not affect the output.  Adams is ROOTED: the R
// wrapper passes each tree on its OWN root.
//
// Fast version: rather than recursing on every block, the block whose leaves lie
// under the heavy child (most leaf descendants) of the root in EVERY tree -- the
// "spine" -- is expanded ITERATIVELY down the centroid paths of all k trees in
// unison; only the off-spine "side" blocks are recursed, and each is <= |L'|/2 (a
// side tree of a centroid path holds < half the leaves), giving O(log n) depth.
//
// Per tree we precompute, for every leaf, the spine level at which it leaves the
// path (branchLevel) and the child it leaves through (sideChild).  The spine walk
// (adamsRecurse) maintains the remaining block B and, per tree, a pointer into
// the leaves sorted by branchLevel: currentPos_j = min branchLevel over B = the
// depth of lca_{T_j}(B) (this is where path compression is handled -- when leaves
// have already split off, a tree's position jumps to the true LCA).  Each step
// splits off every leaf that branches at some tree's currentPos, partitions them
// by the meet key, recurses those side blocks, and continues with the rest.  The
// emitted spine is the nested chain of these step nodes.
//
// Reentrant: no file-scope mutable state; all context is passed by argument.
// VLA-free (std::vector throughout); RAII.

namespace {

using std::vector;

// ---------------------------------------------------------------------------
// InTree: one input tree, built once, immutable through the recursion.  Carries
// an Euler tour + sparse-table RMQ for O(1) lowest-common-ancestor queries, used
// to build the restricted subtree (nub) at each recursion level.
// ---------------------------------------------------------------------------
struct InTree {
  int nNode = 0;
  int root = 0;
  int nTip = 0;
  vector<int> leafTaxon;     // leafTaxon[node] = taxon (1..nTip) or 0 (internal)
  vector<int> leafNode;      // leafNode[taxon] = node id (size nTip + 1)
  vector<int> euler;         // node id at each Euler step (length 2*nNode - 1)
  vector<int> eulerDepth;    // depth at each Euler step
  vector<int> firstVisit;    // firstVisit[node] = first index in euler
  vector<int> lastVisit;     // lastVisit[node]  = last index in euler
  vector<int> logTab;        // logTab[len] = floor(log2(len))
  vector<int> sparse;        // RMQ table, [level * eulerLen + i] -> euler index
  int eulerLen = 0;
  int nLevels = 0;
  vector<int> taxaEulerOrder;  // every taxon, in firstVisit order

  int lca(int u, int v) const {
    int a = firstVisit[u], b = firstVisit[v];
    if (a > b) std::swap(a, b);
    int k = logTab[b - a + 1];
    int i1 = sparse[static_cast<size_t>(k) * eulerLen + a];
    int i2 = sparse[static_cast<size_t>(k) * eulerLen + (b - (1 << k) + 1)];
    return eulerDepth[i1] <= eulerDepth[i2] ? euler[i1] : euler[i2];
  }
  bool isAncestor(int a, int u) const {
    return firstVisit[a] <= firstVisit[u] && firstVisit[u] <= lastVisit[a];
  }
};

InTree buildInTree(const fact::Tree& t) {
  InTree g;
  g.nNode = t.cnt;
  g.root = t.root;
  g.nTip = t.N;
  g.leafTaxon.assign(t.cnt, 0);
  g.leafNode.assign(t.N + 1, 0);
  for (int v = 0; v < t.cnt; ++v) {
    g.leafTaxon[v] = t.leaf[v];
    if (t.leaf[v] > 0) g.leafNode[t.leaf[v]] = v;
  }

  // Euler tour by iterative DFS over the children lists t.G.
  g.firstVisit.assign(t.cnt, -1);
  g.lastVisit.assign(t.cnt, -1);
  g.euler.reserve(static_cast<size_t>(2) * t.cnt);
  g.eulerDepth.reserve(static_cast<size_t>(2) * t.cnt);
  vector<int> dep(t.cnt, 0);
  vector<int> stNode, stIdx;
  stNode.reserve(t.cnt);
  stIdx.reserve(t.cnt);
  stNode.push_back(t.root);
  stIdx.push_back(0);
  while (!stNode.empty()) {
    int v = stNode.back();
    int ci = stIdx.back();
    if (ci == 0) {
      g.firstVisit[v] = static_cast<int>(g.euler.size());
      g.euler.push_back(v);
      g.eulerDepth.push_back(dep[v]);
    }
    if (ci < static_cast<int>(t.G[v].size())) {
      int c = t.G[v][ci];
      stIdx.back() = ci + 1;
      dep[c] = dep[v] + 1;
      stNode.push_back(c);
      stIdx.push_back(0);
    } else {
      g.lastVisit[v] = static_cast<int>(g.euler.size()) - 1;
      stNode.pop_back();
      stIdx.pop_back();
      if (!stNode.empty()) {
        int p = stNode.back();
        g.euler.push_back(p);
        g.eulerDepth.push_back(dep[p]);
      }
    }
  }

  // Sparse table over eulerDepth (returns the euler index of the minimum).
  g.eulerLen = static_cast<int>(g.euler.size());
  g.logTab.assign(g.eulerLen + 1, 0);
  for (int i = 2; i <= g.eulerLen; ++i) g.logTab[i] = g.logTab[i / 2] + 1;
  g.nLevels = g.logTab[g.eulerLen] + 1;
  g.sparse.assign(static_cast<size_t>(g.nLevels) * g.eulerLen, 0);
  for (int i = 0; i < g.eulerLen; ++i) g.sparse[i] = i;
  for (int k = 1; k < g.nLevels; ++k) {
    size_t base = static_cast<size_t>(k) * g.eulerLen;
    size_t prev = static_cast<size_t>(k - 1) * g.eulerLen;
    int half = 1 << (k - 1);
    for (int i = 0; i + (1 << k) <= g.eulerLen; ++i) {
      int a = g.sparse[prev + i];
      int b = g.sparse[prev + i + half];
      g.sparse[base + i] = g.eulerDepth[a] <= g.eulerDepth[b] ? a : b;
    }
  }

  // Taxa in firstVisit order (the Euler tour visits leaves in that order).
  g.taxaEulerOrder.reserve(g.nTip);
  for (int i = 0; i < g.eulerLen; ++i) {
    int v = g.euler[i];
    if (g.leafTaxon[v] > 0 && g.firstVisit[v] == i)
      g.taxaEulerOrder.push_back(g.leafTaxon[v]);
  }
  return g;
}

// ---------------------------------------------------------------------------
// Nub: the restricted tree T_j|L', built fresh per recursion level via the
// consecutive-LCA auxiliary-tree construction.  Node ids are local (0-indexed in
// firstVisit/preorder); root = node 0 = lca_{T_j}(L'); leafTaxon[id] = taxon or
// 0; ch = children lists.
// ---------------------------------------------------------------------------
struct Nub {
  int cnt = 0;
  int root = 0;
  vector<vector<int> > ch;
  vector<int> leafTaxon;
  vector<int> cntLeaves;
};

// `leaves` are taxa in g's firstVisit order.
Nub buildNub(const InTree& g, const vector<int>& leaves) {
  int m = static_cast<int>(leaves.size());
  Nub nub;
  if (m == 1) {
    nub.cnt = 1;
    nub.root = 0;
    nub.ch.assign(1, vector<int>());
    nub.leafTaxon.assign(1, leaves[0]);
    nub.cntLeaves.assign(1, 1);
    return nub;
  }
  // Auxiliary-tree node set: the leaf nodes plus the LCA of each consecutive pair.
  vector<int> nodes;
  nodes.reserve(static_cast<size_t>(2) * m);
  for (int t : leaves) nodes.push_back(g.leafNode[t]);
  for (int i = 1; i < m; ++i)
    nodes.push_back(g.lca(g.leafNode[leaves[i - 1]], g.leafNode[leaves[i]]));
  // Order by firstVisit (preorder) with a stable LSD radix on the bounded key
  // firstVisit < eulerLen.  O(m), not the O(m log m) of a comparison sort, so the
  // recursion as a whole stays within the paper's O(kn log n) bound.  Duplicate
  // node ids share a firstVisit, so they land adjacent and std::unique drops them.
  {
    int ns = static_cast<int>(nodes.size());
    vector<int> tmp(ns);
    for (int shift = 0; (g.eulerLen >> shift) > 0; shift += 8) {
      int cnt[257];
      for (int b = 0; b < 257; ++b) cnt[b] = 0;
      for (int x : nodes) ++cnt[((g.firstVisit[x] >> shift) & 0xFF) + 1];
      for (int b = 0; b < 256; ++b) cnt[b + 1] += cnt[b];
      for (int x : nodes) tmp[cnt[(g.firstVisit[x] >> shift) & 0xFF]++] = x;
      nodes.swap(tmp);
    }
  }
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

  int cnt = static_cast<int>(nodes.size());
  nub.cnt = cnt;
  nub.ch.assign(cnt, vector<int>());
  nub.leafTaxon.assign(cnt, 0);
  nub.cntLeaves.assign(cnt, 0);
  // Stack build over nodes in firstVisit (preorder) order; nub id = index.
  vector<int> stk;
  stk.reserve(cnt);
  for (int i = 0; i < cnt; ++i) {
    while (!stk.empty() && !g.isAncestor(nodes[stk.back()], nodes[i]))
      stk.pop_back();
    if (!stk.empty()) nub.ch[stk.back()].push_back(i);
    stk.push_back(i);
  }
  nub.root = 0;  // smallest firstVisit = root of the induced subtree
  for (int i = 0; i < cnt; ++i) {
    int t = g.leafTaxon[nodes[i]];
    if (t > 0) nub.leafTaxon[i] = t;
  }
  // Leaf counts: children have larger index than their parent (preorder).
  for (int i = cnt - 1; i >= 0; --i) {
    if (nub.leafTaxon[i] > 0) {
      nub.cntLeaves[i] = 1;
    } else {
      int s = 0;
      for (int c : nub.ch[i]) s += nub.cntLeaves[c];
      nub.cntLeaves[i] = s;
    }
  }
  return nub;
}

// Compute, from a nub, the centroid path and per-leaf branch data.
//   branchLevel[taxon] = spine level (1-based) where the leaf leaves the path,
//   sideChild[taxon]   = the nub child it leaves through (a nub node id),
//   heavyChild[w]      = the heavy child of the level-w spine node (nub node id),
//                        for w = 1 .. spineLen - 1.
// branchLevel/sideChild are caller scratch indexed by taxon (only active taxa are
// written); heavyChild is sized to the spine length here.
void processNub(const Nub& nub, vector<int>& branchLevel, vector<int>& sideChild,
                vector<int>& heavyChild) {
  vector<int> spine;
  spine.push_back(nub.root);
  int cur = nub.root;
  while (nub.leafTaxon[cur] == 0) {  // descend the heavy child
    int best = -1, bc = -1;
    for (int c : nub.ch[cur])
      if (nub.cntLeaves[c] > bc) { bc = nub.cntLeaves[c]; best = c; }
    spine.push_back(best);
    cur = best;
  }
  int len = static_cast<int>(spine.size());
  heavyChild.assign(len + 2, -1);
  for (int w = 1; w <= len - 1; ++w) heavyChild[w] = spine[w];

  vector<int> sweep;
  for (int w = 1; w <= len; ++w) {
    int sp = spine[w - 1];
    if (nub.leafTaxon[sp] > 0) {  // the spine bottom leaf
      branchLevel[nub.leafTaxon[sp]] = w;
      sideChild[nub.leafTaxon[sp]] = sp;
      continue;
    }
    int heavy = (w <= len - 1) ? spine[w] : -1;
    for (int c : nub.ch[sp]) {
      if (c == heavy) continue;
      sweep.clear();
      sweep.push_back(c);
      while (!sweep.empty()) {
        int v = sweep.back();
        sweep.pop_back();
        if (nub.leafTaxon[v] > 0) {
          branchLevel[nub.leafTaxon[v]] = w;
          sideChild[nub.leafTaxon[v]] = c;
        } else {
          for (int cc : nub.ch[v]) sweep.push_back(cc);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// OutTree: the consensus tree being assembled, as a node array.
// ---------------------------------------------------------------------------
struct OutTree {
  vector<int> leafTaxon;        // 0 for internal nodes
  vector<vector<int> > ch;
  int addLeaf(int taxon) {
    int id = static_cast<int>(leafTaxon.size());
    leafTaxon.push_back(taxon);
    ch.push_back(vector<int>());
    return id;
  }
  int addInternal() {
    int id = static_cast<int>(leafTaxon.size());
    leafTaxon.push_back(0);
    ch.push_back(vector<int>());
    return id;
  }
  void link(int parent, int child) { ch[parent].push_back(child); }
};

// Per-recursion scratch indexed by taxon, reused down the recursion: each call
// fully consumes its scratch for its own (active) leaves before it recurses, so
// child calls may overwrite them.
struct Scratch {
  int k;
  vector<vector<int> > branchLevel;  // [k][nTip + 1]
  vector<vector<int> > sideChild;    // [k][nTip + 1]
  vector<int> blockOf;               // [nTip + 1], -1 = not in a side block
  vector<char> removed;              // [nTip + 1]
  vector<int> stamp;                 // [nTip + 1], per-step dedup
  int stampCtr;
  vector<int> compStamp;             // [2*nTip + 5], coord-compression epoch
  vector<int> compId;                // [2*nTip + 5], compressed dense id
  int compCtr;
  Scratch(int k_, int nTip)
      : k(k_),
        branchLevel(k_, vector<int>(nTip + 1, 0)),
        sideChild(k_, vector<int>(nTip + 1, 0)),
        blockOf(nTip + 1, -1),
        removed(nTip + 1, 0),
        stamp(nTip + 1, 0),
        stampCtr(0),
        compStamp(static_cast<size_t>(2) * nTip + 5, 0),
        compId(static_cast<size_t>(2) * nTip + 5, 0),
        compCtr(0) {}
};

int makeCherry(OutTree& out, int a, int b) {
  int p = out.addInternal();
  out.link(p, out.addLeaf(a));
  out.link(p, out.addLeaf(b));
  return p;
}

// Returns the OutTree node id for the Adams consensus subtree on the active leaf
// set, which is given as one firstVisit-ordered taxon list per input tree.
int adamsRecurse(const vector<InTree>& trees,
                 const vector<vector<int> >& leavesByTree, OutTree& out,
                 Scratch& sc) {
  int k = sc.k;
  const vector<int>& active = leavesByTree[0];
  int m = static_cast<int>(active.size());
  if (m == 1) return out.addLeaf(active[0]);
  if (m == 2) return makeCherry(out, active[0], active[1]);

  // Per tree: restricted subtree, centroid path, per-leaf branch data, and the
  // active leaves bucketed by branchLevel (ascending) for the spine walk.
  vector<vector<int> > heavyChild(k);
  vector<vector<int> > sortedByBL(k);
  for (int j = 0; j < k; ++j) {
    Nub nub = buildNub(trees[j], leavesByTree[j]);
    processNub(nub, sc.branchLevel[j], sc.sideChild[j], heavyChild[j]);
    // counting-sort active leaves by branchLevel[j]
    int len = static_cast<int>(heavyChild[j].size());  // >= maxBranchLevel + 2
    vector<int> cnt(len + 1, 0);
    for (int t : active) ++cnt[sc.branchLevel[j][t]];
    for (int w = 1; w < len; ++w) cnt[w] += cnt[w - 1];
    vector<int>& sj = sortedByBL[j];
    sj.assign(m, 0);
    // stable fill in firstVisit order is irrelevant; just place by bucket
    for (int t : active) sj[--cnt[sc.branchLevel[j][t]]] = t;
  }

  for (int t : active) {
    sc.removed[t] = 0;
    sc.blockOf[t] = -1;
  }

  // Spine walk: each step splits off the leaves branching at some tree's current
  // position (= the LCA depth of the remaining block in that tree), partitions
  // them by the meet key, and continues with the rest.
  vector<int> ptr(k, 0);
  vector<vector<int> > stepBlocks;   // stepBlocks[i] = side-block ids at step i
  int nBlocks = 0;
  int remaining = m;
  vector<int> Xlist;
  vector<int> idx;
  vector<int> keyFlat;
  vector<int> radixTmp;
  vector<int> radixCnt;
  vector<int> curPos(k, INT_MAX);
  while (remaining >= 2) {
    ++sc.stampCtr;
    Xlist.clear();
    for (int j = 0; j < k; ++j) {
      curPos[j] = INT_MAX;
      vector<int>& sj = sortedByBL[j];
      int& p = ptr[j];
      while (p < m && sc.removed[sj[p]]) ++p;
      if (p == m) continue;
      curPos[j] = sc.branchLevel[j][sj[p]];
      int q = p;
      while (q < m && sc.branchLevel[j][sj[q]] == curPos[j]) {
        int t = sj[q];
        if (!sc.removed[t] && sc.stamp[t] != sc.stampCtr) {
          sc.stamp[t] = sc.stampCtr;
          Xlist.push_back(t);
        }
        ++q;
      }
    }
    int sz = static_cast<int>(Xlist.size());
    // meet key for each branching leaf
    keyFlat.assign(static_cast<size_t>(sz) * k, 0);
    for (int i = 0; i < sz; ++i) {
      int t = Xlist[i];
      for (int j = 0; j < k; ++j)
        keyFlat[static_cast<size_t>(i) * k + j] =
            (sc.branchLevel[j][t] == curPos[j]) ? sc.sideChild[j][t]
                                                : heavyChild[j][curPos[j]];
    }
    idx.assign(sz, 0);
    for (int i = 0; i < sz; ++i) idx[i] = i;
    // Group identical meet-key vectors with an LSD radix over the k coordinates
    // (least-significant coordinate first), each pass a stable counting sort on
    // epoch-stamped, dense-compressed coordinate values.  O(k*sz) per step -- no
    // comparison sort, so the meet stays within the O(kn log n) bound.  Group
    // order is arbitrary but deterministic; the Adams tree is unique only up to
    // sibling order, which the clade-set contract ignores.
    radixTmp.assign(sz, 0);
    for (int j = k - 1; j >= 0; --j) {
      ++sc.compCtr;
      int d = 0;
      for (int i = 0; i < sz; ++i) {
        int v = keyFlat[static_cast<size_t>(idx[i]) * k + j];
        if (sc.compStamp[v] != sc.compCtr) {
          sc.compStamp[v] = sc.compCtr;
          sc.compId[v] = d++;
        }
      }
      radixCnt.assign(d + 1, 0);
      for (int i = 0; i < sz; ++i)
        ++radixCnt[sc.compId[keyFlat[static_cast<size_t>(idx[i]) * k + j]] + 1];
      for (int c = 0; c < d; ++c) radixCnt[c + 1] += radixCnt[c];
      for (int i = 0; i < sz; ++i) {
        int v = keyFlat[static_cast<size_t>(idx[i]) * k + j];
        radixTmp[radixCnt[sc.compId[v]]++] = idx[i];
      }
      idx.swap(radixTmp);
    }
    vector<int> theseBlocks;
    int i = 0;
    while (i < sz) {
      int run = i + 1;
      const int* ki = &keyFlat[static_cast<size_t>(idx[i]) * k];
      while (run < sz) {
        const int* kr = &keyFlat[static_cast<size_t>(idx[run]) * k];
        bool same = true;
        for (int j = 0; j < k; ++j)
          if (ki[j] != kr[j]) { same = false; break; }
        if (!same) break;
        ++run;
      }
      int bid = nBlocks++;
      theseBlocks.push_back(bid);
      for (int r = i; r < run; ++r) sc.blockOf[Xlist[idx[r]]] = bid;
      i = run;
    }
    stepBlocks.push_back(std::move(theseBlocks));
    for (int t : Xlist) sc.removed[t] = 1;
    remaining -= sz;
  }

  // The single leaf (if any) left on the spine bottom.
  int lastLeaf = -1;
  if (remaining == 1)
    for (int t : active)
      if (!sc.removed[t]) { lastLeaf = t; break; }

  // Distribute side-block leaves, preserving each tree's firstVisit order.
  vector<vector<vector<int> > > blockLeaves(
      nBlocks, vector<vector<int> >(k));
  for (int j = 0; j < k; ++j)
    for (int t : leavesByTree[j])
      if (sc.blockOf[t] >= 0) blockLeaves[sc.blockOf[t]][j].push_back(t);

  vector<int> blockNode(nBlocks, -1);
  for (int b = 0; b < nBlocks; ++b) {
    const vector<int>& bl = blockLeaves[b][0];
    int bs = static_cast<int>(bl.size());
    if (bs == 1)
      blockNode[b] = out.addLeaf(bl[0]);
    else if (bs == 2)
      blockNode[b] = makeCherry(out, bl[0], bl[1]);
    else
      blockNode[b] = adamsRecurse(trees, blockLeaves[b], out, sc);
  }

  // Assemble the spine as a nested chain (deepest step first), threading `deeper`.
  // Most step nodes have >= 2 children, but the deepest step can resolve to a
  // single block when no spine-bottom leaf is left (remaining reached 0): that
  // lone block is the whole subtree below, so its degree-1 node is suppressed and
  // passed through.  Rare but reachable (see the n=7/k=2 regression in
  // test-adams.R).  `kids` is never empty -- every spine step records a block.
  int deeper = (lastLeaf >= 0) ? out.addLeaf(lastLeaf) : -1;
  for (int s = static_cast<int>(stepBlocks.size()) - 1; s >= 0; --s) {
    vector<int> kids;
    for (int b : stepBlocks[s]) kids.push_back(blockNode[b]);
    if (deeper != -1) kids.push_back(deeper);
    if (kids.empty()) {
      // # nocov start
      // Unreachable: each spine step records >= 1 side block in stepBlocks[s].
      Rcpp::stop("cons_adams: empty assembly step");
      // # nocov end
    }
    if (kids.size() == 1) {
      deeper = kids[0];  // degree-1 suppression -- load-bearing, see above
    } else {
      int node = out.addInternal();
      for (int c : kids) out.link(node, c);
      deeper = node;
    }
  }
  return deeper;
}

std::string emitNewick(const OutTree& o, int root) {
  std::string s;
  struct Frame { int node; int ci; };
  vector<Frame> st;
  st.push_back(Frame{root, 0});
  while (!st.empty()) {
    int v = st.back().node;
    if (o.leafTaxon[v] > 0) {
      s += std::to_string(o.leafTaxon[v]);
      st.pop_back();
      continue;
    }
    int ci = st.back().ci;
    if (ci == 0) s += '(';
    if (ci < static_cast<int>(o.ch[v].size())) {
      if (ci > 0) s += ',';
      st.back().ci = ci + 1;
      st.push_back(Frame{o.ch[v][ci], 0});
    } else {
      s += ')';
      st.pop_back();
    }
  }
  return s;
}

}  // namespace

// Adams consensus of `edgeList` (each a PREORDER ape edge matrix on the tree's
// OWN root) on `nTip` leaves.  Returns an integer-label Newick string without a
// trailing ';'.  Implements Jansson, Li & Sung (2017), O(kn log n).
// [[Rcpp::export]]
std::string adamsConsensusCpp(Rcpp::List edgeList, int nTip) {
  int k = edgeList.size();
  vector<InTree> trees;
  trees.reserve(k);
  for (int i = 0; i < k; ++i) {
    Rcpp::IntegerMatrix edge = edgeList[i];
    fact::Tree t = fact::buildTreeFromEdge(edge, nTip);
    trees.push_back(buildInTree(t));
  }
  vector<vector<int> > leavesByTree(k);
  for (int j = 0; j < k; ++j) leavesByTree[j] = trees[j].taxaEulerOrder;

  OutTree out;
  Scratch sc(k, nTip);
  int root = adamsRecurse(trees, leavesByTree, out, sc);
  return emitNewick(out, root);
}
