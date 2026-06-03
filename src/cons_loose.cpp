#include <Rcpp.h>
#include "fact_tree.h"

#include <string>
#include <utility>
#include <vector>

// Loose (semi-strict / combinable-component) consensus, ported from FACT
// looseConsensusFast() + looseMerge() (dev/oracle/fact-src/loose.cpp:94) and
// contract() (dev/oracle/fact-src/strict.cpp:53), de-globalised, VLA-free and
// RAII (reusing fact::Tree / buildTreeFromEdge / newick / precompute from
// fact_tree.h).  Implements the asymptotically efficient algorithm of Jansson,
// Shen & Sung (2016): build a one-way-compatible tree by repeatedly merging the
// inputs (looseMerge with op == 1, each an O(n) consecutive-range query rather
// than an O(s^2) pairwise compatibility matrix), then mark which of its clusters
// are compatible with EVERY input (op == 0) and contract away the rest.
//
// This is a purely STRUCTURAL algorithm (Day's relabelling + consecutive-range
// /DEPTH queries); unlike greedy it does no leaf-set bit-packing, so there is no
// BUCKET_SIZE / word-index arithmetic here.
//
// IMPORTANT (rooting): the clusters are extracted as ROOTED clades, so the
// unrooted loose consensus is recovered only if every input tree is rooted
// consistently (FACT roots at the node adjacent to taxon 1).  The R wrapper
// Loose() roots each input at taxon 1 before extracting edges; do not call this
// on arbitrarily rooted trees.
//
// Loose's selection rule is deterministic (a split is kept iff it is compatible
// with every input -- no frequency tie-break), so the result is FACT-exact even
// at large n: any divergence from the fact.exe oracle is a real bug.

namespace {

using fact::Tree;

// Contract A down to its leaves plus the internal clusters flagged in A.good[]
// (good[label] marks the cluster A.precompute() hashed at Day's label `label`),
// rewiring each kept node to its nearest kept ancestor.  Ported verbatim from
// FACT strict.cpp:53; the unused `rooted` constructor argument is dropped and the
// VLAs / byte-correct memsets become std::vector fills.
Tree contract(Tree A) {
  std::vector<char> keep(A.cnt + 5, 0);
  std::vector<int> label(A.cnt + 5, -1), anc(A.cnt + 5, -1);
  int tmp = A.N;
  // Keep every leaf, and every good internal cluster.
  for (int i = 1; i <= A.N; ++i) {
    keep[A.idx[i]] = 1;
    if (A.good[i]) {
      keep[A.H[i]] = 1;
      ++tmp;
    }
  }
  Tree ret(A.N, tmp);
  tmp = 0;
  for (int i = 0; i < A.cnt; ++i) {
    if (keep[i]) label[i] = tmp++;
    if (A.leaf[i]) ret.leaf[label[i]] = A.leaf[i];
  }
  ret.root = label[A.root];
  // DFS: attach each kept vertex to the lowest kept ancestor of its parent.
  std::vector<std::pair<int, int> > S;
  S.push_back(std::make_pair(A.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (t1.first != A.root && keep[t1.first])
        ret.G[label[anc[t1.first]]].push_back(label[t1.first]);
      if (keep[t1.first]) anc[t1.first] = t1.first;
      else anc[t1.first] = anc[anc[t1.first]];
    }
    ++t1.second;
    if (t1.second < static_cast<int>(A.G[t1.first].size())) {
      S.push_back(t1);
      S.push_back(std::make_pair(A.G[t1.first][t1.second], -1));
      anc[A.G[t1.first][t1.second]] = anc[t1.first];
    }
  }
  for (int i = 0; i < ret.cnt; ++i)
    for (size_t c = 0; c < ret.G[i].size(); ++c) ret.parent[ret.G[i][c]] = i;
  ret.precompute();
  return ret;
}

// looseMerge(A, B, op): the consecutive-range compatibility query at the heart of
// the algorithm.  Relabels B's leaves by A's Day's labels and reorders B so each
// of A's subtrees would span a consecutive leaf range; then, for every internal
// cluster of A, decides in O(1) (after the linear setup) whether it is compatible
// with B.  op == 1 returns a new tree holding all of B plus the A-clusters
// compatible with B (the one-way-compatible refinement); op == 0 instead marks
// A.good[] for the A-clusters compatible with B and returns A.  A and B are taken
// BY VALUE -- both are mutated, and the callers' trees must be preserved.
// Ported verbatim from FACT loose.cpp:107 (globals -> arguments, VLAs ->
// std::vector, fixed stack -> std::vector stack).
Tree looseMerge(Tree A, Tree B, bool op, int numTaxas) {
  int newNodes = 0;
  std::vector<int> POS(numTaxas + 5, 0), DEPTH(B.cnt + 5, 0), pid(B.cnt + 5, 0);
  std::vector<char> can(A.cnt + 5, 0);
  // NB: FACT's `memset(A.size,0,(A.cnt+5))` here is a dead, byte-short write --
  // A.precompute() re-zeroes `size` (and everything else) immediately, and
  // nothing reads it in between, so we simply rely on precompute().
  // Relabel A's leaves so each subtree of A is a consecutive label range.
  A.precompute();

  for (int i = 1; i <= numTaxas; ++i) B.label[i] = A.label[i];

  std::vector<std::vector<int> > L(numTaxas + 5), R(numTaxas + 5),
      BEFORE(B.cnt + 5), AFTER(B.cnt + 5);
  std::vector<std::pair<int, int> > S;

  // Relabel B's leaves with A's labels and reorder B into an ordered tree
  // (children sorted by their minimum A-label).
  S.push_back(std::make_pair(B.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (B.leaf[t1.first] > 0) B.minL[t1.first] = B.label[B.leaf[t1.first]];
      else B.minL[t1.first] = numTaxas + 1;
    }
    ++t1.second;
    if (t1.second < static_cast<int>(B.G[t1.first].size())) {
      S.push_back(t1);
      S.push_back(std::make_pair(B.G[t1.first][t1.second], -1));
    } else {
      if (B.parent[t1.first] >= 0) L[B.minL[t1.first]].push_back(t1.first);
      B.G[t1.first].clear();
      if (B.parent[t1.first] >= 0) {
        if (B.minL[t1.first] < B.minL[B.parent[t1.first]])
          B.minL[B.parent[t1.first]] = B.minL[t1.first];
      }
    }
  }
  for (int i = 1; i <= numTaxas; ++i) {
    for (size_t k = 0; k < L[i].size(); ++k)
      B.G[B.parent[L[i][k]]].push_back(L[i][k]);
    L[i].clear();
  }

  // Construct the left (L) and right (R) paths and the DEPTH/POS/pid indices.
  int cnt = 0;
  POS[numTaxas + 1] = numTaxas + 1;
  POS[numTaxas + 2] = -1;
  DEPTH[B.root] = 0;
  S.push_back(std::make_pair(B.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (B.leaf[t1.first] > 0) {
        B.minL[t1.first] = B.maxL[t1.first] = B.leaf[t1.first];
        POS[B.leaf[t1.first]] = ++cnt;
      } else {
        B.minL[t1.first] = numTaxas + 1;
        B.maxL[t1.first] = numTaxas + 2;
      }
    }
    ++t1.second;
    if (t1.second < static_cast<int>(B.G[t1.first].size())) {
      DEPTH[B.G[t1.first][t1.second]] = DEPTH[t1.first] + 1;
      S.push_back(t1);
      S.push_back(std::make_pair(B.G[t1.first][t1.second], -1));
    } else {
      L[B.minL[t1.first]].push_back(t1.first);
      R[B.maxL[t1.first]].push_back(t1.first);
      for (size_t c = 0; c < B.G[t1.first].size(); ++c) {
        BEFORE[t1.first].push_back(0);
        AFTER[t1.first].push_back(0);
      }
      if (B.parent[t1.first] >= 0) {
        int p = B.parent[t1.first];
        if (POS[B.minL[t1.first]] < POS[B.minL[p]]) B.minL[p] = B.minL[t1.first];
        if (POS[B.maxL[t1.first]] > POS[B.maxL[p]]) B.maxL[p] = B.maxL[t1.first];
      }
      int cur = 0;
      for (size_t c = 0; c < B.G[t1.first].size(); ++c) pid[B.G[t1.first][c]] = cur++;
    }
  }

  // Query every subtree of A: is its leaf set a cluster of (the ordered) B?
  S.push_back(std::make_pair(A.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (A.leaf[t1.first]) A.minL[t1.first] = A.maxL[t1.first] = A.leaf[t1.first];
      else {
        A.minL[t1.first] = numTaxas + 1;
        A.maxL[t1.first] = numTaxas + 2;
      }
    }
    ++t1.second;
    if (t1.second < static_cast<int>(A.G[t1.first].size())) {
      S.push_back(t1);
      S.push_back(std::make_pair(A.G[t1.first][t1.second], -1));
    } else {
      int a = A.minL[t1.first], b = A.maxL[t1.first];
      bool skip = false;
      if (POS[b] - POS[a] + 1 != A.size[t1.first]) {
        can[t1.first] = 0;
        skip = true;
      }
      if (!skip) {
        if (DEPTH[L[a].back()] >= DEPTH[R[b].back()]) {
          if (DEPTH[L[a].back()] - DEPTH[B.idx[b]] > 0) {
            can[t1.first] = 0;
            skip = true;
          } else {
            a = L[a].back();
            // Guard is not reachable via the public API (R wrapper validates
            // input), but defence-in-depth against future refactors.
            int ridx = DEPTH[B.idx[b]] - DEPTH[a];
            if (ridx < 0 || ridx >= static_cast<int>(R[b].size()))
              Rcpp::stop("Internal error: node ID %d out of range (depth vector size %d)",
                         ridx, (int)R[b].size());
            b = R[b][ridx];
          }
        } else {
          if (DEPTH[R[b].back()] - DEPTH[B.idx[a]] > 0) {
            can[t1.first] = 0;
            skip = true;
          } else {
            b = R[b].back();
            // Guard is not reachable via the public API (R wrapper validates
            // input), but defence-in-depth against future refactors.
            int lidx = DEPTH[B.idx[a]] - DEPTH[b];
            if (lidx < 0 || lidx >= static_cast<int>(L[a].size()))
              Rcpp::stop("Internal error: node ID %d out of range (depth vector size %d)",
                         lidx, (int)L[a].size());
            a = L[a][lidx];
          }
        }
      }
      if (!skip) {
        if (a == b) {
          can[t1.first] = 1;
        } else if (B.parent[a] == B.parent[b]) {
          can[t1.first] = 1;
          int q1 = pid[a], q2 = pid[b];
          ++newNodes;
          ++BEFORE[B.parent[a]][q1];
          ++AFTER[B.parent[a]][q2];
        } else {
          can[t1.first] = 0;
        }
      }
      if (A.parent[t1.first] >= 0) {
        int p = A.parent[t1.first];
        if (POS[A.minL[t1.first]] < POS[A.minL[p]]) A.minL[p] = A.minL[t1.first];
        if (POS[A.maxL[t1.first]] > POS[A.maxL[p]]) A.maxL[p] = A.maxL[t1.first];
      }
    }
  }

  if (op == 0) {
    // Mark/unmark A's clusters: good[i] stays set only where compatible with B.
    for (int i = 1; i <= numTaxas; ++i)
      if (A.H[i] >= 0) A.good[i] = static_cast<char>(A.good[i] & can[A.H[i]]);
    return A;
  }
  // op == 1: build a new tree holding all of B plus the new (BEFORE/AFTER) nodes.
  if (newNodes == 0) return B;
  Tree ret(numTaxas, B.cnt + newNodes);
  ret.root = B.root;
  ret.parent[ret.root] = -1;
  int next = B.cnt, cur;
  for (int a = 0; a < B.cnt; ++a) {
    ret.leaf[a] = B.leaf[a];
    cur = a;
    for (int i = 0; i < static_cast<int>(B.G[a].size()); ++i) {
      while (BEFORE[a][i]) {
        ret.parent[next] = cur;
        ret.G[cur].push_back(next);
        cur = next++;
        --BEFORE[a][i];
      }
      ret.G[cur].push_back(B.G[a][i]);
      ret.parent[B.G[a][i]] = cur;
      while (AFTER[a][i]) {
        cur = ret.parent[cur];
        --AFTER[a][i];
      }
    }
  }
  // FACT loose.cpp:258 invariant: exactly `newNodes` vertices were inserted, so
  // the running index must have reached the allocated node count.  A mismatch
  // means the BEFORE/AFTER bookkeeping built the wrong shape -- fail loud rather
  // than return a silently mis-built consensus.  Cannot be triggered via a
  // legitimate R call (looseMerge is internal; the check fires only if the C++
  // bookkeeping itself is wrong), so excluded from coverage like analogous
  // post-condition guards in rstar.cpp / bhv.cpp.
  if (next != ret.cnt) {
    // # nocov start
    Rcpp::stop("looseMerge: built %d nodes but allocated %d", next, ret.cnt);
    // # nocov end
  }
  for (int i = 1; i <= numTaxas; ++i) ret.idx[i] = B.idx[i];
  return ret;
}

// looseConsensusFast (loose.cpp:94): merge the inputs into a one-way-compatible
// tree, mark the clusters compatible with every input, then contract.
Tree looseConsensus(const std::vector<Tree>& T, int numTaxas) {
  const int numTrees = static_cast<int>(T.size());
  Tree ret;
  for (int i = 0; i < numTrees; ++i) {
    if (i == 0) ret = T[0];
    else ret = looseMerge(ret, T[i], 1, numTaxas);  // one-way-compatible tree
  }
  ret.precompute();
  for (int i = 1; i <= numTaxas; ++i) ret.good[i] = (ret.H[i] >= 0);
  // Determine which vertices of the one-way-compatible tree lie in the loose
  // consensus (compatible with EVERY input).
  for (int i = 0; i < numTrees; ++i) ret = looseMerge(ret, T[i], 0, numTaxas);
  return contract(ret);
}

}  // namespace

// Loose consensus of `edgeList` (each a PREORDER ape edge matrix, rooted at
// taxon 1 by the caller) on `nTip` leaves.  Returns an integer-label Newick
// string without a trailing ';'.
// [[Rcpp::export]]
std::string looseConsensusCpp(Rcpp::List edgeList, int nTip) {
  int nTree = edgeList.size();
  std::vector<fact::Tree> T;
  T.reserve(nTree);
  for (int i = 0; i < nTree; ++i) {
    Rcpp::IntegerMatrix edge = edgeList[i];
    T.push_back(fact::buildTreeFromEdge(edge, nTip));
  }
  fact::Tree ret = looseConsensus(T, nTip);
  return ret.newick();
}
