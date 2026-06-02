#include <Rcpp.h>
#include "fact_tree.h"

#include <string>
#include <utility>
#include <vector>

// Majority-rule (+) consensus, ported from FACT majorityPlusConsensus()
// (dev/oracle/fact-src/majorityplus.cpp:6) together with its helpers
// updateCounter() and majorityPlusMerge() (same file) and majContract()
// (dev/oracle/fact-src/majority.cpp:169).  De-globalised, VLA-free and RAII, in
// the same style as the Greedy worked example (src/cons_greedy.cpp):
//
//   * NO file-scope globals.  FACT threads `numTaxas`/`numTrees`/`rooted` and
//     `tree *T` through externs; here every routine takes its state by argument
//     (numTaxas is read off Tree::N, the tree list is passed in), so the code is
//     reentrant when called repeatedly from R.
//   * NO C99 variable-length arrays.  FACT's `pair<int,int> S[..]`,
//     `int POS[..]`, `vector<int> L[..]` etc. become std::vector, runtime-sized;
//     the explicit-index stacks become push_back/back/pop_back (as in greedy).
//   * RAII.  fact::Tree owns its arrays and is copyable/movable, so the
//     by-value tree arguments and `ret = ...` reassignments are self-freeing.
//   * `goto end` (skip the compatibility verdict, always run the parent
//     propagation) is rewritten with a `decided` flag.
//   * The unused `bool can[..]` scratch and the dead `memset(A.size,0,A.cnt+5)`
//     (a no-op even in FACT -- it zeroes only cnt+5 *bytes* and is immediately
//     overwritten by precompute()) are dropped; FACT's `assert(cnt==ret.cnt)`
//     becomes a comment (it holds on well-formed input).
//
// Unlike Greedy there is NO bit-packing here (no BUCKET_SIZE): compatibility is
// decided by Day's leaf relabelling (Tree::precompute() supplies `label`,
// `size`, `idx`) plus left/right path queries on an ordered copy of the second
// tree.  MajorityPlus keeps a cluster iff it is displayed by strictly more input
// trees than contradict it -- a deterministic count rule with no frequency
// tie-break, so the result is FACT-exact (the asymptotically efficient O(kn)
// algorithm of Jansson, Shen & Sung 2016).
//
// IMPORTANT (rooting): the algorithm reasons over ROOTED clusters; the unrooted
// majority-rule (+) consensus is recovered because every input is rooted
// consistently.  The R wrapper MajorityPlus() roots each input at taxon 1 before
// extracting edges; do not call this on arbitrarily rooted trees.

namespace {

using fact::Tree;

// FACT majContract (majority.cpp:169): keep every node whose counter
// (`goodLabel`) is positive, renumbering the survivors into a fresh Tree.
Tree majContract(Tree X) {
  const int cnt = X.cnt;
  std::vector<int> label(cnt + 5, -1), anc(cnt + 5, -1);
  int tmp = 0;
  for (int i = 0; i < cnt; ++i) tmp += (X.goodLabel[i] > 0);
  Tree ret(X.N, tmp);
  tmp = 0;
  for (int i = 0; i < cnt; ++i) {
    if (X.goodLabel[i] > 0) {
      label[i] = tmp++;
      ret.goodLabel[label[i]] = X.goodLabel[i];
    }
    if (X.leaf[i] > 0) {
      ret.leaf[label[i]] = X.leaf[i];
      ret.idx[X.leaf[i]] = label[i];
    }
  }
  ret.root = label[X.root];
  std::vector<std::pair<int, int> > S;
  S.push_back(std::make_pair(X.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (X.goodLabel[t1.first] > 0 && anc[t1.first] >= 0)
        ret.G[label[anc[t1.first]]].push_back(label[t1.first]);
      if (X.goodLabel[t1.first] > 0) anc[t1.first] = t1.first;
    }
    ++t1.second;
    if (t1.second < static_cast<int>(X.G[t1.first].size())) {
      S.push_back(t1);
      S.push_back(std::make_pair(X.G[t1.first][t1.second], -1));
      anc[X.G[t1.first][t1.second]] = anc[t1.first];
    }
  }
  for (int i = 0; i < ret.cnt; ++i)
    for (std::vector<int>::iterator it = ret.G[i].begin();
         it != ret.G[i].end(); ++it)
      ret.parent[*it] = i;
  return ret;
}

// FACT updateCounter (majorityplus.cpp:26): for every cluster of A, A's counter
// (A.goodLabel) is incremented if B displays the cluster, decremented if B
// contradicts it, and left unchanged if B is merely compatible.  A is updated in
// place; B is consumed (relabelled and reordered into an ordered tree).
void updateCounter(Tree& A, Tree B) {
  const int numTaxas = A.N;
  std::vector<int> POS(numTaxas + 5, 0);
  std::vector<int> DEPTH(B.cnt + 5, 0);
  // Relabel A's leaves so each subtree spans a consecutive range (Day's).
  A.precompute();
  for (int i = 1; i <= numTaxas; ++i) B.label[i] = A.label[i];

  std::vector<std::vector<int> > L(numTaxas + 5), R(numTaxas + 5);
  std::vector<std::vector<int> > BEFORE(B.cnt + 5), AFTER(B.cnt + 5);
  std::vector<std::pair<int, int> > S;

  // Relabel B's leaves with A's labels and reorder B into an ordered tree.
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
  for (int i = 1; i <= numTaxas; L[i].clear(), ++i)
    for (std::vector<int>::iterator it = L[i].begin(); it != L[i].end(); ++it)
      B.G[B.parent[*it]].push_back(*it);

  // Build the left/right paths, the position index POS and the depths DEPTH.
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
    }
  }

  // Query each subtree of A against B and adjust A's counter.
  S.push_back(std::make_pair(A.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (A.leaf[t1.first]) {
        A.minL[t1.first] = A.maxL[t1.first] = A.leaf[t1.first];
      } else {
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
      bool decided = false;
      if (POS[b] - POS[a] + 1 != A.size[t1.first]) {
        --A.goodLabel[t1.first];
        decided = true;
      }
      if (!decided) {
        if (DEPTH[L[a].back()] >= DEPTH[R[b].back()]) {
          if (DEPTH[L[a].back()] - DEPTH[B.idx[b]] > 0) {
            --A.goodLabel[t1.first];
            decided = true;
          } else {
            a = L[a].back();
            b = R[b][DEPTH[B.idx[b]] - DEPTH[a]];
          }
        } else {
          if (DEPTH[R[b].back()] - DEPTH[B.idx[a]] > 0) {
            --A.goodLabel[t1.first];
            decided = true;
          } else {
            b = R[b].back();
            a = L[a][DEPTH[B.idx[a]] - DEPTH[b]];
          }
        }
      }
      if (!decided) {
        if (a == b) ++A.goodLabel[t1.first];
        else if (B.parent[a] != B.parent[b]) --A.goodLabel[t1.first];
      }
      if (A.parent[t1.first] >= 0) {
        int p = A.parent[t1.first];
        if (POS[A.minL[t1.first]] < POS[A.minL[p]]) A.minL[p] = A.minL[t1.first];
        if (POS[A.maxL[t1.first]] > POS[A.maxL[p]]) A.maxL[p] = A.maxL[t1.first];
      }
    }
  }
}

// FACT majorityPlusMerge (majorityplus.cpp:138): insert into B every cluster of
// A that is compatible with B but not yet present (a refinement), giving each
// inserted node counter 1, then contract.  Counters are NOT touched for clusters
// already in B (the `can`/goodLabel adjustment of majorityMerge is intentionally
// absent here -- updateCounter does the counting separately).
Tree majorityPlusMerge(Tree A, Tree B) {
  const int numTaxas = A.N;
  std::vector<int> POS(numTaxas + 5, 0);
  std::vector<int> DEPTH(B.cnt + 5, 0);
  std::vector<int> pid(B.cnt + 5, 0);
  int newNodes = 0;
  A.precompute();
  for (int i = 1; i <= numTaxas; ++i) B.label[i] = A.label[i];

  std::vector<std::vector<int> > L(numTaxas + 5), R(numTaxas + 5);
  std::vector<std::vector<int> > BEFORE(B.cnt + 5), AFTER(B.cnt + 5);
  std::vector<std::pair<int, int> > S;

  // Relabel B's leaves with A's labels and reorder B into an ordered tree.
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
  for (int i = 1; i <= numTaxas; L[i].clear(), ++i)
    for (std::vector<int>::iterator it = L[i].begin(); it != L[i].end(); ++it)
      B.G[B.parent[*it]].push_back(*it);

  // Build the left/right paths, POS, DEPTH and the per-parent child index pid.
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
      for (size_t c = 0; c < B.G[t1.first].size(); ++c)
        pid[B.G[t1.first][c]] = static_cast<int>(c);
    }
  }

  // Query each subtree of A; for an insertable refinement, record where a new
  // node must be threaded into B (BEFORE/AFTER counts on its parent's children).
  S.push_back(std::make_pair(A.root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.back();
    S.pop_back();
    if (t1.second == -1) {
      if (A.leaf[t1.first]) {
        A.minL[t1.first] = A.maxL[t1.first] = A.leaf[t1.first];
      } else {
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
      bool decided = false;
      if (POS[b] - POS[a] + 1 != A.size[t1.first]) decided = true;
      if (!decided) {
        if (DEPTH[L[a].back()] >= DEPTH[R[b].back()]) {
          if (DEPTH[L[a].back()] - DEPTH[B.idx[b]] > 0) decided = true;
          else {
            a = L[a].back();
            b = R[b][DEPTH[B.idx[b]] - DEPTH[a]];
          }
        } else {
          if (DEPTH[R[b].back()] - DEPTH[B.idx[a]] > 0) decided = true;
          else {
            b = R[b].back();
            a = L[a][DEPTH[B.idx[a]] - DEPTH[b]];
          }
        }
      }
      if (!decided && a != b && B.parent[a] == B.parent[b]) {
        int pa = pid[a], pb = pid[b];
        ++newNodes;
        ++BEFORE[B.parent[a]][pa];
        ++AFTER[B.parent[a]][pb];
      }
      if (A.parent[t1.first] >= 0) {
        int p = A.parent[t1.first];
        if (POS[A.minL[t1.first]] < POS[A.minL[p]]) A.minL[p] = A.minL[t1.first];
        if (POS[A.maxL[t1.first]] > POS[A.maxL[p]]) A.maxL[p] = A.maxL[t1.first];
      }
    }
  }

  Tree ret;
  if (newNodes > 0) {
    ret = Tree(numTaxas, B.cnt + newNodes);
    ret.root = B.root;
    ret.goodLabel[ret.root] = B.goodLabel[B.root];
    ret.parent[ret.root] = -1;
    int cnt2 = B.cnt, cur;
    for (int a = 0; a < B.cnt; ++a) {
      ret.leaf[a] = B.leaf[a];
      cur = a;
      for (int i = 0; i < static_cast<int>(B.G[a].size()); ++i) {
        while (BEFORE[a][i]) {
          ret.parent[cnt2] = cur;
          ret.G[cur].push_back(cnt2);
          ret.goodLabel[cnt2] = 1;
          cur = cnt2++;
          --BEFORE[a][i];
        }
        ret.G[cur].push_back(B.G[a][i]);
        ret.goodLabel[B.G[a][i]] = B.goodLabel[B.G[a][i]];
        ret.parent[B.G[a][i]] = cur;
        while (AFTER[a][i]) {
          cur = ret.parent[cur];
          --AFTER[a][i];
        }
      }
    }
    for (int i = 1; i <= numTaxas; ++i) ret.idx[i] = B.idx[i];
    // FACT asserts cnt2 == ret.cnt here; it holds on well-formed input.
  } else {
    ret = B;
  }
  ret = majContract(ret);
  return ret;
}

// FACT majorityPlusConsensus (majorityplus.cpp:6): build the union of all
// candidate clusters by repeated merge while streaming the per-cluster counter,
// then recount the survivors over every tree and contract to those displayed by
// strictly more trees than contradict them.
Tree majorityPlusConsensus(const std::vector<Tree>& T, int numTaxas) {
  const int numTrees = static_cast<int>(T.size());
  Tree ret;
  for (int i = 0; i < numTrees; ++i) {
    if (i == 0) {
      ret = T[i];
      for (int j = 0; j < ret.cnt; ++j) ret.goodLabel[j] = 1;
    } else {
      updateCounter(ret, T[i]);
      ret = majorityPlusMerge(T[i], ret);
    }
  }
  for (int i = 0; i < ret.cnt; ++i) ret.goodLabel[i] = 0;
  for (int i = 0; i < numTrees; ++i) updateCounter(ret, T[i]);
  ret = majContract(ret);
  (void) numTaxas;  // numTaxas is carried on each Tree::N; kept for parity.
  return ret;
}

}  // namespace

// Majority-rule (+) consensus of `edgeList` (each a PREORDER ape edge matrix,
// rooted at taxon 1 by the caller) on `nTip` leaves.  Returns an integer-label
// Newick string without a trailing ';'.
// [[Rcpp::export]]
std::string majorityPlusConsensusCpp(Rcpp::List edgeList, int nTip) {
  int nTree = edgeList.size();
  std::vector<fact::Tree> T;
  T.reserve(nTree);
  for (int i = 0; i < nTree; ++i) {
    Rcpp::IntegerMatrix edge = edgeList[i];
    T.push_back(fact::buildTreeFromEdge(edge, nTip));
  }
  fact::Tree ret = majorityPlusConsensus(T, nTip);
  return ret.newick();
}
