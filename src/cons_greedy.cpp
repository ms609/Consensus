#include <Rcpp.h>
#include "fact_tree.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

// Greedy (extended majority-rule) consensus, ported from FACT
// greedyConsensusFast() (dev/oracle/fact-src/greedy.cpp:134), de-globalised and
// VLA-free.  Distinct clusters are extracted as packed bit-vectors in one
// post-order sweep per tree, grouped by occurrence count, then added in
// descending-frequency order whenever compatible with the tree built so far --
// the asymptotically efficient algorithm of Jansson, Shen & Sung (2016), in
// place of the previous R O(s^2) compatibility matrix.
//
// IMPORTANT (rooting): the algorithm extracts ROOTED clusters, so the unrooted
// greedy consensus is recovered only if every input tree is rooted consistently
// (FACT roots at the node adjacent to taxon 1).  The R wrapper Greedy() roots
// each input at taxon 1 before extracting edges; do not call this on arbitrarily
// rooted trees.

namespace {

typedef long long ll;
const int BUCKET_SIZE = 60;  // bits packed per word (matches FACT)

using fact::Tree;

fact::Tree greedyConsensus(const std::vector<Tree>& T, int numTaxas) {
  const int numTrees = static_cast<int>(T.size());
  const int LEN = (numTaxas + BUCKET_SIZE - 1) / BUCKET_SIZE;

  // Distinct internal clusters across all trees, as packed leaf-set bit-vectors.
  std::vector<std::vector<ll> > LF;

  // ---- extract every internal cluster from every tree -----------------------
  {
    std::vector<ll> cluster(static_cast<size_t>(2 * numTaxas + 5) * LEN, 0);
    std::vector<std::pair<int, int> > S;
    for (int i = 0; i < numTrees; ++i) {
      const Tree& Ti = T[i];
      for (int j = 0; j < Ti.cnt; ++j)
        for (int k = 0; k < LEN; ++k)
          cluster[static_cast<size_t>(j) * LEN + k] = 0;
      S.clear();
      S.push_back(std::make_pair(Ti.root, -1));
      while (!S.empty()) {
        std::pair<int, int> t1 = S.back();
        S.pop_back();
        ++t1.second;
        if (t1.second < static_cast<int>(Ti.G[t1.first].size())) {
          S.push_back(t1);
          S.push_back(std::make_pair(Ti.G[t1.first][t1.second], -1));
        } else {
          if (Ti.leaf[t1.first] > 0) {
            int c = Ti.leaf[t1.first];
            cluster[static_cast<size_t>(t1.first) * LEN + (c - 1) / BUCKET_SIZE]
                |= (1LL << ((c - 1) % BUCKET_SIZE));
          } else if (t1.first != Ti.root) {
            std::vector<ll> v(LEN);
            for (int k = 0; k < LEN; ++k)
              v[k] = cluster[static_cast<size_t>(t1.first) * LEN + k];
            LF.push_back(std::move(v));
          }
          if (Ti.parent[t1.first] >= 0) {
            int p = Ti.parent[t1.first];
            for (int k = 0; k < LEN; ++k)
              cluster[static_cast<size_t>(p) * LEN + k]
                  |= cluster[static_cast<size_t>(t1.first) * LEN + k];
          }
        }
      }
    }
  }

  const int L = static_cast<int>(LF.size());

  // bit x (1-indexed taxon) set in cluster LF[a]?
  auto bitExist = [&](int a, int x) -> bool {
    return (LF[a][(x - 1) / BUCKET_SIZE] & (1LL << ((x - 1) % BUCKET_SIZE))) != 0;
  };
  auto cmpLeafSet = [&](int a, int b) -> bool {
    for (int k = 0; k < LEN; ++k)
      if (LF[a][k] != LF[b][k]) return LF[a][k] < LF[b][k];
    return false;
  };
  auto sameLeafSet = [&](int a, int b) -> bool {
    for (int k = 0; k < LEN; ++k)
      if (LF[a][k] != LF[b][k]) return false;
    return true;
  };

  // ---- bucket distinct clusters by occurrence count -------------------------
  // After sorting, identical clusters are adjacent; each run's length is the
  // number of trees that display that cluster.  CountingSort[c] lists the
  // (representative indices of) clusters occurring in exactly c trees.
  std::vector<std::vector<int> > CountingSort(numTrees + 5);
  if (L > 0) {
    std::vector<int> A(L);
    for (int i = 0; i < L; ++i) A[i] = i;
    std::sort(A.begin(), A.end(), cmpLeafSet);
    int run = 0;
    for (int i = 0; i < L; ++i) {
      if (sameLeafSet(A[run], A[i])) continue;
      CountingSort[i - run].push_back(A[run]);
      run = i;
    }
    CountingSort[L - run].push_back(A[run]);
  }

  // ---- greedily add compatible clusters, most frequent first ----------------
  Tree ret(numTaxas, numTaxas + 1);
  ret.root = 0;
  for (int i = 1; i <= numTaxas; ++i) {
    ret.leaf[i] = i;
    ret.parent[i] = 0;
    ret.G[0].push_back(i);
  }
  int accepted = 1;
  std::vector<std::pair<int, int> > S;
  for (int a = numTrees; a > 0; --a) {
    for (size_t b = 0; b < CountingSort[a].size(); ++b) {
      int tmp = CountingSort[a][b];
      int sze = 0;
      for (int k = 0; k < LEN; ++k)
        for (int j = 0; j < BUCKET_SIZE; ++j)
          sze += (LF[tmp][k] & (1LL << j)) != 0;
      S.clear();
      S.push_back(std::make_pair(ret.root, -1));
      while (!S.empty()) {
        std::pair<int, int> t1 = S.back();
        S.pop_back();
        if (t1.second == -1) {
          ret.size[t1.first] = ret.minL[t1.first] = 0;
          if (ret.leaf[t1.first] > 0) {
            ++ret.size[t1.first];
            if (bitExist(tmp, ret.leaf[t1.first])) ++ret.minL[t1.first];
          }
        }
        ++t1.second;
        if (t1.second < static_cast<int>(ret.G[t1.first].size())) {
          S.push_back(t1);
          S.push_back(std::make_pair(ret.G[t1.first][t1.second], -1));
        } else {
          if (ret.minL[t1.first] == sze) {
            bool clean = true;
            for (size_t c = 0; c < ret.G[t1.first].size(); ++c) {
              int ch = ret.G[t1.first][c];
              clean = clean && (ret.minL[ch] == 0 || ret.minL[ch] == ret.size[ch]);
            }
            if (clean) {
              Tree newT(numTaxas, ret.cnt + 1);
              newT.root = 0;
              for (int i = 1; i <= numTaxas; ++i) newT.leaf[i] = i;
              for (int i = 0; i < ret.cnt; ++i) {
                if (i == t1.first) continue;
                for (size_t c = 0; c < ret.G[i].size(); ++c) {
                  int ch = ret.G[i][c];
                  newT.G[i].push_back(ch);
                  newT.parent[ch] = i;
                }
              }
              for (size_t c = 0; c < ret.G[t1.first].size(); ++c) {
                int ch = ret.G[t1.first][c];
                if (ret.minL[ch] == ret.size[ch]) {
                  newT.G[ret.cnt].push_back(ch);
                  newT.parent[ch] = ret.cnt;
                } else {
                  newT.G[t1.first].push_back(ch);
                  newT.parent[ch] = t1.first;
                }
              }
              newT.G[t1.first].push_back(ret.cnt);
              newT.parent[ret.cnt] = t1.first;
              ret = newT;
              ++accepted;
              if (accepted == numTaxas) return ret;
            }
            break;
          }
          if (ret.parent[t1.first] >= 0) {
            int p = ret.parent[t1.first];
            ret.minL[p] += ret.minL[t1.first];
            ret.size[p] += ret.size[t1.first];
          }
        }
      }
    }
  }
  return ret;
}

}  // namespace

// Greedy consensus of `edgeList` (each a PREORDER ape edge matrix, rooted at
// taxon 1 by the caller) on `nTip` leaves.  Returns an integer-label Newick
// string without a trailing ';'.
// [[Rcpp::export]]
std::string greedyConsensusCpp(Rcpp::List edgeList, int nTip) {
  int nTree = edgeList.size();
  std::vector<fact::Tree> T;
  T.reserve(nTree);
  for (int i = 0; i < nTree; ++i) {
    Rcpp::IntegerMatrix edge = edgeList[i];
    T.push_back(fact::buildTreeFromEdge(edge, nTip));
  }
  fact::Tree ret = greedyConsensus(T, nTip);
  return ret.newick();
}
