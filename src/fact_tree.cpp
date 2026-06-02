#include "fact_tree.h"

#include <algorithm>
#include <stack>
#include <utility>

namespace fact {

Tree::Tree(int N_, int cnt_) : N(N_), cnt(cnt_) {
  idx.assign(N + 5, 0);
  label.assign(N + 5, 0);
  minH.assign(N + 5, 0);
  maxH.assign(N + 5, 0);
  H.assign(N + 5, -1);
  good.assign(N + 5, 0);
  leaf.assign(cnt + 5, 0);
  size.assign(cnt + 5, 0);
  minL.assign(cnt + 5, 0);
  maxL.assign(cnt + 5, 0);
  parent.assign(cnt + 5, -1);
  goodLabel.assign(cnt + 5, -1);
  G.assign(cnt + 5, std::vector<int>());
}

// Day's algorithm: relabel the leaves by an in-order numbering so that every
// subtree spans a contiguous [minL, maxL] range, and hash each internal cluster
// by (minH, maxH, H).  Ported verbatim from FACT tree::precompute()
// (dev/oracle/fact-src/tree.cpp:86); globals/VLAs removed.
//
// NOTE: precompute() is shared scaffolding for the Loose and MajorityPlus chips
// (Greedy does not use it).  It is first exercised against the FACT oracle by
// those chips; it is shipped here so both chips can call it without each adding
// a colliding definition.
void Tree::precompute() {
  std::stack<std::pair<int, int> > S;
  std::vector<char> labelled(N + 5, 0);
  std::fill(label.begin(), label.end(), 0);
  std::fill(size.begin(), size.end(), 0);
  std::fill(H.begin(), H.end(), -1);
  for (int i = 0; i < cnt; ++i) {
    size[i] = 0;
    minL[i] = cnt + 5;
    maxL[i] = -1;
    if (leaf[i] > 0) idx[leaf[i]] = i;
  }
  int num = 0;
  S.push(std::make_pair(root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.top();
    S.pop();
    if (leaf[t1.first] > 0 && !labelled[leaf[t1.first]]) {
      label[leaf[t1.first]] = minL[t1.first] = maxL[t1.first] = ++num;
      size[t1.first] = 1;
      labelled[leaf[t1.first]] = 1;
    }
    ++t1.second;
    if (t1.second < (int)G[t1.first].size()) {
      S.push(std::make_pair(t1.first, t1.second));
      S.push(std::make_pair(G[t1.first][t1.second], -1));
    } else {
      if (!leaf[t1.first]) {
        if (!S.empty() && S.top().second == 0) {
          minH[maxL[t1.first]] = minL[t1.first];
          maxH[maxL[t1.first]] = maxL[t1.first];
          H[maxL[t1.first]] = t1.first;
        } else {
          minH[minL[t1.first]] = minL[t1.first];
          maxH[minL[t1.first]] = maxL[t1.first];
          H[minL[t1.first]] = t1.first;
        }
      }
      if (!S.empty()) {
        size[S.top().first] += size[t1.first];
        minL[S.top().first] = std::min(minL[S.top().first], minL[t1.first]);
        maxL[S.top().first] = std::max(maxL[S.top().first], maxL[t1.first]);
      }
    }
  }
  // FACT asserts num == N here; on a well-formed tree this always holds.
}

// Emit the tree as an integer-label Newick string WITHOUT a trailing ';'
// (matching src/rstar.cpp / src/local_consensus.cpp; the R side appends ';').
// Ported from FACT tree::printNex() (dev/oracle/fact-src/tree.cpp:136), but
// returns a std::string instead of writing to std::cout.
std::string Tree::newick() const {
  std::stack<std::pair<int, int> > S;
  std::string ans;
  S.push(std::make_pair(root, -1));
  while (!S.empty()) {
    std::pair<int, int> t1 = S.top();
    S.pop();
    if (leaf[t1.first]) {
      ans += std::to_string(leaf[t1.first]);
      continue;
    }
    if (t1.second == -1) ans += '(';
    ++t1.second;
    if (t1.second < (int)G[t1.first].size()) {
      if (t1.second > 0) ans += ',';
      S.push(t1);
      S.push(std::make_pair(G[t1.first][t1.second], -1));
    } else {
      ans += ')';
    }
  }
  return ans;
}

Tree buildTreeFromEdge(const Rcpp::IntegerMatrix& edge, int nTip) {
  int nRow = edge.nrow();
  int nNode = 0;
  for (int r = 0; r < nRow; ++r) {
    if (edge(r, 0) > nNode) nNode = edge(r, 0);
    if (edge(r, 1) > nNode) nNode = edge(r, 1);
  }
  Tree t(nTip, nNode);
  // Tips: ape node v (1..nTip) -> FACT node v-1, carrying taxon v.
  for (int v = 1; v <= nTip; ++v) {
    t.leaf[v - 1] = v;
    t.idx[v] = v - 1;
  }
  for (int r = 0; r < nRow; ++r) {
    int p = edge(r, 0) - 1;
    int c = edge(r, 1) - 1;
    t.G[p].push_back(c);
    t.parent[c] = p;
  }
  // Preorder: the first edge emanates from the root.
  t.root = nRow > 0 ? edge(0, 0) - 1 : 0;
  t.parent[t.root] = -1;
  return t;
}

}  // namespace fact
