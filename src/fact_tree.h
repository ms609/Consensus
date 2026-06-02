#ifndef CONSTREE_FACT_TREE_H
#define CONSTREE_FACT_TREE_H

// Reentrant, allocation-safe C++ tree primitive ported from the FACT toolkit of
// Jansson and colleagues (dev/oracle/fact-src/tree.{h,cpp}; used with
// permission).  It is the shared substrate for the fast split-selection
// consensus methods -- Greedy (this foundation), and Loose / MajorityPlus
// (added by their chips, which port their method-specific machinery into their
// own src/cons_<method>.cpp but reuse Tree, buildTreeFromEdge, newick() and
// precompute() from here).
//
// Three deliberate departures from the FACT original, each required to run the
// code repeatedly and portably from inside an R package:
//   * NO file-scope globals.  FACT threads state through `extern tree *T; extern
//     int numTaxas,...`; here every algorithm takes its state by argument, so
//     the functions are reentrant.
//   * NO C99 variable-length arrays.  All storage is std::vector, runtime-sized
//     (VLAs are non-portable; MSVC and `R CMD check` reject them).
//   * RAII.  The FACT `struct tree` new[]s ~11 raw arrays with no destructor and
//     is copied by value throughout the algorithms; here every array is a
//     std::vector member, so Tree is copyable, movable and self-freeing.
//
// Index conventions inherited verbatim from FACT (so the ported algorithms read
// the same): taxa are 1-indexed, nodes are 0-indexed; idx[taxon] = node,
// leaf[node] = taxon (0 for internal nodes).  Per-taxon arrays are sized N + 5,
// per-node arrays cnt + 5, matching FACT's guard padding.

#include <Rcpp.h>
#include <string>
#include <vector>

namespace fact {

struct Tree {
  int N = 0;     // number of leaves (taxa)
  int cnt = 0;   // number of nodes
  int root = 0;  // root node id

  // Per-taxon arrays (size N + 5), indexed by taxon 1..N.
  std::vector<int> idx, label, minH, maxH, H;
  std::vector<char> good;            // good[t]: taxon t heads a kept split
  // Per-node arrays (size cnt + 5), indexed by node 0..cnt-1.
  std::vector<int> leaf, size, minL, maxL, parent, goodLabel;
  std::vector<std::vector<int>> G;   // children, directed away from the root

  Tree() = default;
  Tree(int N_, int cnt_);

  void precompute();                 // Day's perfect-hash relabelling
  std::string newick() const;        // integer-label Newick, no trailing ';'
};

// Build a Tree from an ape/TreeTools edge matrix that is in PREORDER (1-indexed
// node ids, parent listed before child, tips numbered 1..nTip).  nTip is the
// number of leaves.  The first row's parent is taken as the root (true for a
// preorder edge matrix).
Tree buildTreeFromEdge(const Rcpp::IntegerMatrix& edge, int nTip);

}  // namespace fact

#endif  // CONSTREE_FACT_TREE_H
