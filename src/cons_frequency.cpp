#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Frequency-difference consensus -- near-linear O(kn log n) algorithm of
// Jansson, Sung, Tabatabaee & Yang (2024, STACS, doi:10.4230/LIPIcs.STACS.2024.43),
// ported from their reference implementation freqdiff2.h (+ Tree, taxas_ranges,
// lca_preprocessing, radix_sort, utils) at github.com/tswddd2/FDCT_new (the
// software cited by that paper; used with permission).
//
// Three deliberate departures from the upstream single-shot batch code, each
// required to run repeatedly and portably from inside an R package (mirrors
// src/fact_tree.{h,cpp} / src/cons_greedy.cpp):
//   * NO file-scope globals.  freqdiff2.h threads ~30 raw scratch arrays through
//     file-scope `int* start; ...` globals and a static Tree taxa map; here all
//     scratch lives in a per-call FreqDiff context (std::vector members) and the
//     taxa count `n` is threaded explicitly, so the code is reentrant.
//   * NO raw owning new[]/VLAs.  All scratch is std::vector; owned sub-objects use
//     RAII (destructors / unique_ptr), fixing the upstream leaks (the ~30 arrays
//     freqdiff() never frees, prob_set's members, contracted subtrees, the
//     alloc_int_matrix block, radix_t's node pool).
//   * BOOST-FREE.  freqdiff2.h's only boost use is the `node_bitvec_t` struct,
//     which is DEAD CODE there (the live weight pass is calc_w_knlogn, a radix /
//     divide-and-conquer labelling -- the bitset weight calc lives only in the
//     OLD freqdiff.h).  We simply omit the dead struct; no bitset is needed.
//
// IMPORTANT (rooting): frequency-difference is an UNROOTED method, but the
// algorithm weights ROOTED clusters, so the unrooted consensus is recovered only
// when every input is rooted consistently.  The R wrapper Frequency() roots each
// input at taxon 1 (.FactEdges) before extracting edges, exactly the canonical
// unrooted->rooted reduction; do not call this on inconsistently rooted trees.

namespace {

// ===========================================================================
// utils.h
// ===========================================================================

inline int int_log2(int n) {
  int targetlevel = 0;
  while (n >>= 1) ++targetlevel;
  return targetlevel;
}

// Square int matrix as one block (m[0]) with row pointers, freed by
// free_int_matrix (upstream ~rmq_t leaked the block -- it delete[]'d only the
// row-pointer array, not m[0]).
inline int** alloc_int_matrix(int n, int v) {
  int** m = new int*[n];
  m[0] = new int[static_cast<size_t>(n) * n];
  for (int i = 1; i < n; i++) m[i] = m[0] + static_cast<size_t>(i) * n;
  std::fill(m[0], m[0] + static_cast<size_t>(n) * n, v);
  return m;
}
inline int** alloc_int_matrix(int n) { return alloc_int_matrix(n, 0); }
inline void free_int_matrix(int** m) {
  if (m) {
    delete[] m[0];
    delete[] m;
  }
}

// ===========================================================================
// Tree.h / Tree.cpp  (de-static-ified: no global taxa_ids/taxa_names; taxa are
// the 0-indexed integers supplied by the caller, recovered as taxa + 1 on
// output.  `n` -- the global taxon count, upstream Tree::get_taxas_num() -- is
// threaded explicitly by callers.)
// ===========================================================================

class Tree {
 public:
  static const int NONE = -1;

  class Node {
   public:
    std::vector<Node*> children;
    Node* parent;
    size_t pos_in_parent;
    int id, secondary_id;
    int taxa, weight, orig_w;
    int label, tree_id, spoil;
    size_t leaf_size, node_size;
    int depth;

    explicit Node(int id_)
        : parent(NULL), pos_in_parent(NONE), id(id_), secondary_id(id_),
          taxa(NONE), weight(0), orig_w(0), label(0), tree_id(0), spoil(0),
          leaf_size(0), node_size(0), depth(0) {}
    Node(int id_, int taxa_)
        : parent(NULL), pos_in_parent(NONE), id(id_), secondary_id(id_),
          taxa(taxa_), weight(0), orig_w(0), label(0), tree_id(0), spoil(0),
          leaf_size(0), node_size(0), depth(0) {}

    size_t get_children_num() { return children.size(); }
    bool is_leaf() { return taxa != NONE; }
    bool is_root() { return parent == NULL; }

    void add_child(Node* child) {
      child->parent = this;
      child->pos_in_parent = children.size();
      children.push_back(child);
    }
    void set_child(Node* child, size_t pos) {
      child->parent = this;
      child->pos_in_parent = pos;
      children[pos] = child;
    }
    void null_child(size_t pos) { children[pos] = NULL; }
    void clear_children() { children.clear(); }
    void fix_children() {  // drop NULL children, recompact, refresh pos_in_parent
      size_t curr_pos = 0;
      for (size_t i = 0; i < children.size(); i++) {
        if (children[i] != NULL) {
          children[curr_pos] = children[i];
          children[curr_pos]->pos_in_parent = curr_pos;
          curr_pos++;
        }
      }
      children.resize(curr_pos);
    }
  };

  explicit Tree(size_t nodes_num_hint = 0) : leaves_num(0) {
    if (nodes_num_hint > 0) nodes.reserve(nodes_num_hint);
  }
  // Deep copy (upstream Tree(Tree* other)).
  explicit Tree(Tree* other) : leaves_num(other->leaves_num) {
    nodes.reserve(other->get_nodes_num());
    for (Node* node : other->nodes) {
      Node* newnode = new Node(node->id, node->taxa);
      newnode->weight = node->weight;
      newnode->leaf_size = node->leaf_size;
      newnode->node_size = node->node_size;
      newnode->depth = node->depth;
      nodes.push_back(newnode);
      if (!node->is_root()) nodes[node->parent->id]->add_child(newnode);
      if (newnode->is_leaf()) taxa_to_leaf_map[newnode->taxa] = newnode;
    }
  }
  ~Tree() {
    for (Node* node : nodes) delete node;
  }
  Tree(const Tree&) = delete;
  Tree& operator=(const Tree&) = delete;

  Node* get_node(int i) { return nodes[i]; }
  Node* get_root() { return nodes[0]; }
  size_t get_nodes_num() { return nodes.size(); }
  Node* get_leaf(int taxa) { return taxa_to_leaf_map[taxa]; }
  size_t get_leaves_num() { return leaves_num; }

  Node* add_node(int taxa = -1) {
    Node* newnode = new Node(static_cast<int>(get_nodes_num()), taxa);
    nodes.push_back(newnode);
    if (taxa >= 0) {
      leaves_num++;
      taxa_to_leaf_map[taxa] = newnode;
    }
    return newnode;
  }

  void delete_nodes(bool* to_delete) {
    for (size_t i = 1; i < get_nodes_num(); i++) {
      if (to_delete[i]) {
        for (Node* child : nodes[i]->children) nodes[i]->parent->add_child(child);
        nodes[i]->parent->null_child(nodes[i]->pos_in_parent);
        delete nodes[i];
        nodes[i] = NULL;
      }
    }
    fix_tree();
  }

  // Re-root at `root`, renumber ids in DFS preorder (root id 0) and recompute
  // leaf_size / node_size / depth / pos_in_parent.  CRUCIAL: the whole algorithm
  // assumes preorder ids with each subtree a contiguous [id, id+node_size) range
  // and root == node 0; build-from-edges must call this to establish it.
  void fix_tree(Node* root = NULL) {
    if (root == NULL) root = get_root();
    nodes.clear();
    fix_tree_supp(root);
    nodes[0]->leaf_size = 0;
    nodes[0]->node_size = 0;
    nodes[0]->depth = 0;
    for (size_t i = 1; i < nodes.size(); i++) {
      nodes[i]->leaf_size = 0;
      nodes[i]->node_size = 0;
      nodes[i]->depth = nodes[i]->parent->depth + 1;
    }
    for (int i = static_cast<int>(nodes.size()) - 1; i > 0; i--) {
      if (nodes[i]->is_leaf()) nodes[i]->leaf_size = 1;
      nodes[i]->node_size += 1;
      nodes[i]->parent->leaf_size += nodes[i]->leaf_size;
      nodes[i]->parent->node_size += nodes[i]->node_size;
    }
  }

  void reorder() {  // make the heaviest subtree the first child of every node
    for (Node* node : nodes) {
      if (node->is_leaf()) continue;
      int heaviest = 0;
      for (size_t i = 1; i < node->get_children_num(); i++) {
        if (node->children[i]->leaf_size > node->children[heaviest]->leaf_size) {
          heaviest = static_cast<int>(i);
        }
      }
      Node* heaviest_node = node->children[heaviest];
      node->set_child(node->children[0], heaviest);
      node->set_child(heaviest_node, 0);
    }
  }

 private:
  std::vector<Node*> nodes;
  size_t leaves_num;
  std::unordered_map<int, Node*> taxa_to_leaf_map;

  void fix_tree_supp(Node* curr) {
    curr->id = static_cast<int>(nodes.size());
    nodes.push_back(curr);
    curr->fix_children();
    for (Node* child : curr->children) fix_tree_supp(child);
  }
};

// ===========================================================================
// taxas_ranges.h
// ===========================================================================

struct taxas_ranges_t {
  size_t taxas_num;
  std::vector<int> taxas;
  struct interval_t {
    int start, end;
  };
  std::vector<interval_t> intervals;

  taxas_ranges_t(size_t taxas_num_, size_t nodes_num)
      : taxas_num(0), taxas(taxas_num_), intervals(nodes_num) {}
};

inline void build_taxas_ranges_supp(Tree::Node* node, taxas_ranges_t* tr) {
  if (node->is_leaf()) {
    tr->taxas[tr->taxas_num] = node->taxa;
    tr->intervals[node->id].start = tr->intervals[node->id].end =
        static_cast<int>(tr->taxas_num);
    tr->taxas_num++;
    return;
  }
  for (Tree::Node* child : node->children) build_taxas_ranges_supp(child, tr);
  tr->intervals[node->id].start = tr->intervals[node->children[0]->id].start;
  tr->intervals[node->id].end =
      tr->intervals[node->children[node->get_children_num() - 1]->id].end;
}

inline taxas_ranges_t* build_taxas_ranges(Tree* tree) {
  taxas_ranges_t* tr =
      new taxas_ranges_t(tree->get_leaves_num(), tree->get_nodes_num());
  build_taxas_ranges_supp(tree->get_root(), tr);
  return tr;
}

// ranks[taxa] = left-to-right rank of that taxon's leaf; indexed by taxon id in
// [0, n), so sized with the global taxon count n (upstream get_taxas_num()).
inline int* get_taxas_ranks(taxas_ranges_t* tr, int n) {
  int* ranks = new int[n];
  for (size_t i = 0; i < tr->taxas_num; i++)
    ranks[tr->taxas[i]] = static_cast<int>(i);
  return ranks;
}

// ===========================================================================
// lca_preprocessing.h  (LCA via Euler tour + sparse-table/block RMQ; general
// RMQ via a Cartesian tree reduced to LCA)
// ===========================================================================

struct rmq_t {
  std::vector<std::vector<int> > M;
  std::vector<int> v;
  int block_size;
  std::vector<int> addresses;
  std::vector<int**> prep_blocks;

  rmq_t() : block_size(0) {}
  ~rmq_t() {
    for (int** block : prep_blocks) free_int_matrix(block);
  }
  rmq_t(const rmq_t&) = delete;
  rmq_t& operator=(const rmq_t&) = delete;
};

struct lca_t {
  std::vector<int> E, R;
  std::unique_ptr<rmq_t> rmq_prep;
  lca_t() : rmq_prep(new rmq_t) {}
};

struct gen_rmq_t {
  std::unique_ptr<lca_t> lca_prep;
  std::vector<int> v, pos_to_id, id_to_pos;
  gen_rmq_t() {}
};

inline void resize_to_logmul(std::vector<int>& v) {
  int log_size = std::max(int_log2(static_cast<int>(v.size())), 1);
  if (v.size() % log_size != 0)
    v.resize((v.size() / log_size + 1) * log_size, v[v.size() - 1] + 1);
}

inline void rmq_preprocess(rmq_t* rmq_prep, std::vector<int>& v) {
  resize_to_logmul(v);

  size_t size = v.size();
  int block_size = std::max(int_log2(static_cast<int>(size)), 1);
  int blocks = static_cast<int>(size) / block_size + (size % block_size != 0);

  std::vector<int> Ap(blocks, INT32_MAX);
  std::vector<int> B(blocks, 0);
  for (int i = 0; i < blocks; i++) {
    for (int j = 0; j < block_size; j++) {
      if (Ap[i] > v[i * block_size + j]) {
        Ap[i] = v[i * block_size + j];
        B[i] = i * block_size + j;
      }
    }
  }

  rmq_prep->M.push_back(std::vector<int>(Ap.size(), 0));
  for (size_t j = 0; j < Ap.size(); j++) rmq_prep->M[0][j] = B[j];
  for (int i = 1; (1 << i) <= static_cast<int>(Ap.size()); i++) {
    rmq_prep->M.push_back(std::vector<int>(Ap.size() - (1 << i) + 1, 0));
    for (size_t j = 0; j < Ap.size() - (1 << i) + 1; j++) {
      if (v[rmq_prep->M[i - 1][j]] <= v[rmq_prep->M[i - 1][j + (1 << (i - 1))]]) {
        rmq_prep->M[i][j] = rmq_prep->M[i - 1][j];
      } else {
        rmq_prep->M[i][j] = rmq_prep->M[i - 1][j + (1 << (i - 1))];
      }
    }
  }

  rmq_prep->prep_blocks.resize(size, NULL);
  for (int i = 0; i < blocks; i++) {
    int address = 0;
    for (int j = 1; j < block_size; j++) {
      address <<= 1;
      address |= (v[i * block_size + j] > v[i * block_size + j - 1]);
    }
    rmq_prep->addresses.push_back(address);
    if (rmq_prep->prep_blocks[address] == NULL) {
      rmq_prep->prep_blocks[address] = alloc_int_matrix(block_size);
      for (int j = 0; j < block_size; j++) {
        rmq_prep->prep_blocks[address][j][j] = j;
        for (int kk = j + 1; kk < block_size; kk++) {
          rmq_prep->prep_blocks[address][j][kk] = rmq_prep->prep_blocks[address][j][kk - 1];
          if (v[i * block_size + rmq_prep->prep_blocks[address][j][kk - 1]] >
              v[i * block_size + kk]) {
            rmq_prep->prep_blocks[address][j][kk] = kk;
          }
        }
      }
    }
  }
  rmq_prep->v = v;
  rmq_prep->block_size = block_size;
}

inline int rmq2(rmq_t* rmq_prep, int a, int b) {
  int ba = (a / rmq_prep->block_size) + (a % rmq_prep->block_size != 0);
  int bb = (b / rmq_prep->block_size) -
           (b % rmq_prep->block_size != rmq_prep->block_size - 1);
  int range = bb - ba + 1;
  if (range == 0) return a;
  int k = 0;
  while (range >>= 1) ++k;
  if (rmq_prep->v[rmq_prep->M[k][ba]] <= rmq_prep->v[rmq_prep->M[k][bb - (1 << k) + 1]]) {
    return rmq_prep->M[k][ba];
  } else {
    return rmq_prep->M[k][bb - (1 << k) + 1];
  }
}

inline int rmq(rmq_t* rmq_prep, int a, int b) {
  if (a / rmq_prep->block_size == b / rmq_prep->block_size) {
    int block_idx = a / rmq_prep->block_size;
    return block_idx * rmq_prep->block_size +
           rmq_prep->prep_blocks[rmq_prep->addresses[a / rmq_prep->block_size]]
                                [a % rmq_prep->block_size][b % rmq_prep->block_size];
  }
  int min_pos = rmq2(rmq_prep, a, b);
  if (a % rmq_prep->block_size != 0) {
    int temp_min = rmq_prep->prep_blocks[rmq_prep->addresses[a / rmq_prep->block_size]]
                                        [a % rmq_prep->block_size][rmq_prep->block_size - 1];
    if (rmq_prep->v[min_pos] >
        rmq_prep->v[(a / rmq_prep->block_size) * rmq_prep->block_size + temp_min]) {
      min_pos = (a / rmq_prep->block_size) * rmq_prep->block_size + temp_min;
    }
  }
  if (b % rmq_prep->block_size != rmq_prep->block_size - 1) {
    int temp_min = rmq_prep->prep_blocks[rmq_prep->addresses[b / rmq_prep->block_size]]
                                        [0][b % rmq_prep->block_size];
    if (rmq_prep->v[min_pos] >
        rmq_prep->v[(b / rmq_prep->block_size) * rmq_prep->block_size + temp_min]) {
      min_pos = (b / rmq_prep->block_size) * rmq_prep->block_size + temp_min;
    }
  }
  return min_pos;
}

inline int lca(lca_t* lca_prep, int u, int v) {
  if (lca_prep->R[u] < lca_prep->R[v]) {
    return lca_prep->E[rmq(lca_prep->rmq_prep.get(), lca_prep->R[u], lca_prep->R[v])];
  } else {
    return lca_prep->E[rmq(lca_prep->rmq_prep.get(), lca_prep->R[v], lca_prep->R[u])];
  }
}

inline void eulerian_walk(Tree::Node* node, std::vector<int>& E, std::vector<int>& L,
                          std::vector<int>& R, int depth) {
  E.push_back(node->id);
  L.push_back(depth);
  if (R[node->id] == -1) R[node->id] = static_cast<int>(E.size()) - 1;
  for (Tree::Node* child : node->children) {
    eulerian_walk(child, E, L, R, depth + 1);
    E.push_back(node->id);
    L.push_back(depth);
  }
}

inline lca_t* lca_preprocess(Tree* t) {
  lca_t* lca_prep = new lca_t;
  lca_prep->R.resize(t->get_nodes_num(), -1);
  eulerian_walk(t->get_root(), lca_prep->E, lca_prep->rmq_prep->v, lca_prep->R, 0);
  resize_to_logmul(lca_prep->E);
  resize_to_logmul(lca_prep->rmq_prep->v);
  rmq_preprocess(lca_prep->rmq_prep.get(), lca_prep->rmq_prep->v);
  return lca_prep;
}

inline int general_rmq(gen_rmq_t* gen_rmq_prep, int a, int b) {
  return gen_rmq_prep->id_to_pos[lca(gen_rmq_prep->lca_prep.get(),
                                     gen_rmq_prep->pos_to_id[a],
                                     gen_rmq_prep->pos_to_id[b])];
}

inline void general_rmq_preprocess(gen_rmq_t* gen_rmq_prep) {
  Tree* cartesian = new Tree;
  std::vector<Tree::Node*> orig_pos;
  Tree::Node* start = cartesian->add_node();
  Tree::Node* root = start;
  start->weight = gen_rmq_prep->v[0];
  orig_pos.push_back(start);
  for (size_t i = 1; i < gen_rmq_prep->v.size(); i++) {
    if (gen_rmq_prep->v[i] >= start->weight) {
      Tree::Node* node = cartesian->add_node();
      start->add_child(node);
      start = node;
    } else {
      while (start != NULL && start->weight > gen_rmq_prep->v[i]) start = start->parent;
      if (start == NULL) {
        start = cartesian->add_node();
        start->add_child(root);
        root = start;
      } else {
        Tree::Node* node = cartesian->add_node();
        node->add_child(start->children[start->children.size() - 1]);
        start->set_child(node, start->children.size() - 1);
        start = node;
      }
    }
    start->weight = gen_rmq_prep->v[i];
    orig_pos.push_back(start);
  }
  cartesian->fix_tree(root);
  gen_rmq_prep->id_to_pos.resize(orig_pos.size());
  for (Tree::Node* node : orig_pos) {
    gen_rmq_prep->id_to_pos[node->id] = static_cast<int>(gen_rmq_prep->pos_to_id.size());
    gen_rmq_prep->pos_to_id.push_back(node->id);
  }
  gen_rmq_prep->lca_prep.reset(lca_preprocess(cartesian));
  delete cartesian;
}

// ===========================================================================
// radix_sort.h  (stable O(n+k) counting/radix sort; own a node pool with RAII)
// ===========================================================================

struct radix_node_t {
  int key;
  void* val;
  radix_node_t* prev;
  void assign(int _key, void* _val) {
    key = _key;
    val = _val;
  }
};

struct radix_t {
  int node_num, key_num;
  int n, k;
  std::vector<radix_node_t*> tails;
  std::vector<radix_node_t> nodes_storage, out_storage;
  std::vector<radix_node_t*> nodes, out;

  radix_t(int node_num_, int key_num_)
      : node_num(node_num_), key_num(key_num_), n(0), k(0),
        tails(key_num_ + 1, NULL),
        nodes_storage(node_num_), out_storage(node_num_),
        nodes(node_num_), out(node_num_) {
    for (int i = 0; i < node_num; i++) {
      nodes[i] = &nodes_storage[i];
      out[i] = &out_storage[i];
    }
  }
  radix_t(const radix_t&) = delete;
  radix_t& operator=(const radix_t&) = delete;

  void clear() { n = k = 0; }

  void add(int key, void* val) {
    if (n == node_num) Rcpp::stop("frequencyConsensusCpp: radix sort n overflow.");
    if (key > k) k = key;
    nodes[n++]->assign(key, val);
    if (k > key_num) Rcpp::stop("frequencyConsensusCpp: radix sort k overflow.");
  }

  void sort(bool desc = false) {
    int i, j, key;
    for (i = 0; i <= k; i++) tails[i] = NULL;
    for (i = 0; i < n; i++) nodes[i]->prev = NULL;
    for (i = 0; i < n; i++) {
      key = nodes[i]->key;
      if (tails[key] == NULL) {
        tails[key] = nodes[i];
      } else {
        nodes[i]->prev = tails[key];
        tails[key] = nodes[i];
      }
    }
    j = n - 1;
    radix_node_t* node;
    for (i = 0; i <= k; i++) {
      node = desc ? tails[i] : tails[k - i];
      while (node != NULL) {
        *out[j] = *node;
        node = node->prev;
        j--;
      }
    }
  }

  void quicksort(bool desc = false) {
    for (int i = 0; i < n; i++) *out[i] = *nodes[i];
    if (desc) {
      std::sort(out.begin(), out.begin() + n,
                [](radix_node_t* x, radix_node_t* y) { return x->key > y->key; });
    } else {
      std::sort(out.begin(), out.begin() + n,
                [](radix_node_t* x, radix_node_t* y) { return x->key < y->key; });
    }
  }
};

// ===========================================================================
// freqdiff2.h: weight computation -- calc_w_knlogn (O(kn log n))
// ===========================================================================

struct label_node {
  int id;
  label_node *prev_node, *father;
  int label, lval, rval;
};

struct label_prob_set {
  std::vector<Tree*> trees;
  std::vector<std::unique_ptr<lca_t> > lcas;
  int n, k;
  std::vector<label_node*> vis;
  std::vector<std::vector<label_node> > inp;
  std::unique_ptr<radix_t> radix;

  explicit label_prob_set(std::vector<Tree*>& _trees) {
    k = static_cast<int>(_trees.size());
    n = static_cast<int>(_trees[0]->get_leaves_num());
    trees.resize(k);
    lcas.resize(k);
    inp.resize(k);
    for (int i = 0; i < k; i++) {
      trees[i] = _trees[i];
      inp[i].reserve(2 * n);
      for (int j = 0; j < static_cast<int>(trees[i]->get_nodes_num()); j++)
        inp[i].push_back({j, NULL, NULL, 0, 0, 0});
      for (int j = 0; j < static_cast<int>(trees[i]->get_nodes_num()); j++) {
        Tree::Node* node = trees[i]->get_node(j);
        if (node->parent == NULL)
          inp[i][j].father = NULL;
        else
          inp[i][j].father = &inp[i][node->parent->id];
      }
      lcas[i].reset(lca_preprocess(trees[i]));
    }
    vis.assign(2 * n, NULL);
    radix.reset(new radix_t(2 * n * k, 2 * n * k));
  }
};

// build the left/right label-node sublists for taxa range [st, ed)
inline std::vector<label_node>* build_sublist(label_prob_set* prob,
                                              std::vector<label_node>* lnodes,
                                              int st, int ed) {
  Tree::Node *prev, *node, *lca_node;
  std::vector<label_node*>& vis = prob->vis;
  std::vector<label_node>* res;
  int n = ed - st, k = prob->k;
  int i, j;
  res = new std::vector<label_node>[k];
  for (i = 0; i < k; i++) {
    res[i].reserve(2 * n);
    for (j = 0; j < static_cast<int>(lnodes[i].size()); j++)
      vis[lnodes[i][j].id] = &lnodes[i][j];
    prev = NULL;
    for (j = 0; j < static_cast<int>(lnodes[i].size()); j++) {
      node = prob->trees[i]->get_node(lnodes[i][j].id);
      if (!node->is_leaf()) continue;
      if (!(node->taxa >= st && node->taxa < ed)) continue;
      res[i].push_back({node->id, vis[node->id], NULL, 0, 0, 0});
      if (prev != NULL) {
        lca_node = prob->trees[i]->get_node(lca(prob->lcas[i].get(), node->id, prev->id));
        if (vis[lca_node->id] != NULL) {
          res[i].push_back({lca_node->id, vis[lca_node->id], NULL, 0, 0, 0});
          vis[lca_node->id] = NULL;
        }
      }
      prev = node;
    }
    label_node* plnode;
    for (j = 0; j < static_cast<int>(lnodes[i].size()); j++) vis[lnodes[i][j].id] = NULL;
    for (j = 0; j < static_cast<int>(res[i].size()); j++) vis[res[i][j].id] = &res[i][j];
    for (j = 0; j < static_cast<int>(res[i].size()); j++) {
      plnode = res[i][j].prev_node->father;
      while (plnode != NULL && vis[plnode->id] == NULL) plnode = plnode->father;
      if (plnode == NULL)
        res[i][j].father = NULL;
      else
        res[i][j].father = vis[plnode->id];
    }
  }
  return res;
}

inline void label_nodes(label_prob_set* prob, std::vector<label_node>* lnodes,
                        int t_start, int t_end) {
  int k = prob->k;
  int i, j;
  if (t_end - t_start == 1) {
    for (i = 0; i < k; i++) lnodes[i][0].label = 1;
    return;
  }
  int mid = (t_start + t_end) / 2;
  std::vector<label_node>* llnodes = build_sublist(prob, lnodes, t_start, mid);
  std::vector<label_node>* rlnodes = build_sublist(prob, lnodes, mid, t_end);
  label_nodes(prob, llnodes, t_start, mid);
  label_nodes(prob, rlnodes, mid, t_end);

  label_node* node;
  for (i = 0; i < k; i++) {
    for (j = 0; j < static_cast<int>(llnodes[i].size()); j++)
      llnodes[i][j].prev_node->lval = llnodes[i][j].label;
    for (j = 0; j < static_cast<int>(rlnodes[i].size()); j++)
      rlnodes[i][j].prev_node->rval = rlnodes[i][j].label;
    for (j = 0; j < static_cast<int>(lnodes[i].size()); j++) {
      node = &lnodes[i][j];
      if (node->lval != 0)
        while (node->father != NULL && node->father->lval == 0) {
          node->father->lval = node->lval;
          node = node->father;
        }
      node = &lnodes[i][j];
      if (node->rval != 0)
        while (node->father != NULL && node->father->rval == 0) {
          node->father->rval = node->rval;
          node = node->father;
        }
    }
  }

  radix_t* radix = prob->radix.get();
  radix->clear();
  for (i = 0; i < k; i++)
    for (j = 0; j < static_cast<int>(lnodes[i].size()); j++)
      radix->add(lnodes[i][j].lval, &lnodes[i][j]);
  radix->sort();

  int radix_n = radix->n;
  radix->clear();
  for (i = 0; i < radix_n; i++) {
    node = (label_node*)radix->out[i]->val;
    radix->add(node->rval, node);
  }
  radix->sort();

  label_node* prev = NULL;
  int cnt = 0;
  for (i = 0; i < radix->n; i++) {
    node = (label_node*)radix->out[i]->val;
    if (prev == NULL || prev->lval != node->lval || prev->rval != node->rval) cnt++;
    node->label = cnt;
    prev = node;
  }

  delete[] llnodes;
  delete[] rlnodes;
}

// node->weight = number of input trees whose cluster set contains node's cluster
inline void calc_w_knlogn(std::vector<Tree*>& trees, int n) {
  size_t k = trees.size();
  std::vector<int> count(2 * k * n, 0);
  std::unique_ptr<label_prob_set> prob(new label_prob_set(trees));
  label_nodes(prob.get(), prob->inp.data(), 0, n);
  for (size_t i = 0; i < trees.size(); i++)
    for (int j = 0; j < static_cast<int>(trees[i]->get_nodes_num()); j++)
      count[prob->inp[i][j].label]++;
  for (size_t i = 0; i < trees.size(); i++)
    for (int j = 0; j < static_cast<int>(trees[i]->get_nodes_num()); j++) {
      Tree::Node* node = trees[i]->get_node(j);
      node->weight = count[prob->inp[i][j].label];
    }
}

// ===========================================================================
// freqdiff2.h: subpath max-weight queries on a (centroid-decomposed) tree
// ===========================================================================

struct subpath_query_info_t {
  size_t nodes_num;
  std::vector<Tree::Node*> cp_roots;
  std::vector<std::unique_ptr<gen_rmq_t> > cp_rmqs;
  std::vector<int> depths;
  explicit subpath_query_info_t(size_t nodes_num_)
      : nodes_num(nodes_num_), cp_roots(nodes_num_, NULL),
        cp_rmqs(nodes_num_), depths(nodes_num_, 0) {}
};

inline subpath_query_info_t* preprocess_subpaths_queries(Tree* tree, int n) {
  subpath_query_info_t* info = new subpath_query_info_t(static_cast<size_t>(n) * 2);
  std::vector<Tree::Node*>& cp_roots = info->cp_roots;
  cp_roots[0] = tree->get_root();
  for (size_t i = 1; i < tree->get_nodes_num(); i++) {
    Tree::Node* node = tree->get_node(static_cast<int>(i));
    if (node->pos_in_parent == 0)
      cp_roots[i] = cp_roots[node->parent->id];
    else
      cp_roots[i] = node;
  }
  for (size_t i = 0; i < tree->get_nodes_num(); i++) {
    if (info->cp_rmqs[cp_roots[i]->id] == NULL)
      info->cp_rmqs[cp_roots[i]->id].reset(new gen_rmq_t);
    info->cp_rmqs[cp_roots[i]->id]->v.push_back(-tree->get_node(static_cast<int>(i))->weight);
  }
  for (size_t i = 0; i < tree->get_nodes_num(); i++) {
    if (info->cp_rmqs[i] != NULL && info->cp_rmqs[i]->v.size() > 1)
      general_rmq_preprocess(info->cp_rmqs[i].get());
  }
  std::vector<int>& depths = info->depths;
  depths[0] = 0;
  for (size_t i = 1; i < tree->get_nodes_num(); i++)
    depths[i] = 1 + depths[tree->get_node(static_cast<int>(i))->parent->id];
  return info;
}

inline int max_subpath_query(subpath_query_info_t* subpq_info, Tree::Node* ancestor,
                             Tree::Node* descendant) {
  Tree::Node* curr = descendant;
  int res = 0;
  while (subpq_info->cp_roots[curr->id]->id != subpq_info->cp_roots[ancestor->id]->id) {
    gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->id]->id].get();
    int query_endp =
        subpq_info->depths[curr->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
    if (query_endp == 0)
      res = std::min(res, currpath_rmq->v[0]);
    else
      res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, 0, query_endp)]);
    curr = subpq_info->cp_roots[curr->id]->parent;
  }
  gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->id]->id].get();
  int query_startp =
      subpq_info->depths[ancestor->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
  int query_endp =
      subpq_info->depths[curr->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
  if (query_startp < query_endp)
    res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, query_startp + 1, query_endp)]);
  return -res;
}

struct prob_set {
  Tree *tree1, *tree2;
  taxas_ranges_t *t1_tr, *t2_tr;
  lca_t* t2_lcas;
  subpath_query_info_t* subpq_t2;
  prob_set(Tree* t1, Tree* t2, int n)
      : tree1(t1), tree2(t2), t1_tr(build_taxas_ranges(t1)),
        t2_tr(build_taxas_ranges(t2)), t2_lcas(lca_preprocess(t2)),
        subpq_t2(preprocess_subpaths_queries(t2, n)) {}
  ~prob_set() {
    delete t1_tr;
    delete t2_tr;
    delete t2_lcas;
    delete subpq_t2;
  }
  prob_set(const prob_set&) = delete;
  prob_set& operator=(const prob_set&) = delete;
};

inline int unionset_find(int* ancest, int id) {
  int res = id, next;
  while (ancest[res] != res) res = ancest[res];
  while (id != res) {
    next = ancest[id];
    ancest[id] = res;
    id = next;
  }
  return res;
}

struct Interval {
  int l, r;
};

// ===========================================================================
// freqdiff2.h: the FreqDiff context -- all the per-call scratch upstream kept in
// file-scope globals, plus the cluster-filtering / tree-merging core.
// ===========================================================================

class FreqDiff {
 public:
  FreqDiff(int n_, int k_) : n(n_), k(k_) {
    start.assign(2 * n, 0);
    stop.assign(2 * n, 0);
    e.assign(n, 0);
    m.assign(2 * n, 0);
    rsort_lists.assign(n, std::vector<Tree::Node*>());
    _left.assign(n, NULL);
    _right.assign(n, NULL);
    orig_pos_in_parent.assign(2 * n, 0);
    leaf_p_index.assign(n, 0);
    vleft.assign(2 * n, 0);
    vright.assign(2 * n, 0);
    pointer.assign(2 * n, 0);
    levels.assign(2 * n, 0);
    ids.assign(2 * n, 0);
    parent.assign(2 * n, 0);
    exists.assign(2 * n, 0);
    tree_nodes.assign(2 * n, NULL);
    minv.assign(3 * n, 0);
    maxv.assign(3 * n, 0);
    fillv.assign(3 * n, 0);
    ancest.assign(n, 0);
    inext.assign(2 * n, 0);
    iid.assign(2 * n, 0);
    str_n.assign(2 * n, 0);
    intervals.assign(3 * n, Interval());
    origw_to_w.assign(k + 1, 0);
    subtree_cnt.assign(2 * n, 0);
    radix.reset(new radix_t(5 * n, k));
  }

  // Consume `trees` (deletes them) and return the freq-diff consensus tree.
  Tree* run(std::vector<Tree*>& trees);

 private:
  int n, k;
  std::vector<int> start, stop, e, m;
  std::vector<std::vector<Tree::Node*> > rsort_lists;
  std::vector<Tree::Node*> _left, _right;
  std::vector<size_t> orig_pos_in_parent;
  std::vector<int> leaf_p_index;
  std::vector<int> vleft, vright, pointer, levels, ids, parent;
  std::vector<char> exists;
  std::vector<Tree::Node*> tree_nodes;
  std::vector<int> minv, maxv, fillv, ancest, inext, iid, str_n;
  std::vector<Interval> intervals;
  std::vector<int> origw_to_w, subtree_cnt;
  std::unique_ptr<radix_t> radix;

  void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks) {
    for (int i = static_cast<int>(tree1->get_nodes_num()) - 1; i >= 0; i--) {
      Tree::Node* node = tree1->get_node(i);
      if (node->is_leaf()) {
        start[i] = stop[i] = t2_leaves_ranks[node->taxa];
      } else {
        start[i] = INT32_MAX;
        stop[i] = 0;
        for (Tree::Node* child : node->children) {
          if (start[i] > start[child->id]) start[i] = start[child->id];
          if (stop[i] < stop[child->id]) stop[i] = stop[child->id];
        }
      }
    }
  }

  void compute_m(Tree::Node* node) {
    if (node->is_leaf()) m[node->id] = e[node->taxa];
    for (Tree::Node* child : node->children) {
      compute_m(child);
      if (m[node->id] > m[child->id]) m[node->id] = m[child->id];
    }
    if (!node->is_root()) rsort_lists[m[node->id]].push_back(node);
  }

  Tree* contract_tree_fast(Tree* tree, lca_t* lcas, std::vector<int>& marked,
                           subpath_query_info_t* subpq);
  void filter_clusters_nlogn(prob_set* prob, Tree::Node* t1_root, Tree* tree2, bool* to_del);
  void filter(Tree* tree1, Tree* tree2, bool* to_del);
  void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas);
};

// Contract `tree` down to the leaves in `marked` (a left-to-right ordered taxa
// subset), splicing in special weight nodes for the subpaths between kept nodes.
Tree* FreqDiff::contract_tree_fast(Tree* tree, lca_t* lcas, std::vector<int>& marked,
                                   subpath_query_info_t* subpq) {
  if (marked.empty()) return NULL;

  int count = 0;
  levels[count] = tree->get_leaf(marked[0])->depth;
  ids[count] = tree->get_leaf(marked[0])->id;
  count++;
  for (size_t i = 1; i < marked.size(); i++) {
    int lca_id = lca(lcas, tree->get_leaf(marked[i - 1])->id, tree->get_leaf(marked[i])->id);
    levels[count] = tree->get_node(lca_id)->depth;
    ids[count] = tree->get_node(lca_id)->id;
    count++;
    levels[count] = tree->get_leaf(marked[i])->depth;
    ids[count] = tree->get_leaf(marked[i])->id;
    count++;
  }

  std::stack<int> Sl;
  for (int i = 0; i < count; i++) {
    while (!Sl.empty() && levels[i] <= levels[Sl.top()]) Sl.pop();
    vleft[i] = Sl.empty() ? -1 : Sl.top();
    Sl.push(i);
  }
  std::stack<int> Sr;
  for (int i = count - 1; i >= 0; i--) {
    while (!Sr.empty() && levels[i] <= levels[Sr.top()]) Sr.pop();
    vright[i] = Sr.empty() ? -1 : Sr.top();
    Sr.push(i);
  }

  int root_pos = -1;
  std::fill(exists.begin(), exists.begin() + count, 1);
  for (int i = 0; i < count; i++) pointer[i] = i;
  for (int i = 0; i < count; i++) {
    parent[i] = -1;
    if (vleft[i] == -1 && vright[i] == -1) {
      if (root_pos == -1) root_pos = i;
    } else if (vleft[i] == -1) {
      parent[i] = vright[i];
    } else if (vright[i] == -1) {
      parent[i] = pointer[vleft[i]];
    } else {
      if (levels[vleft[i]] >= levels[vright[i]]) {
        parent[i] = pointer[vleft[i]];
        if (levels[vleft[i]] == levels[vright[i]]) {
          pointer[vright[i]] = pointer[vleft[i]];
          exists[vright[i]] = 0;
        }
      } else {
        parent[i] = pointer[vright[i]];
      }
    }
  }

  Tree* new_tree = new Tree(count * 2);
  std::fill(tree_nodes.begin(), tree_nodes.begin() + count, (Tree::Node*)NULL);
  for (int i = 0; i < count; i++) {
    if (i % 2 == 0 || exists[i]) {
      if (tree_nodes[i] == NULL) {
        if (i % 2 == 0)
          tree_nodes[i] = new_tree->add_node(marked[i / 2]);
        else
          tree_nodes[i] = new_tree->add_node();
      }
      tree_nodes[i]->weight = tree->get_node(ids[i])->weight;
      tree_nodes[i]->orig_w = tree->get_node(ids[i])->orig_w;
      tree_nodes[i]->secondary_id = ids[i];
    }
    if (tree_nodes[i] != NULL && parent[i] != -1) {
      if (tree_nodes[parent[i]] == NULL) tree_nodes[parent[i]] = new_tree->add_node();
      tree_nodes[parent[i]]->add_child(tree_nodes[i]);
    }
  }

  if (subpq != NULL) {
    size_t newtree_nodes = new_tree->get_nodes_num();
    for (size_t i = 0; i < newtree_nodes; i++) {
      Tree::Node* curr_node = new_tree->get_node(static_cast<int>(i));
      if (curr_node->is_root()) continue;
      Tree::Node* desc_par = tree->get_node(curr_node->secondary_id)->parent;
      Tree::Node* anc_par = tree->get_node(curr_node->parent->secondary_id);
      int msq = max_subpath_query(subpq, anc_par, desc_par);
      if (msq > 0) {
        Tree::Node* sp_node = new_tree->add_node();
        sp_node->weight = sp_node->orig_w = msq;
        sp_node->secondary_id = anc_par->id;
        curr_node->parent->null_child(curr_node->pos_in_parent);
        curr_node->parent->add_child(sp_node);
        sp_node->add_child(curr_node);
      }
    }
  }

  new_tree->fix_tree(tree_nodes[root_pos]);
  for (int i = 0; i < static_cast<int>(new_tree->get_nodes_num()); i++) {
    Tree::Node* node = new_tree->get_node(i);
    Tree::Node* fnode = tree->get_node(node->secondary_id);
    if (fnode->spoil || fnode->leaf_size > node->leaf_size) node->spoil = 1;
  }
  return new_tree;
}

// Mark for deletion every cluster of tree1 strictly out-weighted by a conflicting
// cluster of tree2 (the heart of the freq-diff filter; recursion on the centroid
// path with tree contraction + a max-Manhattan-skyline sweep).
void FreqDiff::filter_clusters_nlogn(prob_set* prob, Tree::Node* t1_root, Tree* tree2,
                                     bool* to_del) {
  Tree* tree1 = prob->tree1;
  taxas_ranges_t* t1_tr = prob->t1_tr;
  int taxa_st = t1_tr->intervals[t1_root->id].start;
  int t1_st = t1_root->id;
  int t1_ed = t1_st + static_cast<int>(t1_root->node_size);

  std::vector<Tree::Node*> path;
  Tree::Node* node = t1_root;
  while (!node->is_leaf()) node = node->children[0];
  while (true) {
    path.push_back(node);
    if (node == t1_root) break;
    node = node->parent;
  }

  int nn = static_cast<int>(tree2->get_leaves_num());
  int mm = 0, pn = static_cast<int>(path.size());
  std::vector<Tree::Node*> s1_roots;

  leaf_p_index[path[0]->taxa] = -1;
  Tree::Node* pnode;
  int pi, i, j, nid;
  for (pi = 0; pi < pn; pi++) {
    pnode = path[pi];
    pnode->tree_id = -1;
    for (i = 1; i < static_cast<int>(pnode->children.size()); i++) {
      node = pnode->children[i];
      nid = node->id;
      for (j = t1_tr->intervals[nid].start; j <= t1_tr->intervals[nid].end; j++)
        leaf_p_index[t1_tr->taxas[j]] = mm;
      mm++;
      s1_roots.push_back(node);
    }
  }

  std::vector<int> layer_t1_w(t1_ed - t1_st);
  for (i = t1_st; i < t1_ed; i++) layer_t1_w[i - t1_st] = tree1->get_node(i)->weight;
  for (i = 0; i < static_cast<int>(tree2->get_nodes_num()); i++)
    origw_to_w[tree2->get_node(i)->orig_w] = tree2->get_node(i)->weight;

  std::vector<Tree*> sub_t2(mm, NULL);
  std::vector<std::vector<int> > marked(mm);
  taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
  for (i = 0; i < nn; i++)
    if (leaf_p_index[t2_tr->taxas[i]] > -1)
      marked[leaf_p_index[t2_tr->taxas[i]]].push_back(t2_tr->taxas[i]);
  delete t2_tr;

  Tree::Node* snode;
  radix->clear();
  for (i = 0; i < mm; i++) {
    node = s1_roots[i];
    sub_t2[i] = contract_tree_fast(prob->tree2, prob->t2_lcas, marked[i], prob->subpq_t2);
    nid = node->id;
    for (j = nid; j < nid + static_cast<int>(node->node_size); j++) {
      snode = tree1->get_node(j);
      snode->tree_id = i;
      radix->add(snode->weight, snode);
    }
    if (sub_t2[i] != NULL) {
      for (j = 0; j < static_cast<int>(sub_t2[i]->get_nodes_num()); j++) {
        snode = sub_t2[i]->get_node(j);
        snode->weight = origw_to_w[snode->orig_w];
        snode->tree_id = i;
        radix->add(snode->weight, snode);
      }
    }
  }
  radix->sort();

  std::fill(subtree_cnt.begin(), subtree_cnt.begin() + mm, 0);
  Tree::Node* prev = NULL;
  for (i = 0; i < radix->n; i++) {
    node = (Tree::Node*)radix->out[i]->val;
    if (prev == NULL || prev->tree_id != node->tree_id || prev->orig_w != node->orig_w)
      subtree_cnt[node->tree_id]++;
    node->weight = subtree_cnt[node->tree_id];
    prev = node;
  }

  for (i = 0; i < mm; i++) {
    if (sub_t2[i] != NULL) filter_clusters_nlogn(prob, s1_roots[i], sub_t2[i], to_del);
  }

  for (i = t1_st; i < t1_ed; i++) tree1->get_node(i)->weight = layer_t1_w[i - t1_st];

  for (i = 0; i < pn; i++) str_n[i] = t1_tr->intervals[path[i]->id].end - taxa_st;

  radix->clear();
  int sid, anc_id;
  for (i = 0; i < nn; i++) ancest[i] = i;
  for (i = static_cast<int>(tree2->get_nodes_num()) - 1; i >= 0; i--) {
    node = tree2->get_node(i);
    if (node->is_leaf()) {
      minv[i] = maxv[i] = t1_tr->intervals[tree1->get_leaf(node->taxa)->id].start - taxa_st;
      fillv[i] = -1;
    } else {
      for (j = 0; j < static_cast<int>(node->children.size()); j++) {
        sid = node->children[j]->id;
        if (j == 0) {
          minv[i] = minv[sid];
          maxv[i] = maxv[sid];
          fillv[i] = fillv[sid];
        } else {
          minv[i] = std::min(minv[i], minv[sid]);
          maxv[i] = std::max(maxv[i], maxv[sid]);
          fillv[i] = std::max(fillv[i], fillv[sid]);
        }
      }
      anc_id = unionset_find(ancest.data(), minv[node->children[0]->id]);
      for (j = 1; j < static_cast<int>(node->children.size()); j++) {
        sid = node->children[j]->id;
        ancest[unionset_find(ancest.data(), minv[sid])] = anc_id;
      }
    }
    anc_id = unionset_find(ancest.data(), minv[i]);
    while (fillv[i] < nn - 1 && unionset_find(ancest.data(), fillv[i] + 1) == anc_id) fillv[i]++;
    intervals[i].l = std::max(minv[i], fillv[i] + 1);
    intervals[i].r = (node->spoil ? nn - 1 : maxv[i] - 1);
    if (intervals[i].l <= intervals[i].r) radix->add(node->weight, &intervals[i]);
  }
  radix->sort(true);

  int ti = 0;
  for (i = 0; i < nn; i++) {
    inext[i] = str_n[ti];
    if (i == str_n[ti]) {
      iid[i] = path[ti]->id;
      ti++;
    } else {
      iid[i] = -1;
    }
  }

  Interval* interval;
  inext[nn] = nn;
  for (i = 0; i < radix->n; i++) {
    interval = (Interval*)radix->out[i]->val;
    j = unionset_find(inext.data(), interval->l);
    while (j <= interval->r) {
      if (iid[j] >= 0 && radix->out[i]->key >= tree1->get_node(iid[j])->weight)
        to_del[iid[j]] = true;
      inext[j] = unionset_find(inext.data(), j + 1);
      j = inext[j];
    }
  }

  for (i = 0; i < mm; i++) delete sub_t2[i];
}

void FreqDiff::filter(Tree* tree1, Tree* tree2, bool* to_del) {
  for (int i = 0; i < static_cast<int>(tree1->get_nodes_num()); i++)
    tree1->get_node(i)->orig_w = tree1->get_node(i)->weight;
  for (int i = 0; i < static_cast<int>(tree2->get_nodes_num()); i++)
    tree2->get_node(i)->orig_w = tree2->get_node(i)->weight;

  std::unique_ptr<prob_set> prob(new prob_set(tree1, tree2, n));

  // compress node weights in n log n (in case k >> n)
  if (radix->k > radix->n) {
    radix->clear();
    for (int i = 0; i < static_cast<int>(tree1->get_nodes_num()); i++)
      radix->add(tree1->get_node(i)->weight, tree1->get_node(i));
    for (int i = 0; i < static_cast<int>(tree2->get_nodes_num()); i++)
      radix->add(tree2->get_node(i)->weight, tree2->get_node(i));
    radix->quicksort();
    int cnt = 0;
    for (int i = 0; i < radix->n; i++) {
      Tree::Node* node = (Tree::Node*)radix->out[i]->val;
      if (i == 0 || radix->out[i - 1]->key != radix->out[i]->key) cnt++;
      node->weight = cnt;
    }
  }

  filter_clusters_nlogn(prob.get(), tree1->get_root(), tree2, to_del);

  for (int i = 0; i < static_cast<int>(tree1->get_nodes_num()); i++)
    tree1->get_node(i)->weight = tree1->get_node(i)->orig_w;
  for (int i = 0; i < static_cast<int>(tree2->get_nodes_num()); i++)
    tree2->get_node(i)->weight = tree2->get_node(i)->orig_w;
}

// Graft tree1's surviving clusters into tree2 (Section 2.4 of the paper).
void FreqDiff::merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas) {
  for (int jj = 0; jj < n; jj++) e[t1_tr->taxas[jj]] = jj;
  std::fill(m.begin(), m.begin() + tree2->get_nodes_num(), INT32_MAX);
  for (int i = 0; i < n; i++) rsort_lists[i].clear();

  compute_m(tree2->get_root());

  for (size_t i = 0; i < tree2->get_nodes_num(); i++)
    tree2->get_node(static_cast<int>(i))->clear_children();
  for (int i = 0; i < n; i++)
    for (Tree::Node* it : rsort_lists[i]) it->parent->add_child(it);

  taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
  int* t2_ranks = get_taxas_ranks(t2_tr, n);
  compute_start_stop(tree1, tree2, t2_ranks);
  delete[] t2_ranks;

  for (int i = 0; i < n; i++) {
    Tree::Node* curr = tree2->get_leaf(i);
    Tree::Node* par = curr->parent;
    while (par != NULL && *(par->children.begin()) == curr) {
      curr = par;
      par = curr->parent;
    }
    _left[i] = curr;
    curr = tree2->get_leaf(i);
    par = curr->parent;
    while (par != NULL && *(par->children.rbegin()) == curr) {
      curr = par;
      par = curr->parent;
    }
    _right[i] = curr;
  }

  for (size_t i = 0; i < tree2->get_nodes_num(); i++)
    orig_pos_in_parent[i] = tree2->get_node(static_cast<int>(i))->pos_in_parent;

  for (int i = static_cast<int>(tree1->get_nodes_num()) - 1; i >= 1; i--) {
    Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
    Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
    if (a == b) continue;
    Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->id, b->id));
    Tree::Node* a_left = _left[a->taxa];
    Tree::Node* b_right = _right[b->taxa];
    size_t du_pos = (a_left->depth > ru->depth) ? orig_pos_in_parent[a_left->id] : 0;
    size_t eu_pos = (b_right->depth > ru->depth) ? orig_pos_in_parent[b_right->id]
                                                 : ru->get_children_num() - 1;
    if (du_pos == 0 && eu_pos == ru->get_children_num() - 1) continue;
    Tree::Node* newnode = tree2->add_node();
    newnode->weight = tree1->get_node(i)->weight;
    newnode->leaf_size = tree1->get_node(i)->leaf_size;
    for (size_t jj = du_pos; jj <= eu_pos; jj++) {
      if (ru->children[jj] != NULL) {
        newnode->add_child(ru->children[jj]);
        ru->null_child(jj);
      }
    }
    ru->set_child(newnode, du_pos);
  }
  tree2->fix_tree();
  delete t2_tr;
}

Tree* FreqDiff::run(std::vector<Tree*>& trees) {
  calc_w_knlogn(trees, n);

  for (size_t i = 0; i < trees.size(); i++) trees[i]->reorder();

  std::vector<char> to_del_t(2 * n), to_del_ti(2 * n);

  Tree* T = new Tree(trees[0]);
  for (size_t i = 1; i < trees.size(); i++) {
    Tree* Ti = new Tree(trees[i]);

    std::fill(to_del_ti.begin(), to_del_ti.begin() + Ti->get_nodes_num(), 0);
    filter(Ti, T, reinterpret_cast<bool*>(to_del_ti.data()));

    std::fill(to_del_t.begin(), to_del_t.begin() + T->get_nodes_num(), 0);
    filter(T, Ti, reinterpret_cast<bool*>(to_del_t.data()));

    Ti->delete_nodes(reinterpret_cast<bool*>(to_del_ti.data()));
    T->delete_nodes(reinterpret_cast<bool*>(to_del_t.data()));

    lca_t* lca_T = lca_preprocess(T);
    taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

    merge_trees(Ti, T, tr_Ti, lca_T);

    delete lca_T;
    delete tr_Ti;

    for (size_t j = 0; j < T->get_nodes_num(); j++)
      T->get_node(static_cast<int>(j))->secondary_id = T->get_node(static_cast<int>(j))->id;

    delete Ti;
  }

  std::fill(to_del_t.begin(), to_del_t.begin() + T->get_nodes_num(), 0);
  for (size_t i = 0; i < trees.size(); i++) {
    filter(T, trees[i], reinterpret_cast<bool*>(to_del_t.data()));
    delete trees[i];
  }
  T->delete_nodes(reinterpret_cast<bool*>(to_del_t.data()));

  return T;
}

// ===========================================================================
// Marshalling
// ===========================================================================

// Build an FDCT Tree from a PREORDER ape/TreeTools edge matrix (1-indexed node
// ids, parent before child, tips 1..nTip, rooted at taxon 1 by the caller).
// fix_tree() then re-roots/renumbers to the preorder, contiguous-subtree,
// root==0 invariant the algorithm relies on.
Tree* buildFdctTreeFromEdge(const Rcpp::IntegerMatrix& edge, int nTip) {
  int nRow = edge.nrow();
  int nNode = nTip;
  for (int r = 0; r < nRow; r++) {
    if (edge(r, 0) > nNode) nNode = edge(r, 0);
    if (edge(r, 1) > nNode) nNode = edge(r, 1);
  }
  Tree* t = new Tree(static_cast<size_t>(nNode) + 1);
  std::vector<Tree::Node*> apeToNode(nNode + 1, NULL);
  for (int v = 1; v <= nTip; v++) apeToNode[v] = t->add_node(v - 1);  // taxa 0-indexed
  for (int v = nTip + 1; v <= nNode; v++) apeToNode[v] = t->add_node();
  for (int r = 0; r < nRow; r++) apeToNode[edge(r, 0)]->add_child(apeToNode[edge(r, 1)]);
  Tree::Node* root = nRow > 0 ? apeToNode[edge(0, 0)] : apeToNode[1];
  t->fix_tree(root);
  return t;
}

// Integer-label Newick (no trailing ';'), taxa emitted 1-indexed to match the R
// decode labels[as.integer(tip.label)].
void newickInto(Tree::Node* node, std::string& out) {
  if (node->is_leaf()) {
    out += std::to_string(node->taxa + 1);
    return;
  }
  out += '(';
  for (size_t i = 0; i < node->get_children_num(); i++) {
    if (i > 0) out += ',';
    newickInto(node->children[i], out);
  }
  out += ')';
}

}  // namespace

// Frequency-difference consensus of `edgeList` (each a PREORDER ape edge matrix,
// rooted at taxon 1 by the caller) on `nTip` leaves.  Returns an integer-label
// Newick string without a trailing ';'.
// [[Rcpp::export]]
std::string frequencyConsensusCpp(Rcpp::List edgeList, int nTip) {
  int nTree = edgeList.size();
  std::vector<Tree*> trees;
  trees.reserve(nTree);
  for (int i = 0; i < nTree; i++) {
    Rcpp::IntegerMatrix edge = edgeList[i];
    trees.push_back(buildFdctTreeFromEdge(edge, nTip));
  }

  FreqDiff solver(nTip, nTree);
  Tree* consensus = solver.run(trees);  // consumes (deletes) `trees`

  std::string out;
  newickInto(consensus->get_root(), out);
  delete consensus;
  return out;
}
