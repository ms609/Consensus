/* Quartet.cpp
 *
 * C++ core for the Quartet algorithm (Takazawa et al. 2026).
 *
 * Finds the tree that maximizes the net concordant quartet information shared
 * with a set of input trees, via a greedy add-and-prune heuristic that may
 * also drop rogue taxa.
 *
 * For each quartet the consensus resolves, the contribution is
 *   agree - penalty * disagree,
 * where, among the input trees that resolve the quartet, `agree` is the number
 * resolving it the consensus's way and `disagree` the number resolving it
 * differently (the "abstain" convention: input polytomies neither agree nor
 * disagree).  Unresolved quartets, and quartets involving a dropped taxon,
 * contribute 0.  The objective is the absolute sum over quartets (no quartet-
 * count normalization), so the star tree scores 0 and the score is calibrated
 * for both resolution and leaf count.  With penalty = 1 (a majority threshold)
 * and binary input trees this recovers the symmetric-quartet-distance median.
 */

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::RawMatrix;
using Rcpp::RawVector;

typedef int_fast16_t int16;
typedef int_fast32_t int32;

static const int16 QC_MAX_TIPS = 100;
static const int16 SPLIT_CHUNK = 8;

// Combinatorial lookup tables (triangular, tetrahedral and hyper-tetrahedral
// numbers) used to index quartets.  In 'Quartet' these live in AllQuartets.cpp
// and are shared via `extern`; here we keep a self-contained copy so this
// translation unit has no external dependencies.
static int32 tri_num[QC_MAX_TIPS + 1];
static int32 tet_num[QC_MAX_TIPS + 1];
static int32 hyp_num[QC_MAX_TIPS + 1];

__attribute__((constructor))  // Construction avoids floating point worries
static void qc_initialize_triangles() {
  tri_num[0] = 0;
  tet_num[0] = 0;
  hyp_num[0] = 0;
  for (int16 i = 0; i != QC_MAX_TIPS; ++i) {
    const int16 nxt = i + 1;
    tri_num[nxt] = tri_num[i] + nxt;
    tet_num[nxt] = tet_num[i] + tri_num[nxt];
    hyp_num[nxt] = hyp_num[i] + tet_num[nxt];
  }
}


// ============================================================================
// Quartet index computation (reuses AllQuartets.cpp formula)
// ============================================================================

// Given 0-based indices a < b < c < d, return the quartet index.
inline int32 quartet_index(int16 a, int16 b, int16 c, int16 d,
                           int16 n_tips) {
  const int16
    choices1 = n_tips - 3,
    choices2 = n_tips - a - 3,
    choices3 = n_tips - b - 2,
    chosen1  = a,
    chosen2  = b - a - 1,
    chosen3  = c - b - 1,
    chosen4  = d - c - 1;
  return (hyp_num[choices1] - hyp_num[choices1 - chosen1])
       + (tet_num[choices2] - tet_num[choices2 - chosen2])
       + (tri_num[choices3] - tri_num[choices3 - chosen3])
       + chosen4;
}

static inline int32 n_quartets(int16 n_tips) {
  return hyp_num[n_tips - 3];
}


// ============================================================================
// Quartet state from a single split for a single quartet
// ============================================================================

// For quartet {i,j,k,l} (i < j < k < l), determine the quartet state
// (0 = unresolved, 1 = il|jk, 2 = jl|ik, 3 = kl|ij) given which side
// of the split each tip is on.
//
// The state encoding matches AllQuartets.cpp:
// - State 1: i and l on same side (ad|bc pattern)
// - State 2: j and l on same side (bd|ac pattern)
// - State 3: k and l on same side (cd|ab pattern)
//
// Returns 0 if the split doesn't have exactly 2 on each side.
inline int quartet_state_from_sides(bool si, bool sj, bool sk, bool sl) {
  int sum = si + sj + sk + sl;
  if (sum != 2) return 0;
  if (si == sl) return 1;  // i,l same side → il|jk
  if (sj == sl) return 2;  // j,l same side → jl|ik
  // sk == sl must be true  → kl|ij
  return 3;
}


// ============================================================================
// Pooled split representation
// ============================================================================

struct PooledSplits {
  int n_splits;
  int n_bytes;   // number of bytes per split in raw representation
  int n_tips;

  // Raw split data: n_splits * n_bytes elements.
  // Split i occupies [i*n_bytes .. (i+1)*n_bytes).
  std::vector<unsigned char> data;

  // Per-split metadata
  std::vector<int> count;       // how many trees contain this split
  std::vector<int> light_side;  // min(popcount, n_tips - popcount)

  // Per-split: list of tip indices on the "1" side (canonical)
  std::vector<std::vector<int16>> tips_on_side1;

  // Tree membership
  std::vector<std::vector<int>> tree_members;

  const unsigned char* split(int i) const {
    return &data[i * n_bytes];
  }
  unsigned char* split(int i) {
    return &data[i * n_bytes];
  }
};


// ============================================================================
// FNV-1a hash for canonical split arrays
// ============================================================================

struct SplitHash {
  int n_bytes;
  explicit SplitHash(int nb) : n_bytes(nb) {}
  SplitHash() : n_bytes(0) {}

  std::size_t operator()(const unsigned char* sp) const {
    std::size_t h = 14695981039346656037ULL;
    for (int i = 0; i < n_bytes; ++i) {
      h ^= static_cast<std::size_t>(sp[i]);
      h *= 1099511628211ULL;
    }
    return h;
  }
};

struct SplitEqual {
  int n_bytes;
  explicit SplitEqual(int nb) : n_bytes(nb) {}
  SplitEqual() : n_bytes(0) {}

  bool operator()(const unsigned char* a, const unsigned char* b) const {
    for (int i = 0; i < n_bytes; ++i) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }
};


// ============================================================================
// pool_splits: deduplicate and canonicalise all splits from all trees
// ============================================================================

static PooledSplits pool_splits(const List& splits_list, int n_tips) {
  const int n_tree = splits_list.size();
  const unsigned char bitmask[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};

  const RawMatrix first_mat = Rcpp::as<RawMatrix>(splits_list[0]);
  const int n_bytes = first_mat.ncol();

  // Mask for the last byte
  const int used_bits = ((n_tips - 1) % 8) + 1;
  const unsigned char last_mask =
    static_cast<unsigned char>((1U << used_bits) - 1U);

  SplitHash hasher(n_bytes);
  SplitEqual eq(n_bytes);
  std::unordered_map<const unsigned char*, int, SplitHash, SplitEqual>
    split_map(64, hasher, eq);

  if (n_bytes < 1) { // # nocov start
    Rcpp::stop("Internal error: n_bytes < 1 in pool_splits (n_tips = %d).",
               n_tips);
  } // # nocov end
  std::vector<unsigned char> canon_buf(n_bytes);

  PooledSplits pool;
  pool.n_tips = n_tips;
  pool.n_bytes = n_bytes;
  pool.n_splits = 0;
  pool.tree_members.resize(n_tree);

  // Reserve pool.data so it never reallocates.  split_map stores raw pointers
  // into this buffer, so reallocation would create dangling keys.
  size_t total_splits = 0;
  for (int t = 0; t < n_tree; ++t) {
    const RawMatrix mat_t = Rcpp::as<RawMatrix>(splits_list[t]);
    total_splits += mat_t.nrow();
  }
  pool.data.reserve(total_splits * n_bytes);

  for (int t = 0; t < n_tree; ++t) {
    const RawMatrix mat = Rcpp::as<RawMatrix>(splits_list[t]);
    const int n_sp = mat.nrow();
    std::vector<int>& members = pool.tree_members[t];
    members.reserve(n_sp);

    for (int s = 0; s < n_sp; ++s) {
      // Copy raw bytes
      for (int b = 0; b < n_bytes; ++b) {
        canon_buf[b] = static_cast<unsigned char>(mat(s, b));
      }
      canon_buf[n_bytes - 1] &= last_mask;

      // Canonicalise: if bit 0 is set, flip
      if (canon_buf[0] & 1U) {
        for (int b = 0; b < n_bytes; ++b) {
          canon_buf[b] = ~canon_buf[b];
        }
        canon_buf[n_bytes - 1] &= last_mask;
      }

      auto it = split_map.find(canon_buf.data());
      int idx;
      if (it != split_map.end()) {
        idx = it->second;
        pool.count[idx]++;
      } else {
        idx = pool.n_splits++;
        const size_t old_sz = pool.data.size();
        pool.data.resize(old_sz + n_bytes);
        std::copy(canon_buf.begin(), canon_buf.end(),
                  pool.data.begin() + old_sz);
        // Popcount
        int pc = 0;
        for (int b = 0; b < n_bytes; ++b) {
          unsigned char byte = canon_buf[b];
          while (byte) { pc += byte & 1; byte >>= 1; }
        }
        pool.count.push_back(1);
        pool.light_side.push_back(std::min(pc, n_tips - pc));

        // Build tip list for side 1 (canonical side: bit 0 is OFF)
        std::vector<int16> tips1;
        for (int16 tip = 0; tip < n_tips; ++tip) {
          if (canon_buf[tip / 8] & bitmask[tip % 8]) {
            tips1.push_back(tip);
          }
        }
        pool.tips_on_side1.push_back(std::move(tips1));

        split_map[pool.split(idx)] = idx;
      }

      // Record unique membership per tree
      bool found = false;
      for (int m : members) {
        if (m == idx) { found = true; break; }
      }
      if (!found) members.push_back(idx);
    }
  }

  return pool;
}


// ============================================================================
// build_quartet_profile: for each quartet, count how many input trees
// resolve it as each of the 3 states (or leave unresolved).
// ============================================================================

// Returns a flat array of 4 * n_q ints:
// profile[q * 4 + s] = count of trees with state s for quartet q.
static std::vector<int> build_quartet_profile(
    const List& splits_list,
    int n_tips
) {
  const unsigned char bitmask[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};
  const int32 n_q = n_quartets(static_cast<int16>(n_tips));
  const int n_tree = splits_list.size();

  std::vector<int> profile(n_q * 4, 0);

  for (int t = 0; t < n_tree; ++t) {
    Rcpp::checkUserInterrupt();
    const RawMatrix splits = Rcpp::as<RawMatrix>(splits_list[t]);
    const int n_sp = splits.nrow();

    // Compute quartet states for this tree (same as quartet_states())
    int32 q = 0;
    for (int16 a = 0; a < n_tips - 3; ++a) {
      const int16 a_mask = bitmask[a % SPLIT_CHUNK];
      const int16 a_chunk = a / SPLIT_CHUNK;
      for (int16 b = a + 1; b < n_tips - 2; ++b) {
        const int16 b_mask = bitmask[b % SPLIT_CHUNK];
        const int16 b_chunk = b / SPLIT_CHUNK;
        for (int16 c = b + 1; c < n_tips - 1; ++c) {
          const int16 c_mask = bitmask[c % SPLIT_CHUNK];
          const int16 c_chunk = c / SPLIT_CHUNK;
          for (int16 d = c + 1; d < n_tips; ++d) {
            const int16 d_mask = bitmask[d % SPLIT_CHUNK];
            const int16 d_chunk = d / SPLIT_CHUNK;
            int state = 0;
            for (int sp = 0; sp < n_sp; ++sp) {
              const bool
                a_state = static_cast<unsigned char>(splits(sp, a_chunk))
                          & a_mask,
                b_state = static_cast<unsigned char>(splits(sp, b_chunk))
                          & b_mask,
                c_state = static_cast<unsigned char>(splits(sp, c_chunk))
                          & c_mask,
                d_state = static_cast<unsigned char>(splits(sp, d_chunk))
                          & d_mask;
              if (a_state) {
                if (b_state) {
                  if (!c_state && !d_state) { state = 3; break; }
                } else {
                  if (c_state) {
                    if (!d_state) { state = 2; break; }
                  } else {
                    if (d_state) { state = 1; break; }
                  }
                }
              } else {
                if (b_state) {
                  if (c_state) {
                    if (!d_state) { state = 1; break; }
                  } else if (d_state) { state = 2; break; }
                } else {
                  if (c_state && d_state) { state = 3; break; }
                }
              }
            }
            profile[q * 4 + state]++;
            q++;
          }
        }
      }
    }
  }

  return profile;
}


// ============================================================================
// compat_mat: pairwise compatibility between pooled splits
// ============================================================================

static std::vector<uint8_t> compat_mat(const PooledSplits& pool) {
  const int M = pool.n_splits;
  const int nb = pool.n_bytes;
  const int n_tips = pool.n_tips;
  const int used_bits = ((n_tips - 1) % 8) + 1;
  const unsigned char last_mask =
    static_cast<unsigned char>((1U << used_bits) - 1U);

  std::vector<uint8_t> compat(M * M, 1);

  for (int i = 0; i < M; ++i) {
    const unsigned char* a = pool.split(i);
    for (int j = i + 1; j < M; ++j) {
      const unsigned char* b = pool.split(j);
      bool ab = false, anb = false, nab = false, nanb = false;
      for (int byte_idx = 0; byte_idx < nb; ++byte_idx) {
        unsigned char mask = (byte_idx == nb - 1) ? last_mask : 0xFF;
        unsigned char a_bin = a[byte_idx] & mask;
        unsigned char b_bin = b[byte_idx] & mask;
        if (!ab)   ab   = (a_bin & b_bin) != 0;
        if (!anb)  anb  = (a_bin & ~b_bin & mask) != 0;
        if (!nab)  nab  = (~a_bin & b_bin & mask) != 0;
        if (!nanb) nanb = (~a_bin & ~b_bin & mask) != 0;
        if (ab && anb && nab && nanb) break;
      }
      bool comp = !ab || !anb || !nab || !nanb;
      compat[i * M + j] = comp ? 1 : 0;
      compat[j * M + i] = comp ? 1 : 0;
    }
  }
  return compat;
}


// ============================================================================
// Greedy state for quartet consensus
// ============================================================================

struct QCGreedyState {
  int M;            // number of pooled splits
  int n_tips;       // total tips (original)
  int n_q;          // C(n_tips, 4) — total quartet slots in profile
  int k;            // number of input trees
  int n_incl;       // currently included splits

  std::vector<uint8_t>  incl;            // which splits are included
  std::vector<int>&     profile;         // quartet profile [n_q * 4]
  const PooledSplits&   pool;
  const std::vector<uint8_t>& compat;

  // Per-quartet state in the current consensus
  std::vector<int>  consensus_state;     // 0 = unresolved, 1-3 = resolved
  std::vector<int>  resolve_count;       // how many included splits resolve it

  // Cached per-split: incompatibility counts with included splits
  std::vector<int>  n_incompat;

  // ---- Objective: net concordant quartet information (signed similarity) ----
  // For each quartet the consensus resolves as state j, the contribution is
  //   agree - penalty * disagree,
  // where agree = count_j (input trees resolving it the consensus's way) and
  // disagree = k_r - count_j (trees resolving it differently), with
  // k_r = k - count_0 the number of trees that resolve the quartet at all (the
  // "abstain" convention: input polytomies neither agree nor disagree).
  // Quartets left unresolved, or involving a dropped taxon, contribute 0.  The
  // score is the absolute sum (no quartet-count normalization), so the star
  // tree scores 0.  penalty (= b/a) is the misinformation penalty: a quartet
  // is worth resolving when count_j / k_r > penalty / (1 + penalty); the
  // default penalty = 1 gives a majority threshold (count_j > k_r / 2).
  double penalty;
  int64_t score_agree;                   // Sum of count_j over resolved quartets
  int64_t score_disagree;                // Sum of (k_r - count_j), resolved q.
  bool can_drop;                         // true when taxon dropping is enabled

  // ---- Taxon-dropping support ----
  std::vector<uint8_t> active_tip;       // 1 = active, 0 = dropped
  std::vector<uint8_t> never_drop;       // 1 = protected
  int n_active;                          // count of active tips
  std::vector<uint8_t> quartet_active;   // 1 = all 4 tips active
  std::vector<int> dropped_tips;         // tip indices in drop order
  std::vector<double> drop_scores;       // objective score after each drop

  QCGreedyState(
      std::vector<int>& profile_,
      const PooledSplits& pool_,
      const std::vector<uint8_t>& compat_,
      int M_, int n_tips_, int k_,
      double penalty_ = 1.0,
      bool can_drop_ = false,
      const std::vector<uint8_t>& never_drop_ = std::vector<uint8_t>()
  ) : M(M_), n_tips(n_tips_),
      n_q(n_quartets(static_cast<int16>(n_tips_))),
      k(k_), n_incl(0),
      incl(M_, 0), profile(profile_), pool(pool_), compat(compat_),
      consensus_state(n_q, 0), resolve_count(n_q, 0),
      n_incompat(M_, 0),
      penalty(penalty_), score_agree(0), score_disagree(0),
      can_drop(can_drop_),
      active_tip(n_tips_, 1),
      never_drop(never_drop_.empty()
                   ? std::vector<uint8_t>(n_tips_, 0)
                   : never_drop_),
      n_active(n_tips_),
      quartet_active(n_q, 1)
  {}

  // Objective score = net concordant quartet information (higher = better).
  double score() const {
    return static_cast<double>(score_agree)
         - penalty * static_cast<double>(score_disagree);
  }

  bool is_compatible(int idx) const {
    return n_incompat[idx] == 0 && n_incl < n_active - 3;
  }

  // ------------------------------------------------------------------
  // Helpers for active-tip-filtered split sides
  // ------------------------------------------------------------------

  // Get active tips on each side of split c.
  void active_split_sides(int c,
                          std::vector<int16>& a_tips0,
                          std::vector<int16>& a_tips1) const {
    const auto& raw_tips1 = pool.tips_on_side1[c];
    a_tips1.clear();
    for (int16 tip : raw_tips1) {
      if (active_tip[tip]) a_tips1.push_back(tip);
    }
    a_tips0.clear();
    int idx1 = 0;
    const int am1 = static_cast<int>(a_tips1.size());
    for (int16 tip = 0; tip < n_tips; ++tip) {
      if (!active_tip[tip]) continue;
      if (idx1 < am1 && a_tips1[idx1] == tip) {
        idx1++;
      } else {
        a_tips0.push_back(tip);
      }
    }
  }

  // Check whether split c is trivial among active tips.
  bool is_trivial_active(int c) const {
    int active_on_1 = 0;
    for (int16 tip : pool.tips_on_side1[c]) {
      if (active_tip[tip]) active_on_1++;
    }
    return active_on_1 < 2 || (n_active - active_on_1) < 2;
  }

  // ------------------------------------------------------------------
  // Benefit/execute for adding a split
  // ------------------------------------------------------------------

  double add_benefit(int c) const {
    const unsigned char bitmask[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};
    const unsigned char* sp = pool.split(c);

    std::vector<int16> a_tips0, a_tips1;
    active_split_sides(c, a_tips0, a_tips1);
    const int m0 = static_cast<int>(a_tips0.size());
    const int m1 = static_cast<int>(a_tips1.size());
    if (m0 < 2 || m1 < 2) return -1e30;

    int64_t agree_sum = 0, disagree_sum = 0;

    for (int ai = 0; ai < m0 - 1; ++ai) {
      for (int bi = ai + 1; bi < m0; ++bi) {
        for (int ci = 0; ci < m1 - 1; ++ci) {
          for (int di = ci + 1; di < m1; ++di) {
            int16 t[4] = {a_tips0[ai], a_tips0[bi], a_tips1[ci], a_tips1[di]};
            for (int x = 1; x < 4; ++x) {
              int16 key = t[x];
              int y = x - 1;
              while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
              t[y + 1] = key;
            }

            int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                     static_cast<int16>(n_tips));
            if (resolve_count[qi] > 0) continue;

            bool sa = sp[t[0] / 8] & bitmask[t[0] % 8];
            bool sb = sp[t[1] / 8] & bitmask[t[1] % 8];
            bool sc = sp[t[2] / 8] & bitmask[t[2] % 8];
            bool sd = sp[t[3] / 8] & bitmask[t[3] % 8];
            int state = quartet_state_from_sides(sa, sb, sc, sd);
            if (state == 0) continue;

            int count_j = profile[qi * 4 + state];
            int count_0 = profile[qi * 4 + 0];
            agree_sum += count_j;
            disagree_sum += (k - count_0) - count_j;
          }
        }
      }
    }

    return static_cast<double>(agree_sum)
         - penalty * static_cast<double>(disagree_sum);
  }

  void do_add(int c) {
    const unsigned char bitmask[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};
    incl[c] = 1;
    n_incl++;

    for (int j = 0; j < M; ++j) {
      if (!compat[j * M + c]) n_incompat[j]++;
    }

    const unsigned char* sp = pool.split(c);
    std::vector<int16> a_tips0, a_tips1;
    active_split_sides(c, a_tips0, a_tips1);
    const int m0 = static_cast<int>(a_tips0.size());
    const int m1 = static_cast<int>(a_tips1.size());

    for (int ai = 0; ai < m0 - 1; ++ai) {
      for (int bi = ai + 1; bi < m0; ++bi) {
        for (int ci = 0; ci < m1 - 1; ++ci) {
          for (int di = ci + 1; di < m1; ++di) {
            int16 t[4] = {a_tips0[ai], a_tips0[bi], a_tips1[ci], a_tips1[di]};
            for (int x = 1; x < 4; ++x) {
              int16 key = t[x];
              int y = x - 1;
              while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
              t[y + 1] = key;
            }

            int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                     static_cast<int16>(n_tips));

            bool sa = sp[t[0] / 8] & bitmask[t[0] % 8];
            bool sb = sp[t[1] / 8] & bitmask[t[1] % 8];
            bool sc = sp[t[2] / 8] & bitmask[t[2] % 8];
            bool sd = sp[t[3] / 8] & bitmask[t[3] % 8];
            int state = quartet_state_from_sides(sa, sb, sc, sd);
            if (state == 0) continue;

            resolve_count[qi]++;
            if (resolve_count[qi] == 1) {
              consensus_state[qi] = state;
              int count_j = profile[qi * 4 + state];
              int count_0 = profile[qi * 4 + 0];
              score_agree += count_j;
              score_disagree += (k - count_0) - count_j;
            }
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // Benefit/execute for removing a split
  // ------------------------------------------------------------------

  double remove_benefit(int c) const {
    std::vector<int16> a_tips0, a_tips1;
    active_split_sides(c, a_tips0, a_tips1);
    const int m0 = static_cast<int>(a_tips0.size());
    const int m1 = static_cast<int>(a_tips1.size());

    int64_t agree_sum = 0, disagree_sum = 0;

    for (int ai = 0; ai < m0 - 1; ++ai) {
      for (int bi = ai + 1; bi < m0; ++bi) {
        for (int ci = 0; ci < m1 - 1; ++ci) {
          for (int di = ci + 1; di < m1; ++di) {
            int16 t[4] = {a_tips0[ai], a_tips0[bi], a_tips1[ci], a_tips1[di]};
            for (int x = 1; x < 4; ++x) {
              int16 key = t[x];
              int y = x - 1;
              while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
              t[y + 1] = key;
            }

            int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                     static_cast<int16>(n_tips));
            if (resolve_count[qi] != 1) continue;
            if (consensus_state[qi] == 0) continue;

            int state = consensus_state[qi];
            int count_j = profile[qi * 4 + state];
            int count_0 = profile[qi * 4 + 0];
            agree_sum += count_j;
            disagree_sum += (k - count_0) - count_j;
          }
        }
      }
    }

    // Benefit of removing = -(value of the quartets that become unresolved)
    return penalty * static_cast<double>(disagree_sum)
         - static_cast<double>(agree_sum);
  }

  void do_remove(int c) {
    const unsigned char bitmask[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};
    incl[c] = 0;
    n_incl--;

    for (int j = 0; j < M; ++j) {
      if (!compat[j * M + c]) n_incompat[j]--;
    }

    const unsigned char* sp = pool.split(c);
    std::vector<int16> a_tips0, a_tips1;
    active_split_sides(c, a_tips0, a_tips1);
    const int m0 = static_cast<int>(a_tips0.size());
    const int m1 = static_cast<int>(a_tips1.size());

    for (int ai = 0; ai < m0 - 1; ++ai) {
      for (int bi = ai + 1; bi < m0; ++bi) {
        for (int ci = 0; ci < m1 - 1; ++ci) {
          for (int di = ci + 1; di < m1; ++di) {
            int16 t[4] = {a_tips0[ai], a_tips0[bi], a_tips1[ci], a_tips1[di]};
            for (int x = 1; x < 4; ++x) {
              int16 key = t[x];
              int y = x - 1;
              while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
              t[y + 1] = key;
            }

            int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                     static_cast<int16>(n_tips));

            bool sa = sp[t[0] / 8] & bitmask[t[0] % 8];
            bool sb = sp[t[1] / 8] & bitmask[t[1] % 8];
            bool sc = sp[t[2] / 8] & bitmask[t[2] % 8];
            bool sd = sp[t[3] / 8] & bitmask[t[3] % 8];
            int state = quartet_state_from_sides(sa, sb, sc, sd);
            if (state == 0) continue;

            resolve_count[qi]--;
            if (resolve_count[qi] == 0) {
              int count_j = profile[qi * 4 + consensus_state[qi]];
              int count_0 = profile[qi * 4 + 0];
              score_agree -= count_j;
              score_disagree -= (k - count_0) - count_j;
              consensus_state[qi] = 0;
            }
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // Taxon dropping
  // ------------------------------------------------------------------

  // Benefit of dropping a single tip = the increase in the objective score.
  // Dropping removes every quartet containing the tip; each such quartet that
  // is currently resolved contributes (agree - penalty * disagree) to the
  // score, so dropping changes the score by minus that sum.  Dropping helps
  // exactly when the tip's resolved quartets are, on balance, misleading
  // (their net contribution is negative).  Quartets the tip leaves unresolved
  // contribute 0 and so neither favour nor oppose dropping.
  double drop_benefit(int16 tip) const {
    if (!can_drop || !active_tip[tip] || never_drop[tip]) return -1e30;
    if (n_active <= 4) return -1e30;

    // Build list of other active tips
    std::vector<int16> others;
    others.reserve(n_active - 1);
    for (int16 t = 0; t < n_tips; ++t) {
      if (active_tip[t] && t != tip) others.push_back(t);
    }
    const int nO = static_cast<int>(others.size());

    int64_t agree_sum = 0, disagree_sum = 0;
    for (int ai = 0; ai < nO - 2; ++ai) {
      for (int bi = ai + 1; bi < nO - 1; ++bi) {
        for (int ci = bi + 1; ci < nO; ++ci) {
          int16 t[4] = {tip, others[ai], others[bi], others[ci]};
          for (int x = 1; x < 4; ++x) {
            int16 key = t[x];
            int y = x - 1;
            while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
            t[y + 1] = key;
          }
          int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                   static_cast<int16>(n_tips));
          if (consensus_state[qi] == 0) continue;
          int count_j = profile[qi * 4 + consensus_state[qi]];
          int count_0 = profile[qi * 4 + 0];
          agree_sum += count_j;
          disagree_sum += (k - count_0) - count_j;
        }
      }
    }

    return penalty * static_cast<double>(disagree_sum)
         - static_cast<double>(agree_sum);
  }

  // Execute: drop a tip.  Subtract every resolved quartet containing it from
  // the score, mark those quartets inactive, then prune any included split
  // that has become trivial among the remaining active tips.
  void do_drop(int16 tip) {
    active_tip[tip] = 0;

    std::vector<int16> others;
    others.reserve(n_active - 1);
    for (int16 t = 0; t < n_tips; ++t) {
      if (active_tip[t] && t != tip) others.push_back(t);
    }
    const int nO = static_cast<int>(others.size());

    for (int ai = 0; ai < nO - 2; ++ai) {
      for (int bi = ai + 1; bi < nO - 1; ++bi) {
        for (int ci = bi + 1; ci < nO; ++ci) {
          int16 t[4] = {tip, others[ai], others[bi], others[ci]};
          for (int x = 1; x < 4; ++x) {
            int16 key = t[x];
            int y = x - 1;
            while (y >= 0 && t[y] > key) { t[y + 1] = t[y]; y--; }
            t[y + 1] = key;
          }
          int32 qi = quartet_index(t[0], t[1], t[2], t[3],
                                   static_cast<int16>(n_tips));
          if (consensus_state[qi] != 0) {
            int count_j = profile[qi * 4 + consensus_state[qi]];
            int count_0 = profile[qi * 4 + 0];
            score_agree -= count_j;
            score_disagree -= (k - count_0) - count_j;
          }
          quartet_active[qi] = 0;
          resolve_count[qi] = 0;
          consensus_state[qi] = 0;
        }
      }
    }

    n_active--;

    // Remove included splits that became trivial among the remaining tips
    for (int s = 0; s < M; ++s) {
      if (incl[s] && is_trivial_active(s)) {
        do_remove(s);
      }
    }

    dropped_tips.push_back(tip);
    drop_scores.push_back(score());
  }
};


// ============================================================================
// Greedy "best" strategy
// ============================================================================

// Move type for the greedy loop
enum MoveType { SPLIT_ADD, SPLIT_REMOVE, TIP_DROP };

// Tolerance for benefit comparisons.  Benefits are now absolute (un-
// normalized) sums of the form agree - penalty * disagree, so for the default
// penalty = 1 (and any dyadic penalty) they are exact integers/half-integers
// computed from int64 sub-sums; the smallest real improvement is therefore
// >= a small fraction, while floating noise from the single penalty multiply
// stays well below 1e-6.
static const double BENEFIT_TOL = 1e-6;

static void greedy_best(QCGreedyState& st,
                        const std::vector<int>& sort_ord) {
  while (true) {
    Rcpp::checkUserInterrupt();
    double best_ben = 0.0;
    int best_idx = -1;
    MoveType best_type = SPLIT_ADD;

    // Evaluate split adds/removes
    for (int si = 0; si < st.M; ++si) {
      int idx = sort_ord[si];
      if (st.incl[idx]) {
        double ben = st.remove_benefit(idx);
        if (ben > best_ben) {
          best_ben = ben;
          best_idx = idx;
          best_type = SPLIT_REMOVE;
        }
      } else {
        if (!st.is_compatible(idx)) continue;
        double ben = st.add_benefit(idx);
        if (ben > best_ben) {
          best_ben = ben;
          best_idx = idx;
          best_type = SPLIT_ADD;
        }
      }
    }

    // Evaluate taxon drops
    if (st.can_drop && st.n_active > 4) {
      for (int16 tip = 0; tip < st.n_tips; ++tip) {
        if (!st.active_tip[tip] || st.never_drop[tip]) continue;
        double ben = st.drop_benefit(tip);
        if (ben > best_ben) {
          best_ben = ben;
          best_idx = tip;
          best_type = TIP_DROP;
        }
      }
    }

    if (best_ben <= BENEFIT_TOL || best_idx < 0) break;

    switch (best_type) {
      case SPLIT_ADD:    st.do_add(best_idx); break;
      case SPLIT_REMOVE: st.do_remove(best_idx); break;
      case TIP_DROP:     st.do_drop(static_cast<int16>(best_idx)); break;
    }
  }
}


// ============================================================================
// Greedy "first" strategy
// ============================================================================

static void greedy_first(QCGreedyState& st,
                         const std::vector<int>& sort_ord) {
  bool improving = true;
  while (improving) {
    Rcpp::checkUserInterrupt();
    improving = false;

    // Try split moves first
    for (int si = 0; si < st.M; ++si) {
      int idx = sort_ord[si];
      if (st.incl[idx]) {
        if (st.remove_benefit(idx) > BENEFIT_TOL) {
          st.do_remove(idx);
          improving = true;
          break;
        }
      } else {
        if (!st.is_compatible(idx)) continue;
        if (st.add_benefit(idx) > BENEFIT_TOL) {
          st.do_add(idx);
          improving = true;
          break;
        }
      }
    }
    if (improving) continue;

    // Try taxon drops
    if (st.can_drop && st.n_active > 4) {
      for (int16 tip = 0; tip < st.n_tips; ++tip) {
        if (!st.active_tip[tip] || st.never_drop[tip]) continue;
        if (st.drop_benefit(tip) > BENEFIT_TOL) {
          st.do_drop(tip);
          improving = true;
          break;
        }
      }
    }
  }
}


// ============================================================================
// Main exported function
// ============================================================================

//' Quartet consensus (C++ implementation)
//'
//' @param splits_list List of raw matrices (one per tree), from as.Splits().
//' @param n_tips Number of tips.
//' @param init_majority Logical: TRUE to start from majority-rule splits.
//' @param init_extended Logical: TRUE to start from extended majority splits.
//' @param greedy_best_flag Logical: TRUE for "best", FALSE for "first".
//' @param never_drop_r Integer vector (1-based) of tip indices that must not
//'   be dropped, or integer(0) to allow all drops.  If NULL, taxon dropping
//'   is disabled.
//' @param penalty_r Double: the misinformation penalty b/a.  A quartet is
//'   resolved when its support among resolving input trees exceeds
//'   penalty / (1 + penalty); the default 1 gives a majority threshold.
//'
//' @return A list with `splits` (raw matrix of non-trivial splits remapped
//'   to active tips), `n_active` (integer), `active_tips` (logical),
//'   `dropped_tips` (integer, 1-based), and `drop_scores` (double).
//' @keywords internal
// [[Rcpp::export]]
List cpp_quartet_consensus(
    const List& splits_list,
    const int n_tips,
    const bool init_majority,
    const bool init_extended,
    const bool greedy_best_flag,
    Rcpp::Nullable<Rcpp::IntegerVector> never_drop_r = R_NilValue,
    const double penalty_r = 1.0
) {
  if (n_tips > QC_MAX_TIPS) {
    Rcpp::stop("Quartet supports at most %d tips.", QC_MAX_TIPS);
  }
  if (n_tips < 4) {
    Rcpp::stop("Need at least 4 tips for quartet consensus.");
  }

  const int n_tree = splits_list.size();

  // ---- Taxon dropping configuration ----
  bool dropping = never_drop_r.isNotNull();

  std::vector<uint8_t> never_drop_vec(n_tips, 0);
  if (dropping) {
    Rcpp::IntegerVector nd(never_drop_r);
    for (int i = 0; i < nd.size(); ++i) {
      int idx = nd[i] - 1;  // convert from 1-based R to 0-based C++
      if (idx >= 0 && idx < n_tips) never_drop_vec[idx] = 1;
    }
  }

  // ---- Pool unique splits ----
  PooledSplits pool = pool_splits(splits_list, n_tips);
  const int M = pool.n_splits;

  if (M == 0) {
    return List::create(
      Rcpp::Named("splits") = RawMatrix(0, 0),
      Rcpp::Named("n_active") = n_tips,
      Rcpp::Named("dropped_tips") = IntegerVector(0),
      Rcpp::Named("drop_scores") = Rcpp::NumericVector(0),
      Rcpp::Named("active_tips") = LogicalVector(n_tips, true)
    );
  }

  // ---- Build quartet profile ----
  std::vector<int> profile = build_quartet_profile(splits_list, n_tips);

  // ---- Compatibility matrix ----
  std::vector<uint8_t> compat = compat_mat(pool);

  // ---- Sort order (by count descending) ----
  std::vector<int> sort_ord(M);
  std::iota(sort_ord.begin(), sort_ord.end(), 0);
  std::sort(sort_ord.begin(), sort_ord.end(),
            [&](int a, int b) { return pool.count[a] > pool.count[b]; });

  // ---- Initialize greedy state ----
  QCGreedyState st(profile, pool, compat, M, n_tips, n_tree,
                   penalty_r, dropping, never_drop_vec);

  // ---- Init from majority or extended majority ----
  if (init_majority || init_extended) {
    double half = n_tree / 2.0;

    if (init_extended) {
      for (int si = 0; si < M; ++si) {
        int idx = sort_ord[si];
        if (pool.count[idx] > half) {
          st.do_add(idx);
        } else if (st.is_compatible(idx)) {
          st.do_add(idx);
        }
      }
    } else {
      for (int i = 0; i < M; ++i) {
        if (pool.count[i] > half) {
          st.do_add(i);
        }
      }
    }
  }

  // ---- Greedy loop ----
  if (greedy_best_flag) {
    greedy_best(st, sort_ord);
  } else {
    greedy_first(st, sort_ord);
  }

  // ---- Build output: remap splits to active tips ----

  // Map original tip indices to active-only indices (0-based)
  std::vector<int> tip_remap(n_tips, -1);
  int n_active_out = 0;
  for (int i = 0; i < n_tips; ++i) {
    if (st.active_tip[i]) {
      tip_remap[i] = n_active_out++;
    }
  }
  const int n_bytes_new = (n_active_out + 7) / 8;

  // Collect included splits, remap bitvectors, drop trivials
  std::vector<std::vector<unsigned char>> out_splits;
  out_splits.reserve(M);

  const unsigned char bitmask_out[8] = {1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U};

  for (int i = 0; i < M; ++i) {
    if (!st.incl[i]) continue;
    const unsigned char* src = pool.split(i);

    // Repack into active-tip-only bitvector
    std::vector<unsigned char> row(n_bytes_new, 0);
    int popcount = 0;
    for (int tip = 0; tip < n_tips; ++tip) {
      if (tip_remap[tip] < 0) continue;
      if (src[tip / 8] & bitmask_out[tip % 8]) {
        int new_idx = tip_remap[tip];
        row[new_idx / 8] |= bitmask_out[new_idx % 8];
        popcount++;
      }
    }

    // Keep only non-trivial splits (>= 2 on each side)
    if (popcount >= 2 && (n_active_out - popcount) >= 2) {
      out_splits.push_back(std::move(row));
    }
  }

  const int n_out = static_cast<int>(out_splits.size());
  RawMatrix splits_r(n_out, n_bytes_new);
  for (int i = 0; i < n_out; ++i) {
    for (int j = 0; j < n_bytes_new; ++j) {
      splits_r(i, j) = Rbyte(out_splits[i][j]);
    }
  }

  // Dropped tips (convert to 1-based for R)
  IntegerVector dropped_r(st.dropped_tips.size());
  for (size_t i = 0; i < st.dropped_tips.size(); ++i) {
    dropped_r[i] = st.dropped_tips[i] + 1;
  }
  Rcpp::NumericVector scores_r(st.drop_scores.begin(), st.drop_scores.end());

  LogicalVector active_r(n_tips);
  for (int i = 0; i < n_tips; ++i) active_r[i] = st.active_tip[i] != 0;

  return List::create(
    Rcpp::Named("splits") = splits_r,
    Rcpp::Named("n_active") = n_active_out,
    Rcpp::Named("dropped_tips") = dropped_r,
    Rcpp::Named("drop_scores") = scores_r,
    Rcpp::Named("active_tips") = active_r
  );
}
