// Owen-Provan (2011) GTP geodesic distance and interpolation in BHV tree
// space, plus the Sturm/Miller iterative Frechet mean (Brown & Owen 2020).
//
// A tree is represented by its interior splits, each stored as a *clade*
// bitset canonicalised to EXCLUDE the reference tip (column 0), together with
// that split's length, plus a per-tip vector of leaf (pendant) edge lengths.
// Because compatible cross-tree edges are separated by the vertex-cover step,
// and interpolation only ever rescales existing splits, no split bitset is
// ever newly constructed beyond those present in the input trees.
//
// The algorithm here is validated against the hand-derived 15*sqrt(2) oracle
// of Owen & Provan (2011, Fig. 1 / Example 2) and a suite of invariants; see
// tests/testthat/test-BHV.R.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

namespace bhv {

typedef std::vector<uint64_t> Clade;     // bit w*64 + b set => tip (w*64+b) in clade
static const double TOL = 1e-10;

struct Split { Clade clade; double len; };

struct Tree {
  std::vector<Split> splits;
  std::vector<double> leaf;
  int nTip = 0;
  int nBins = 0;
};

inline bool cl_equal(const Clade& a, const Clade& b) {
  for (size_t w = 0; w < a.size(); ++w) if (a[w] != b[w]) return false;
  return true;
}
inline bool cl_subset(const Clade& a, const Clade& b) {     // a subseteq b
  for (size_t w = 0; w < a.size(); ++w) if (a[w] & ~b[w]) return false;
  return true;
}
inline bool cl_disjoint(const Clade& a, const Clade& b) {
  for (size_t w = 0; w < a.size(); ++w) if (a[w] & b[w]) return false;
  return true;
}
// # nocov start
// Only used by the dead common-edge branch of split_on_common() below (see the
// note there); unreachable for inputs produced by build_geodesic().
inline bool cl_proper_subset(const Clade& a, const Clade& b) {
  return cl_subset(a, b) && !cl_equal(a, b);
}
// # nocov end
inline bool compatible(const Clade& a, const Clade& b) {
  return cl_subset(a, b) || cl_subset(b, a) || cl_disjoint(a, b);
}

inline double vnorm(const std::vector<double>& v) {
  double s = 0; for (double x : v) s += x * x; return std::sqrt(s);
}

// ---- min-weight vertex cover of a bipartite incompatibility graph ----------
// incid[i*nb + j] true iff a_i incompatible with b_j (an edge needing cover).
// Vertex weights are squared edge lengths, normalised so each side sums to 1.
// Returns covers via aCover/bCover and the cover weight (= max-flow value).
struct VC { std::vector<char> aCover, bCover; double weight; };

static VC min_weight_vc(const std::vector<char>& incid,
                        const std::vector<double>& wA,
                        const std::vector<double>& wB) {
  int na = (int)wA.size(), nb = (int)wB.size();
  double sA = 0, sB = 0;
  for (double w : wA) sA += w;
  for (double w : wB) sB += w;
  // Node ids: 0 = source, 1..na = A, na+1..na+nb = B, na+nb+1 = sink
  int S = 0, T = na + nb + 1, N = T + 1;
  const double INF = 1e18, eps = 1e-13;
  std::vector<double> cap((size_t)N * N, 0.0), flow((size_t)N * N, 0.0);
  for (int i = 0; i < na; ++i) cap[(size_t)S * N + (1 + i)] = wA[i] / sA;
  for (int j = 0; j < nb; ++j) cap[(size_t)(1 + na + j) * N + T] = wB[j] / sB;
  for (int i = 0; i < na; ++i)
    for (int j = 0; j < nb; ++j)
      if (incid[(size_t)i * nb + j]) cap[(size_t)(1 + i) * N + (1 + na + j)] = INF;

  std::vector<int> pred(N);
  while (true) {                                   // Edmonds-Karp
    std::fill(pred.begin(), pred.end(), -1);
    pred[S] = S;
    std::vector<int> q; q.push_back(S);
    size_t head = 0; bool found = false;
    while (head < q.size()) {
      int u = q[head++];
      if (u == T) { found = true; break; }
      for (int v = 0; v < N; ++v)
        if (pred[v] < 0 && cap[(size_t)u * N + v] - flow[(size_t)u * N + v] > eps) {
          pred[v] = u; q.push_back(v);
        }
    }
    if (!found) break;
    double b = INF;
    for (int v = T; v != S; v = pred[v]) {
      int u = pred[v];
      b = std::min(b, cap[(size_t)u * N + v] - flow[(size_t)u * N + v]);
    }
    for (int v = T; v != S; v = pred[v]) {
      int u = pred[v];
      flow[(size_t)u * N + v] += b;
      flow[(size_t)v * N + u] -= b;
    }
  }
  std::vector<char> reach(N, 0); reach[S] = 1;
  std::vector<int> q; q.push_back(S); size_t head = 0;
  while (head < q.size()) {
    int u = q[head++];
    for (int v = 0; v < N; ++v)
      if (!reach[v] && cap[(size_t)u * N + v] - flow[(size_t)u * N + v] > eps) {
        reach[v] = 1; q.push_back(v);
      }
  }
  VC out;
  out.aCover.resize(na); out.bCover.resize(nb);
  double weight = 0;
  for (int i = 0; i < na; ++i) { out.aCover[i] = !reach[1 + i];
    if (out.aCover[i]) weight += wA[i] / sA; }
  for (int j = 0; j < nb; ++j) { out.bCover[j] = reach[1 + na + j];
    if (out.bCover[j]) weight += wB[j] / sB; }
  out.weight = weight;
  return out;
}

// ---- one geodesic ----------------------------------------------------------
// A Ratio holds the dropping (e, from tree A) and growing (f, from tree B)
// edges of one support pair, as indices into the group's clade/length arrays.
struct Ratio {
  std::vector<Clade> eC; std::vector<double> eL;   // A side (shrinks)
  std::vector<Clade> fC; std::vector<double> fL;   // B side (grows)
  double eNorm() const { return vnorm(eL); }
  double fNorm() const { return vnorm(fL); }
  double ratio() const {
    double e = eNorm(), f = fNorm();
    if (f == 0) return R_PosInf;
    if (e == 0) return 0.0;
    return e / f;
  }
};

struct CommonEdge { Clade clade; double lenA, lenB; };

struct Geodesic {
  std::vector<Ratio> ratios;          // sorted by increasing ratio value
  std::vector<CommonEdge> common;
  std::vector<double> leafA, leafB;
};

// GTP on a group of mutually non-common edges.
static void gtp_no_common(const std::vector<Clade>& Ac, const std::vector<double>& Al,
                          const std::vector<Clade>& Bc, const std::vector<double>& Bl,
                          std::vector<Ratio>& out) {
  struct QR { std::vector<int> a, b; };
  std::vector<QR> queue;
  QR init;
  for (int i = 0; i < (int)Ac.size(); ++i) init.a.push_back(i);
  for (int j = 0; j < (int)Bc.size(); ++j) init.b.push_back(j);
  queue.push_back(init);
  while (!queue.empty()) {
    QR r = queue.front(); queue.erase(queue.begin());
    int na = (int)r.a.size(), nb = (int)r.b.size();
    bool finalRatio = (na == 0 || nb == 0);
    VC cov;
    if (!finalRatio) {
      std::vector<char> incid((size_t)na * nb);
      for (int ii = 0; ii < na; ++ii)
        for (int jj = 0; jj < nb; ++jj)
          incid[(size_t)ii * nb + jj] = !compatible(Ac[r.a[ii]], Bc[r.b[jj]]);
      std::vector<double> wA(na), wB(nb);
      for (int ii = 0; ii < na; ++ii) wA[ii] = Al[r.a[ii]] * Al[r.a[ii]];
      for (int jj = 0; jj < nb; ++jj) wB[jj] = Bl[r.b[jj]] * Bl[r.b[jj]];
      cov = min_weight_vc(incid, wA, wB);
      // Owen & Provan (2011), Lemma 3.2: a support pair is part of the geodesic
      // iff its min-weight vertex cover has weight >= 1; otherwise an extension
      // (split) shortens the path.  This includes a weight-0 (empty) cover when
      // the two sides are mutually compatible: the pair then splits into
      // (empty, B) and (A, empty), giving the shared-orthant distance
      // sqrt(||A||^2 + ||B||^2) rather than the cone ||A|| + ||B||.
      if (cov.weight >= 1 - TOL) finalRatio = true;
    }
    if (finalRatio) {
      Ratio ro;
      for (int i : r.a) { ro.eC.push_back(Ac[i]); ro.eL.push_back(Al[i]); }
      for (int j : r.b) { ro.fC.push_back(Bc[j]); ro.fL.push_back(Bl[j]); }
      out.push_back(ro);
      continue;
    }
    QR r1, r2;
    for (int ii = 0; ii < na; ++ii) (cov.aCover[ii] ? r1.a : r2.a).push_back(r.a[ii]);
    for (int jj = 0; jj < nb; ++jj) (cov.bCover[jj] ? r2.b : r1.b).push_back(r.b[jj]);
    queue.insert(queue.begin(), r2);
    queue.insert(queue.begin(), r1);
  }
}

// Recursively partition non-common edges into groups separated by common
// edges, running GTP on each group.
static void split_on_common(const std::vector<Clade>& Ac, const std::vector<double>& Al,
                            const std::vector<Clade>& Bc, const std::vector<double>& Bl,
                            std::vector<Ratio>& out) {
  int ci = -1;
  for (int i = 0; i < (int)Ac.size() && ci < 0; ++i)
    for (int j = 0; j < (int)Bc.size(); ++j)
      if (cl_equal(Ac[i], Bc[j])) { ci = i; break; }
  if (ci < 0) {
    if (!Ac.empty() || !Bc.empty()) gtp_no_common(Ac, Al, Bc, Bl, out);
    return;
  }
  // # nocov start
  // Unreachable: build_geodesic() extracts every exact-equal clade into the
  // common-edge list *before* calling split_on_common(), so the Ac/Bc passed
  // here share no clade and `ci` is always < 0 (handled above).  Kept for
  // algorithmic completeness, mirroring the reference GTP recursion.
  const Clade& cl = Ac[ci];
  std::vector<Clade> blA, abA, blB, abB;
  std::vector<double> blAl, abAl, blBl, abBl;
  for (int i = 0; i < (int)Ac.size(); ++i) {
    if (cl_proper_subset(Ac[i], cl)) { blA.push_back(Ac[i]); blAl.push_back(Al[i]); }
    else if (!cl_equal(Ac[i], cl))   { abA.push_back(Ac[i]); abAl.push_back(Al[i]); }
  }
  for (int j = 0; j < (int)Bc.size(); ++j) {
    if (cl_proper_subset(Bc[j], cl)) { blB.push_back(Bc[j]); blBl.push_back(Bl[j]); }
    else if (!cl_equal(Bc[j], cl))   { abB.push_back(Bc[j]); abBl.push_back(Bl[j]); }
  }
  split_on_common(blA, blAl, blB, blBl, out);
  split_on_common(abA, abAl, abB, abBl, out);
  // # nocov end
}

static Geodesic build_geodesic(const Tree& A, const Tree& B) {
  Geodesic g;
  g.leafA = A.leaf; g.leafB = B.leaf;
  std::vector<char> matchedB(B.splits.size(), 0), keepA(A.splits.size(), 1);
  for (size_t i = 0; i < A.splits.size(); ++i)
    for (size_t j = 0; j < B.splits.size(); ++j)
      if (!matchedB[j] && cl_equal(A.splits[i].clade, B.splits[j].clade)) {
        CommonEdge ce; ce.clade = A.splits[i].clade;
        ce.lenA = A.splits[i].len; ce.lenB = B.splits[j].len;
        g.common.push_back(ce);
        keepA[i] = 0; matchedB[j] = 1; break;
      }
  std::vector<Clade> Ac, Bc; std::vector<double> Al, Bl;
  for (size_t i = 0; i < A.splits.size(); ++i)
    if (keepA[i]) { Ac.push_back(A.splits[i].clade); Al.push_back(A.splits[i].len); }
  for (size_t j = 0; j < B.splits.size(); ++j)
    if (!matchedB[j]) { Bc.push_back(B.splits[j].clade); Bl.push_back(B.splits[j].len); }
  split_on_common(Ac, Al, Bc, Bl, g.ratios);
  std::sort(g.ratios.begin(), g.ratios.end(),
            [](const Ratio& x, const Ratio& y) { return x.ratio() < y.ratio(); });
  return g;
}

static double geo_dist(const Geodesic& g) {
  double s = 0;
  for (const Ratio& r : g.ratios) { double v = r.eNorm() + r.fNorm(); s += v * v; }
  for (const CommonEdge& c : g.common) { double d = c.lenA - c.lenB; s += d * d; }
  for (size_t i = 0; i < g.leafA.size(); ++i) {
    double d = g.leafA[i] - g.leafB[i]; s += d * d;
  }
  return std::sqrt(s);
}

// Tree on the geodesic at parameter lambda in [0, 1].
static Tree tree_at(const Geodesic& g, double lambda, int nTip, int nBins) {
  Tree t; t.nTip = nTip; t.nBins = nBins;
  t.leaf.resize(nTip);
  for (int i = 0; i < nTip; ++i)
    t.leaf[i] = (1 - lambda) * g.leafA[i] + lambda * g.leafB[i];
  for (const CommonEdge& c : g.common) {
    double L = (1 - lambda) * c.lenA + lambda * c.lenB;
    if (L > TOL) t.splits.push_back(Split{c.clade, L});
  }
  for (const Ratio& r : g.ratios) {
    double eN = r.eNorm(), fN = r.fNorm();
    double time = (eN + fN > 0) ? eN / (eN + fN) : 0.0;
    if (time >= lambda) {                          // A-side edges present
      double scale = (eN > 0) ? ((1 - lambda) * eN - lambda * fN) / eN : 0.0;
      if (scale > 0)
        for (size_t k = 0; k < r.eC.size(); ++k) {
          double L = scale * r.eL[k];
          if (L > TOL) t.splits.push_back(Split{r.eC[k], L});
        }
    } else {                                        // B-side edges present
      double scale = (fN > 0) ? (lambda * fN - (1 - lambda) * eN) / fN : 0.0;
      if (scale > 0)
        for (size_t k = 0; k < r.fC.size(); ++k) {
          double L = scale * r.fL[k];
          if (L > TOL) t.splits.push_back(Split{r.fC[k], L});
        }
    }
  }
  return t;
}

// Build a Tree from R inputs: membership (nSplit x nTip, 1/0), lengths, leaf.
static Tree from_r(const IntegerMatrix& mem, const NumericVector& len,
                   const NumericVector& leaf) {
  Tree t;
  t.nTip = leaf.size();
  t.nBins = (t.nTip + 63) / 64;
  t.leaf.assign(leaf.begin(), leaf.end());
  int nSplit = mem.nrow();
  for (int s = 0; s < nSplit; ++s) {
    if (len[s] <= TOL) continue;   // a zero-length interior edge is an absent split
    Clade c(t.nBins, 0);
    bool bit0 = (mem(s, 0) != 0);
    for (int tip = 0; tip < t.nTip; ++tip) {
      bool in = (mem(s, tip) != 0);
      if (bit0) in = !in;                          // canonicalise: exclude tip 0
      if (in) c[tip >> 6] |= (uint64_t)1 << (tip & 63);
    }
    t.splits.push_back(Split{c, len[s]});
  }
  return t;
}

// Serialise a Tree's clades back to an R membership matrix (nSplit x nTip).
static List to_r(const Tree& t) {
  int nSplit = (int)t.splits.size();
  IntegerMatrix mem(nSplit, t.nTip);
  NumericVector len(nSplit);
  for (int s = 0; s < nSplit; ++s) {
    const Clade& c = t.splits[s].clade;
    for (int tip = 0; tip < t.nTip; ++tip)
      mem(s, tip) = (c[tip >> 6] >> (tip & 63)) & 1;
    len[s] = t.splits[s].len;
  }
  return List::create(_["membership"] = mem, _["lengths"] = len,
                      _["leaf"] = NumericVector(t.leaf.begin(), t.leaf.end()));
}

} // namespace bhv

// [[Rcpp::export]]
double cpp_bhv_distance(IntegerMatrix memA, NumericVector lenA, NumericVector leafA,
                        IntegerMatrix memB, NumericVector lenB, NumericVector leafB) {
  bhv::Tree A = bhv::from_r(memA, lenA, leafA);
  bhv::Tree B = bhv::from_r(memB, lenB, leafB);
  return bhv::geo_dist(bhv::build_geodesic(A, B));
}

// [[Rcpp::export]]
List cpp_bhv_tree_at(IntegerMatrix memA, NumericVector lenA, NumericVector leafA,
                     IntegerMatrix memB, NumericVector lenB, NumericVector leafB,
                     double lambda) {
  bhv::Tree A = bhv::from_r(memA, lenA, leafA);
  bhv::Tree B = bhv::from_r(memB, lenB, leafB);
  bhv::Geodesic g = bhv::build_geodesic(A, B);
  return bhv::to_r(bhv::tree_at(g, lambda, A.nTip, A.nBins));
}

// Sturm/Miller iterative Frechet mean. `mems`, `lens`, `leaves` are parallel
// lists, one entry per input tree. `tol` is a convergence threshold relative
// to the sample standard deviation. Uses R's RNG for the permutation schedule.
// [[Rcpp::export]]
List cpp_bhv_mean(List mems, List lens, List leaves, int nTip,
                  int maxIter, double tol, int cauchyLength) {
  int r = mems.size();
  std::vector<bhv::Tree> trees(r);
  for (int i = 0; i < r; ++i)
    trees[i] = bhv::from_r(as<IntegerMatrix>(mems[i]), as<NumericVector>(lens[i]),
                           as<NumericVector>(leaves[i]));
  int nBins = (nTip + 63) / 64;
  if (r == 1) return bhv::to_r(trees[0]);

  RNGScope scope;
  std::vector<int> perm(r);
  for (int i = 0; i < r; ++i) perm[i] = i;
  auto shuffle = [&]() {                           // Fisher-Yates with R RNG
    for (int i = r - 1; i > 0; --i) {
      int j = (int)std::floor(R::unif_rand() * (i + 1));
      if (j > i) j = i;
      std::swap(perm[i], perm[j]);
    }
  };
  // One Sturm walk of `iters` steps, starting from a fresh shuffle.
  auto sturmWalk = [&](int iters) -> bhv::Tree {
    shuffle();
    bhv::Tree m = trees[perm[0]];
    for (int i = 1; i < iters; ++i) {
      if (i % r == 0) shuffle();
      m = bhv::tree_at(bhv::build_geodesic(m, trees[perm[i % r]]),
                       1.0 / (i + 1), nTip, nBins);
    }
    return m;
  };

  // Scale the convergence threshold to the data: epsilon = tol * stdDev,
  // estimated from a short 5-round pre-pass (cf. Owen's getEpsilon).
  bhv::Tree test = sturmWalk(5 * r + 1);
  double var = 0;
  for (int i = 0; i < r; ++i) {
    double d = bhv::geo_dist(bhv::build_geodesic(test, trees[i]));
    var += d * d;
  }
  double scale = std::sqrt(var / r);
  double epsilon = tol * scale;
  if (!(epsilon > 0)) epsilon = tol;               // degenerate: identical trees

  // Main run: converge once `cauchyLength` consecutive steps move < epsilon.
  // A step to tree_at(g, lambda) travels exactly lambda * length(g), so the
  // step distance is available without an extra geodesic computation.
  shuffle();
  bhv::Tree mean = trees[perm[0]];
  int consecutive = 0;
  bool converged = false;
  int i = 1;
  for (; i < maxIter && !converged; ++i) {
    if (i % r == 0) shuffle();
    bhv::Geodesic g = bhv::build_geodesic(mean, trees[perm[i % r]]);
    double lambda = 1.0 / (i + 1);
    double stepDist = lambda * bhv::geo_dist(g);
    mean = bhv::tree_at(g, lambda, nTip, nBins);
    consecutive = (stepDist <= epsilon) ? consecutive + 1 : 0;
    if (consecutive >= cauchyLength) converged = true;
  }
  List out = bhv::to_r(mean);
  out["iterations"] = i - 1;
  out["converged"] = converged;
  out["epsilon"] = epsilon;
  return out;
}
