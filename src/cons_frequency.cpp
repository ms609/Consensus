#include <Rcpp.h>
#include <string>

// Frequency-difference consensus -- FAST PATH PENDING.
//
// The Frequency chip ports the near-linear FDCT_new freqdiff implementation
// (freqdiff2.h + lca_preprocessing.h + radix_sort.h + taxas_ranges.h from
// github.com/tswddd2/FDCT_new), making it boost-free (replace dynamic_bitset
// with the vector<ll> word-packing FACT uses), INTO THIS FILE and fills in the
// body below.  The export is pre-registered in the foundation so RcppExports.*
// is generated once and the chip never adds an export tag (run
// roxygen2::roxygenise() only).  Validate against dev/oracle/freqdiff/.
//
// [[Rcpp::export]]
std::string frequencyConsensusCpp(Rcpp::List edgeList, int nTip) {
  (void) edgeList;
  (void) nTip;
  Rcpp::stop("frequencyConsensusCpp: fast path not yet implemented (foundation stub).");
}
