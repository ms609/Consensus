#include <Rcpp.h>
#include "fact_tree.h"
#include <string>

// Majority-rule (+) consensus -- FAST PATH PENDING.
//
// The MajorityPlus chip ports FACT majorityPlusConsensus + updateCounter +
// majorityPlusMerge + majContract (dev/oracle/fact-src/majorityplus.cpp) INTO
// THIS FILE, reusing fact::Tree, buildTreeFromEdge(), newick() and precompute()
// from fact_tree.h, and fills in the body below.  The export is pre-registered
// in the foundation so RcppExports.* is generated once and the chip never adds
// an export tag (run roxygen2::roxygenise() only).
//
// [[Rcpp::export]]
std::string majorityPlusConsensusCpp(Rcpp::List edgeList, int nTip) {
  (void) edgeList;
  (void) nTip;
  Rcpp::stop("majorityPlusConsensusCpp: fast path not yet implemented (foundation stub).");
}
