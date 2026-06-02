#include <Rcpp.h>
#include "fact_tree.h"
#include <string>

// Loose (semi-strict / combinable-component) consensus -- FAST PATH PENDING.
//
// The Loose chip ports FACT looseConsensusFast + looseMerge + contract
// (dev/oracle/fact-src/loose.cpp) INTO THIS FILE, reusing fact::Tree,
// buildTreeFromEdge(), newick() and precompute() from fact_tree.h, and fills in
// the body below.  The export is pre-registered in the foundation so that
// RcppExports.* is generated exactly once and the chip never adds an export tag
// (run roxygen2::roxygenise() only -- never compileAttributes()/document()).
//
// [[Rcpp::export]]
std::string looseConsensusCpp(Rcpp::List edgeList, int nTip) {
  (void) edgeList;
  (void) nTip;
  Rcpp::stop("looseConsensusCpp: fast path not yet implemented (foundation stub).");
}
