#include <Rcpp.h>
#include <string>

// Adams consensus -- FAST PATH PENDING.
//
// The Adams chip implements the O(kn log n) algorithm of Jansson, Li & Sung
// (2017) from the paper (no portable reference C++ exists; FACT's adams.cpp is
// the slow recursion and its fast binary mis-prints), replacing the pure-R
// recursion in R/adams.R, INTO THIS FILE, and fills in the body below.  Adams is
// a ROOTED method: the R wrapper passes each tree on its OWN root (it does NOT
// re-root at taxon 1 as the split methods do).  Validate against the classical
// slow Adams clade oracle (dev/oracle: FACT rule 512, rooted=1).  The export is
// pre-registered in the foundation so RcppExports.* is generated once and the
// chip never adds an export tag (run roxygen2::roxygenise() only).
//
// [[Rcpp::export]]
std::string adamsConsensusCpp(Rcpp::List edgeList, int nTip) {
  (void) edgeList;
  (void) nTip;
  Rcpp::stop("adamsConsensusCpp: fast path not yet implemented (foundation stub).");
}
