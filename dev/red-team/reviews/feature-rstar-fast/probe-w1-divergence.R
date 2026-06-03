# W1 sharpened: prove the subset-leaf case is a SILENT WRONG ANSWER, not just
# "computes on inconsistent input". A subset tree missing leaf X should ABSTAIN
# on every triple involving X. Instead X's id aliases an internal node (~root,
# depth 0), so the tree fabricates an "X is outermost" vote. These differ
# whenever X is NOT outermost in the trees that do mention it.
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/rstar/.agent-cons", .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })
cladeVec <- function(tree) {
  tree <- TreeTools::Preorder(tree); edge <- tree[["edge"]]; nTip <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) { below <- TreeTools::DescendantEdges(edge[,1], edge[,2], edge=e)
    tips <- edge[below,2]; paste(sort(tl[tips[tips <= nTip]]), collapse=",") }, character(1))
  sizes <- lengths(strsplit(cs, ",")); sort(unique(cs[sizes > 1 & sizes < nTip])) }
cs <- function(t) { v <- cladeVec(t); if (length(v)) paste(v, collapse="  ") else "(star)" }

# A: the ONLY tree mentioning t4 puts it sister to t1 (t4 is INNER, not outermost).
A <- read.tree(text = "(((t1,t4),t2),t3);")
# B: a subset tree on {t1,t2,t3} -- no t4. No tree contradicts (t1,t4).
B_subset <- read.tree(text = "((t1,t2),t3);")
# Correct behaviour: t4-triples decided by A alone -> {t1,t4} must survive.
cat("A (only tree with t4) clades        :", cs(A), "\n")
cat("RStar(A, B_subset)  [BUGGY PATH]     :", cs(RStar(list(A, B_subset))), "\n")

# Control: make B explicitly contain t4 as the OUTERMOST taxon. If the bug ==
# 't4 outermost' semantics, this reproduces the buggy output exactly.
B_t4_outer <- read.tree(text = "(((t1,t2),t3),t4);")
cat("RStar(A, B_t4_outermost) [CONTROL]   :", cs(RStar(list(A, B_t4_outer))), "\n")

# Control 2: B explicitly abstains is impossible to express, but A duplicated is
# the 'B carries no conflicting t4 info' limit -> {t1,t4} kept.
cat("RStar(A, A)  [t4 info uncontested]   :", cs(RStar(list(A, A))), "\n")
