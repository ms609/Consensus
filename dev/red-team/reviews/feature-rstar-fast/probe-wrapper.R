# Continuation of the R* review: probe the R WRAPPER and DEGENERATE INPUTS.
# The core C++ is exhaustively validated; the sibling reviewer (same author)
# found every real bug in the wrapper / degenerate inputs, so that is where the
# residual risk lives.  We do NOT pre-align labels here (real users won't).
# Each case reports: clean error (GOOD - fail loud), sensible output, or
# SILENT-WRONG / crash (BAD).
.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/rstar/.agent-cons", .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })

cladeSet <- function(tree) {
  if (is.null(tree) || !inherits(tree, "phylo")) return(NA_character_)
  tree <- TreeTools::Preorder(tree); edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]
  cs <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]; paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sizes <- lengths(strsplit(cs, ",")); paste(sort(unique(cs[sizes > 1 & sizes < nTip])), collapse=" ")
}
descr <- function(x) {
  if (inherits(x, "probe_err"))  return(paste0("ERROR: ", as.character(x)))
  if (is.null(x))                return("NULL")
  if (inherits(x, "phylo"))      return(sprintf("phylo n=%d rooted=%s clades=[%s]",
                                  NTip(x), is.rooted(x), cladeSet(x)))
  paste0("OTHER: ", class(x)[1])
}
run <- function(label, expr) {
  x <- tryCatch(force(expr),
                error = function(e) structure(conditionMessage(e), class = "probe_err"))
  cat(sprintf("%-34s -> %s\n", label, descr(x)))
  invisible(x)
}
nl <- function() cat("\n")

# Reference topology, two label vectors (same SET) and a shuffled-vector copy.
t_abcd <- read.tree(text = "(((a,b),c),d);")            # rooted, 4 tips

cat("==== k < 2 and n < 3 ====\n")
run("k=1 (single tree)",        RStar(list(t_abcd)))
run("k=0 (empty list)",         RStar(list()))
run("k=0 (multiPhylo len 0)",   RStar(structure(list(), class = "multiPhylo")))
run("n=2 (two 2-tip trees)",    RStar(list(read.tree(text="(a,b);"), read.tree(text="(b,a);"))))
nl()

cat("==== normal multi-tree: same label SET, DIFFERENT tip-vector order ====\n")
# Build 3 identical-topology trees but with permuted tip.label vectors; the
# wrapper relabels each to trees[[1]]'s order. Correct answer: the shared clades.
set.seed(1)
mk_perm <- function(tmpl) { tr <- tmpl; p <- sample(NTip(tmpl)); tr$tip.label <- tmpl$tip.label[p]
  tr$edge[tr$edge[,2] <= NTip(tmpl), 2] <- match(seq_len(NTip(tmpl)), p)[tr$edge[tr$edge[,2] <= NTip(tmpl),2]]; tr }
big <- RandomTree(paste0("t", 1:9), root = TRUE)
perm3 <- list(big, mk_perm(big), mk_perm(big))
run("3x same tree, permuted vectors", RStar(perm3))   # expect clades == cladeSet(big)
cat(sprintf("%-34s    expect clades=[%s]\n", "", cladeSet(big)))
nl()

cat("==== MISMATCHED label sets across trees (no pre-align) ====\n")
t1 <- read.tree(text = "(((t1,t2),t3),t4);")
t2_extra   <- read.tree(text = "(((t1,t2),t3),zz);")   # 'zz' not in tree 1; 't4' missing
t2_more    <- read.tree(text = "((((t1,t2),t3),t4),t5);") # 5 tips vs 4
t2_fewer   <- read.tree(text = "((t1,t2),t3);")          # 3 tips vs 4
run("label set differs (t4 vs zz)", RStar(list(t1, t2_extra)))
run("tree2 has MORE tips (5 vs 4)", RStar(list(t1, t2_more)))
run("tree2 has FEWER tips (3 vs 4)", RStar(list(t1, t2_fewer)))
nl()

cat("==== duplicate tip labels ====\n")
t_dup <- read.tree(text = "(((a,a),b),c);")
run("tree with duplicate label 'a'", RStar(list(t_dup, t_dup)))
nl()

cat("==== unrooted / star / basal polytomy input ====\n")
set.seed(2)
unr <- lapply(1:4, function(i) ape::rtree(8, rooted = FALSE))   # unrooted (basal trifurcation)
unr <- lapply(unr, function(tr) { tr$edge.length <- NULL; tr })
run("4x unrooted 8-tip trees", RStar(unr))
star <- read.tree(text = "(a,b,c,d,e,f,g,h);")
run("3x star (root polytomy)", RStar(list(star, star, star)))
basal <- read.tree(text = "((a,b),c,d,e);")             # basal polytomy + one cherry
run("3x basal-polytomy tree", RStar(list(basal, basal, basal)))
nl()

cat("==== input container variants ====\n")
run("multiPhylo input", RStar(structure(list(t_abcd, t_abcd), class = "multiPhylo")))
run("list with a NULL element", RStar(list(t_abcd, NULL, t_abcd)))
run("non-list (a number)", RStar(42))
nl()
cat("==== done ====\n")
