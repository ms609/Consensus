# Probe: resolve rooted-clade vs unrooted-split semantics and test the public
# building blocks for an all-R / hybrid implementation of the selection methods.
.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages(library(TreeTools))

trees <- ape::as.phylo(0:5, 8)

tt50 <- Consensus(trees, p = 0.5)
ac50 <- ape::consensus(trees, p = 0.5)

cat("== rooted/unrooted equivalence ==\n")
cat("TreeTools Consensus NSplits:", NSplits(tt50),
    "rooted:", ape::is.rooted(tt50), "\n")
cat("ape consensus     NSplits:", NSplits(ac50),
    "rooted:", ape::is.rooted(ac50), "\n")
cat("RF distance TT vs ape:", TreeDist::RobinsonFoulds(tt50, ac50), "\n")

cat("\n== as.Splits shapes ==\n")
spTree <- as.Splits(trees[[1]])
cat("as.Splits(phylo): class", class(spTree), "length", length(spTree), "\n")
spMulti <- as.Splits(trees)
cat("as.Splits(multiPhylo): class", class(spMulti),
    "length", if (is.list(spMulti)) length(spMulti) else NA, "\n")

cat("\n== reconstruction round-trip ==\n")
recon <- as.phylo(as.Splits(tt50))
cat("as.phylo(as.Splits(tt50)) RF to tt50:",
    TreeDist::RobinsonFoulds(recon, tt50), "\n")

cat("\n== CompatibleSplits shape ==\n")
cs <- CompatibleSplits(spTree)
cat("CompatibleSplits(splits): class", class(cs),
    "dim", paste(dim(as.matrix(cs)), collapse = "x"), "\n")

cat("\n== pooled-splits-with-counts via public API ==\n")
pooled <- do.call(c, lapply(trees, as.Splits, tipLabels = TipLabels(trees[[1]])))
cat("pooled length:", length(pooled), "\n")
dup <- duplicated(pooled)
cat("unique splits:", sum(!dup), "\n")
