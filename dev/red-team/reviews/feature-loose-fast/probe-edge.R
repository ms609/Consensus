.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/loose-fast/.rt-lib", .libPaths()))
suppressMessages(library(ConsTree)); suppressMessages(library(TreeTools))

cat("== star-output (fully conflicting) survives read.tree + .RootLikeFirst ==\n")
trees <- c(
  ape::read.tree(text = "((t1,t2),(t3,t4),t5);"),
  ape::read.tree(text = "((t1,t3),(t2,t4),t5);"))
l <- Loose(trees)
cat(sprintf("  class=%s nTip=%d nSplits=%d\n", paste(class(l),collapse=","), NTip(l), NSplits(l)))

cat("\n== 4-leaf minimal (boundary of .PrepareTrees < 4 guard) ==\n")
t4 <- c(ape::read.tree(text="((t1,t2),(t3,t4));"), ape::read.tree(text="((t1,t3),(t2,t4));"))
l4 <- Loose(t4)
cat(sprintf("  nTip=%d nSplits=%d\n", NTip(l4), NSplits(l4)))

cat("\n== single tree / trivial guards (parity with Greedy) ==\n")
one <- ape::read.tree(text="((t1,t2),(t3,t4),t5);")
g1 <- tryCatch(class(Loose(one)), error=function(e) paste("ERR", conditionMessage(e)))
cat(sprintf("  Loose(single phylo) -> %s\n", paste(g1, collapse=",")))
cat(sprintf("  Loose(list(one)) nTip=%d\n", NTip(Loose(list(one)))))

cat("\n== Loose vs Greedy on identical inputs (sanity: same leaf set) ==\n")
trees <- ape::as.phylo(0:20, 9)
labels <- TipLabels(trees[[1]])
cat(sprintf("  Loose nSplits=%d  Greedy nSplits=%d\n", NSplits(Loose(trees)), NSplits(Greedy(trees))))
