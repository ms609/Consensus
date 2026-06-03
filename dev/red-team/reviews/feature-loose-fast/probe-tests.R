.libPaths(c("C:/Users/pjjg18/GitHub/worktrees/Consensus/loose-fast/.rt-lib", .libPaths()))
suppressMessages(library(ConsTree)); suppressMessages(library(TreeTools))

cat("== test 'Every loose split compatible' — NSplits per spec (next-guard probe) ==\n")
for (spec in list(0:20, c(0,0,0,1,2,53,99), 1:15)) {
  trees <- ape::as.phylo(spec, 9)
  labels <- TipLabels(trees[[1]])
  ls <- as.Splits(Loose(trees), tipLabels = labels)
  cat(sprintf("  spec len=%-2d  NSplits=%d  %s\n", length(spec), NSplits(ls),
              if (NSplits(ls) == 0L) "<-- STAR: compat assertion SKIPPED (vacuous)" else ""))
}

cat("\n== test 1 specs (strict<=loose) — NSplits ==\n")
for (spec in list(0:20, 0:30, c(0,0,0,1,2,53,99), 1:15)) {
  trees <- ape::as.phylo(spec, 9)
  labels <- TipLabels(trees[[1]])
  l <- Loose(trees)
  s <- Strict(trees)
  cat(sprintf("  spec len=%-2d  loose=%d strict=%d\n", length(spec), NSplits(l), NSplits(s)))
}

cat("\n== test 'drops contradicted split' — does loose actually become a star here? ==\n")
trees <- c(
  ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
  ape::read.tree(text = "((t1, t2), (t3, t4), t5);"),
  ape::read.tree(text = "((t1, t3), (t2, t4), t5);"))
labels <- TipLabels(trees[[1]])
l <- Loose(trees)
cat(sprintf("  Loose nSplits=%d  splits=%s\n", NSplits(l),
            paste(as.character(as.Splits(l, tipLabels=labels)), collapse=" | ")))
