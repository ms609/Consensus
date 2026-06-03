.libPaths(c('C:/Users/pjjg18/GitHub/Consensus/.agent-cons', .libPaths()))
suppressMessages({library(ConsTree); library(TreeTools); library(ape)})

runcase <- function(label, trees) {
  trees <- lapply(trees, function(x){x$edge.length<-NULL;x}); class(trees)<-"multiPhylo"
  r <- tryCatch({fr<-Frequency(trees); sprintf("OK splits=%d", length(as.Splits(fr)))},
                error=function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("[%s] %s\n", label, r))
}

# Many independent random trees, large n, many trees -> radix 5n adds per layer stress
for (rep in 1:30) {
  set.seed(rep)
  n <- sample(c(60,80,100,150,200,250), 1)
  k <- sample(c(2,3,5,10,20,40), 1)
  trees <- replicate(k, RandomTree(n, root=FALSE), simplify=FALSE)
  runcase(sprintf("seed=%d n=%d k=%d", rep, n, k), trees)
}

# Worst case for special nodes: pectinate t1 with a single deep side per path node,
# crossed with a wildly different pectinate. Repeat at exact powers and +-1.
for (n in c(15,16,17,31,32,33,63,64,65,85,86,128,170,171,200)) {
  perm <- sample(n)
  t1 <- ape::stree(n, "left",  tip.label=paste0("t",seq_len(n)))
  t2 <- ape::stree(n, "left",  tip.label=paste0("t",perm))
  t3 <- as.phylo(ape::stree(n,"balanced", tip.label=paste0("t",sample(n))))
  runcase(sprintf("pect-perm n=%d", n), list(t1,t2,t3,t2,t3))
}
cat("STRESS COMPLETE\n")
