source('C:/Users/pjjg18/GitHub/worktrees/Consensus/adams/dev/red-team/reviews/feature-adams-fast/brute_adams.R')

# Paper Fig.1 sanity: T1,T2,T3 on {a,b,c,d,e,f}; Adams = (b,c,d),(a,e,f)? 
# Fig.1 right tree clades: {b,c,d} and {a,e,f}... let's just self-check brute is
# internally sane on a tiny known case from the docstring example.
t1 <- ape::read.tree(text = "(((a, b), c), d);")
t2 <- ape::read.tree(text = "(((a, b), d), c);")
cat("docstring example brute:", paste(bruteAdamsClades(list(t1,t2)), collapse=" ; "), "\n")
cat("docstring example C++  :", paste(cladeSetLab(Adams(structure(list(t1,t2),class='multiPhylo'))), collapse=" ; "), "\n\n")

# Random ROOTED binary fuzz (matches in-suite). Compare C++ vs INDEPENDENT brute.
set.seed(1)
fails <- 0; tot <- 0
for (n in c(4,5,6,7,8,10,12,15)) {
  for (k in c(2,3,4,5,8)) {
    for (rep in 1:8) {
      labs <- paste0("t", seq_len(n))
      trees <- structure(lapply(seq_len(k), function(i)
        TreeTools::RandomTree(labs, root = TRUE)), class = "multiPhylo")
      cpp <- tryCatch(cladeSetLab(Adams(trees)), error = function(e) paste("ERR:", conditionMessage(e)))
      bru <- bruteAdamsClades(trees)
      tot <- tot + 1
      if (!identical(cpp, bru)) {
        fails <- fails + 1
        if (fails <= 6) {
          cat(sprintf("MISMATCH n=%d k=%d rep=%d\n", n, k, rep))
          cat("  cpp  :", paste(cpp, collapse=" ; "), "\n")
          cat("  brute:", paste(bru, collapse=" ; "), "\n")
          cat("  trees:\n")
          for (tt in trees) cat("    ", ape::write.tree(tt), "\n")
        }
      }
    }
  }
}
cat(sprintf("\nBINARY ROOTED: %d/%d mismatches\n", fails, tot))
