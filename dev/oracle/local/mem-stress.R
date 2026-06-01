# Memory-safety stress for the compiled Local() core.
# phangorn is NEVER loaded here, so any crash isolates a heap bug in
# localConsensus() from the phangorn/R-devel ABI segfault.
# Env: MEM_NS = comma-separated leaf counts (default "10,12,14,16");
#      MEM_GCTORTURE = "1" to run under gctorture (surfaces PROTECT bugs);
#      MEM_REPS = inner repetitions (default 3).
.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
suppressMessages({library(Consensus); library(TreeTools)})

useGc <- identical(Sys.getenv("MEM_GCTORTURE"), "1")
nsEnv <- Sys.getenv("MEM_NS")
ns <- if (nzchar(nsEnv)) as.integer(strsplit(nsEnv, ",")[[1]]) else c(10L, 12L, 14L, 16L)
reps <- as.integer(Sys.getenv("MEM_REPS", "3"))
if (useGc) gctorture(TRUE)

# Adversarial multiPhylo batteries on n leaves, all sharing labels t1..tn.
batteries <- function(n) {
  cat0 <- ape::as.phylo(0L, n)
  star <- ape::stree(n, "star"); star[["tip.label"]] <- paste0("t", seq_len(n))
  list(
    identical = structure(list(cat0, cat0, cat0), class = "multiPhylo"),
    random    = ape::as.phylo(0:9, n),
    conflict  = ape::as.phylo(c(0L, 1L, 7L, 13L, 41L), n),
    star      = structure(list(star, star), class = "multiPhylo")
  )
}

count <- 0L
for (r in seq_len(reps)) {
  for (n in ns) {
    bs <- batteries(n)
    for (nm in names(bs)) {
      for (ty in c("rooted", "induced")) {
        res <- Local(bs[[nm]], ty)
        stopifnot(inherits(res, "phylo"), TreeTools::NTip(res) == n)
        count <- count + 1L
      }
    }
  }
}
cat(sprintf("OK: %d Local() calls, ns=%s, gctorture=%s\n",
            count, paste(ns, collapse = ","), useGc))
