# check-vs-legacy.R  (dev-only; run AFTER the rewrite + reinstall)
# Asserts the NEW fast RStar reproduces the OLD build's clade set EXACTLY over the
# deterministic legacy-grid (n = 5..200).  R* clade sets are unique (Jansson et
# al. 2016, Lemma 1.1), so new == old is an exact regression gate -- and it covers
# the medium-n, many-cluster regime the brute-force strong-cluster oracle
# (n <= 12) cannot reach.  Exits non-zero on any mismatch.
#
# Run from the package root:  Rscript.exe dev/oracle/rstar/check-vs-legacy.R

.libPaths(c(Sys.getenv("CONSTREE_LIB", "C:/Users/pjjg18/GitHub/Consensus/.agent-cons"), .libPaths()))
suppressMessages({ library(ConsTree); library(TreeTools); library(ape) })
source("dev/oracle/oracle.R")             # CladeSet()
source("dev/oracle/rstar/legacy-grid.R")  # legacyTrials()

# --- gate: confirm this is the NEW (cap-removed) build ------------------------
# Fingerprint: the shared .agent-cons is reinstalled by sibling chip sessions, so
# verify by BEHAVIOUR -- only this chip's build runs R* above 200 leaves.
ok <- tryCatch({
  t <- ape::rtree(210L, rooted = TRUE); t[["edge.length"]] <- NULL
  length(CladeSet(RStar(list(t, t)))) > 0L
}, error = function(e) FALSE)
if (!isTRUE(ok)) {
  stop("check-vs-legacy.R: installed build still caps at <= 200 leaves (version ",
       as.character(packageVersion("ConsTree")),
       "). The shared library was likely re-clobbered; reinstall the new build ",
       "and rerun. Aborting.")
}
cat("New build confirmed (version ", as.character(packageVersion("ConsTree")),
    "; R* runs at n > 200).\n", sep = "")

snap   <- readRDS("dev/oracle/rstar/legacy-clades.rds")
trials <- legacyTrials()
if (!setequal(names(snap), names(trials)))
  stop("legacy-grid drift: snapshot keys != current trial keys. Regenerate the snapshot.")

fail <- 0L
for (key in names(trials)) {
  newcl <- sort(CladeSet(RStar(trials[[key]])))
  if (!setequal(newcl, snap[[key]])) {
    fail <- fail + 1L
    cat("MISMATCH ", key, "\n", sep = "")
    cat("   extra in new:   ", paste(setdiff(newcl, snap[[key]]), collapse = " ; "), "\n")
    cat("   missing in new: ", paste(setdiff(snap[[key]], newcl), collapse = " ; "), "\n")
  }
}
cat(sprintf("\nnew-vs-old: %d / %d trials match.\n", length(trials) - fail, length(trials)))
if (fail > 0L) quit(status = 1L)
cat("ALL new == old (exact clade-set match to n = 200).\n")
