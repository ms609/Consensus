# Cross-validate the R consensus methods against the reference FACT binary.
# Run with: Rscript dev/oracle/check-oracle.R
# Resolve everything relative to THIS script so the oracle uses this worktree's
# OWN build, never a sibling worktree's.  Each worktree installs into its own
# .agent-cons (gitignored by the rooted /.agent* rule); worktrees never share a
# library, which removes the clobber/self-mask hazard at the source.  Install:
#   R CMD INSTALL --no-multiarch --library=<worktree>/.agent-cons <worktree>
.cohArgs <- commandArgs(FALSE)
.cohFile <- sub("^--file=", "", grep("^--file=", .cohArgs, value = TRUE))
.cohDir  <- if (length(.cohFile)) dirname(normalizePath(.cohFile)) else getwd()
.cohRoot <- normalizePath(file.path(.cohDir, "..", ".."))
.cohLib  <- file.path(.cohRoot, ".agent-cons")
if (!dir.exists(.cohLib)) {
  stop("Validation library not found: ", .cohLib,
       "\n  Install first:  R CMD INSTALL --no-multiarch --library=\"", .cohLib,
       "\" \"", .cohRoot, "\"")
}
.libPaths(c(.cohLib, .libPaths()))
suppressMessages(library(ConsTree))
suppressMessages(library(TreeTools))

# --- self-guard ------------------------------------------------------------
# Even with per-worktree libraries, guard against a stale/foreign build (e.g. an
# install that silently partial-failed under a DLL lock): the installed version
# must match THIS worktree's source, and (this branch's deliverable) MajorityPlus
# must be on the C++ path.  See the `agent-cons-install-can-silently-fail` note.
local({
  want <- read.dcf(file.path(.cohRoot, "DESCRIPTION"), fields = "Version")[1, 1]
  have <- as.character(utils::packageVersion("ConsTree"))
  if (have != want) {
    stop(sprintf(paste0("[self-guard] Installed ConsTree %s != this worktree's %s",
                        " -- reinstall this worktree (into its own .agent-cons)",
                        " before trusting the oracle."),
                 have, want))
  }
  if (!any(grepl("majorityPlusConsensusCpp", deparse(body(ConsTree::MajorityPlus))))) {
    stop("[self-guard] MajorityPlus is not on the C++ path -- installed build",
         " predates the port.")
  }
  cat(sprintf("[self-guard] ConsTree %s (this worktree), MajorityPlus on C++ path: OK\n",
              have))
})
# ---------------------------------------------------------------------------

source(file.path(.cohDir, "oracle.R"))

cmp <- function(mine, fact, labels) {
  setequal(SplitSet(mine, labels), SplitSet(fact, labels))
}

datasets <- list(
  "random  n9  k21" = ape::as.phylo(0:20, 9),
  "random  n10 k31" = ape::as.phylo(0:30, 10),
  "conflict n8  k7" = ape::as.phylo(c(0, 0, 0, 1, 2, 53, 99), 8)
)
methods <- list(strict = Strict, majority = Majority, greedy = Greedy,
                loose = Loose, majorityPlus = MajorityPlus)

trees1 <- datasets[[1]]
labels1 <- TipLabels(trees1[[1]])
cat("== Determine FACT rooted flag (strict) ==\n")
for (rt in c(0L, 1L)) {
  ok <- cmp(Strict(trees1), FactConsensus(trees1, "strict", rooted = rt), labels1)
  cat(sprintf("  strict vs FACT rooted=%d : %s\n", rt, if (ok) "MATCH" else "differ"))
}

for (rt in c(0L, 1L)) {
  cat(sprintf("\n== Cross-validation (rooted=%d) ==\n", rt))
  for (dn in names(datasets)) {
    trees <- datasets[[dn]]
    labels <- TipLabels(trees[[1]])
    cat("--", dn, "--\n")
    for (mn in names(methods)) {
      mine <- methods[[mn]](trees)
      fact <- FactConsensus(trees, mn, rooted = rt)
      ok <- cmp(mine, fact, labels)
      cat(sprintf("  %-13s mine=%2d fact=%2d  %s\n",
                  mn, NSplits(mine), NSplits(fact),
                  if (ok) "MATCH" else "*** DIFFER ***"))
    }
  }
}

# Adams is a ROOTED method: validate against the classical (slow) Adams with
# rooted = 1 (each input tree's own root), comparing ROOTED CLADES, not splits.
cat("\n== Adams cross-validation (rooted=1, clade comparison) ==\n")
for (dn in names(datasets)) {
  trees <- datasets[[dn]]
  mine <- Adams(trees)
  fact <- FactConsensus(trees, "adams", rooted = 1L)
  cm <- CladeSet(mine)
  cf <- CladeSet(fact)
  cat(sprintf("  %-16s mine=%2d fact=%2d  %s\n", dn, length(cm), length(cf),
              if (setequal(cm, cf)) "MATCH" else "*** DIFFER ***"))
}

# Multi-word bitset path (n > 60: BUCKET_SIZE = 60, so LEN > 1 -- the word-index
# arithmetic (c-1)/60, %60 and the LEN-length OR-up/compare loops).  The datasets
# above all have LEN = 1, so the multi-word packing needs its own check.
#
# Greedy is tie-break sensitive, and FACT's tie-break depends on its exact
# internal clade representation (it roots at taxon 1's neighbour; the R wrapper
# roots on taxon 1's edge), so an exact match on tie-heavy random input is NOT
# expected -- greedy is "FACT-match up to tie-break" (as the previous R
# implementation also was; AGENTS.md).  Validate the multi-word path two ways
# that ARE exact: (a) idempotence Greedy(list(t,t,t)) == t exercises the pack/
# unpack with no ties; (b) with every tree identically rooted at taxon 1, mine
# and FACT see the same clades, isolating algorithm faithfulness from the
# rooting/tie-break artefact.
cat("\n== Multi-word bitset path (n > 60) ==\n")
for (n in c(80L, 137L)) {
  set.seed(n)
  labs <- paste0("t", seq_len(n))
  base <- RandomTree(labs, root = TRUE)
  idem <- setequal(SplitSet(Greedy(structure(list(base, base, base),
                                             class = "multiPhylo")), labs),
                   SplitSet(base, labs))
  trees <- structure(lapply(1:15, function(i)
                       RootTree(RandomTree(labs, root = TRUE), labs[[1]])),
                     class = "multiPhylo")
  ok <- cmp(Greedy(trees), FactConsensus(trees, "greedy", rooted = 1L), labs)
  cat(sprintf("  n=%-3d LEN=%d  idempotent: %-5s   FACT-exact (same rooting): %s\n",
              n, (n + 59L) %/% 60L, idem,
              if (ok) "MATCH" else "*** DIFFER ***"))
}

# Large-n majorityPlus.  Unlike Greedy, majorityPlus does NOT bit-pack (there is
# no BUCKET_SIZE), so "multi-word" here is a misnomer: n > 60 instead stresses
# the Day's leaf-relabelling / left-right path-query machinery at scale.  Because
# majorityPlus is a deterministic count rule -- keep a clade iff it is displayed
# by strictly more trees than contradict it, with NO frequency tie-break -- it
# must be FACT-EXACT at every n; divergence is a real bug, not a tie-break
# artefact.  Independent random trees would collapse to a star (essentially every
# non-trivial split is contradicted more often than displayed), giving a vacuous
# star-vs-star match, so drive it with CONGRUENT input (one base topology, each
# replicate perturbed by a few tip swaps) and assert a non-trivial result.
cat("\n== majorityPlus at n > 60 (exact) ==\n")
mpPass <- logical(0)
for (n in c(80L, 137L)) {
  set.seed(n + 1000L)
  labs <- paste0("t", seq_len(n))
  base <- RootTree(RandomTree(labs, root = TRUE), labs[[1]])
  idem <- setequal(SplitSet(MajorityPlus(structure(list(base, base, base),
                                                   class = "multiPhylo")), labs),
                   SplitSet(base, labs))
  swap <- function(tr) {
    for (s in 1:3) {
      ij <- sample.int(n, 2L)
      tr[["tip.label"]][ij] <- tr[["tip.label"]][rev(ij)]
    }
    RootTree(tr, labs[[1]])
  }
  trees <- structure(c(list(base), lapply(1:14, function(i) swap(base))),
                     class = "multiPhylo")
  mine <- MajorityPlus(trees)
  ok <- cmp(mine, FactConsensus(trees, "majorityPlus", rooted = 1L), labs)
  mpPass <- c(mpPass, idem, ok, NSplits(mine) > 0L)
  cat(sprintf("  n=%-3d LEN=%d  idempotent: %-5s  splits=%2d  FACT-exact: %s\n",
              n, (n + 59L) %/% 60L, idem, NSplits(mine),
              if (ok) "MATCH" else "*** DIFFER ***"))
}
# Hard assertion: idempotent, FACT-exact, and non-trivial (not a star) at all n.
stopifnot(all(mpPass))
