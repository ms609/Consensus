# check-local.R
# Cross-validation driver: Local(trees, type) vs LocalOracle(trees, type).
# Compares using CladeSet (rooted clade sets) for both "rooted" and "induced".
# Run from the package root:
#   Rscript.exe dev/oracle/local/check-local.R

.libPaths(c("C:/Users/pjjg18/GitHub/Consensus/.agent-cons", .libPaths()))
library(Consensus)
library(TreeTools)
source("dev/oracle/local/oracle_local.R")
source("dev/oracle/oracle.R")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
cladeSet <- function(tree) CladeSet(tree)

checkFixture <- function(label, trees, type) {
  mine   <- Local(trees, type)
  oracle <- tryCatch(LocalOracle(trees, type), error = function(e) NULL)

  if (is.null(oracle)) {
    cat(sprintf("  %-10s  %-8s  SKIPPED (oracle NULL)\n", label, type))
    return(invisible(NULL))
  }

  myClades  <- cladeSet(mine)
  orcClades <- cladeSet(oracle)
  ok        <- setequal(myClades, orcClades)
  status    <- if (ok) "MATCH" else "DIFFER"
  cat(sprintf("  %-10s  %-8s  %s\n", label, type, status))

  if (!ok) {
    extra   <- setdiff(myClades,  orcClades)
    missing <- setdiff(orcClades, myClades)
    if (length(extra))   cat("    Extra in mine:    ", paste(extra,   collapse="; "), "\n")
    if (length(missing)) cat("    Missing in mine:  ", paste(missing, collapse="; "), "\n")
  }
  ok
}

# ---------------------------------------------------------------------------
# Fixture battery
# ---------------------------------------------------------------------------
# Returns list(label, trees, types) ready for checking.
# multi2list: safely converts a multiPhylo to a plain list of phylo objects.
multi2list <- function(mp) lapply(seq_along(mp), function(i) mp[[i]])

fixture <- function(label, trees, types = c("rooted", "induced")) {
  list(label = label, trees = trees, types = types)
}

fixtures <- list(
  # --- smoke case (n=4) ---
  fixture("smoke-n4",
    list(ape::read.tree(text = "(1,((2,3),4));"),
         ape::read.tree(text = "(1,((2,4),3));"))),

  # --- identical trees: result should equal the tree itself ---
  fixture("ident-n4",
    list(ape::as.phylo(0, 4), ape::as.phylo(0, 4), ape::as.phylo(0, 4))),

  fixture("ident-n6",
    list(ape::as.phylo(3, 6), ape::as.phylo(3, 6), ape::as.phylo(3, 6))),

  # --- small n (n=5, various k) ---
  fixture("n5-k3",    multi2list(ape::as.phylo(0:2, 5))),
  fixture("n5-k5",    multi2list(ape::as.phylo(0:4, 5))),
  fixture("n5-k10",   multi2list(ape::as.phylo(0:9, 5))),

  # --- n=6 ---
  fixture("n6-k5",    multi2list(ape::as.phylo(0:4, 6))),
  fixture("n6-k10",   multi2list(ape::as.phylo(0:9, 6))),

  # --- conflict cases n=8 ---
  fixture("n8-conf1", multi2list(ape::as.phylo(c(0,0,0,1,2,53,99), 8))),
  fixture("n8-k4",    multi2list(ape::as.phylo(0:3, 8))),
  fixture("n8-k8",    multi2list(ape::as.phylo(0:7, 8))),

  # --- n=10 ---
  fixture("n10-k5",   multi2list(ape::as.phylo(0:4, 10))),
  fixture("n10-k10",  multi2list(ape::as.phylo(0:9, 10))),

  # --- n=12 (low conflict — only a few trees) ---
  fixture("n12-k3",   multi2list(ape::as.phylo(0:2, 12))),
  fixture("n12-k5",   multi2list(ape::as.phylo(0:4, 12))),

  # --- n=14 (low conflict) ---
  fixture("n14-k3",   multi2list(ape::as.phylo(0:2, 14))),
  fixture("n14-k5",   multi2list(ape::as.phylo(0:4, 14))),

  # --- n=16 (low conflict only — high conflict would be too slow) ---
  fixture("n16-k3",   multi2list(ape::as.phylo(0:2, 16))),

  # --- discriminating fixture: rooted != induced (seed=42, trial 52) ---
  # Confirms the two modes genuinely diverge and each matches its oracle.
  fixture("disc-n8k4", {
    set.seed(42)
    trs <- NULL
    for (trial in 1:100) {
      n <- sample(5:8, 1); k <- sample(3:8, 1)
      trees <- lapply(1:k, function(i) ape::rtree(n, rooted=TRUE))
      labs <- trees[[1]]$tip.label
      trees <- lapply(trees, function(tr) TreeTools::RenumberTips(tr, labs))
      cr <- tryCatch(CladeSet(Local(trees, "rooted")), error = function(e) NULL)
      ci <- tryCatch(CladeSet(Local(trees, "induced")), error = function(e) NULL)
      if (!is.null(cr) && !is.null(ci) && !setequal(cr, ci)) { trs <- trees; break }
    }
    trs
  })
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
cat("====  Local consensus oracle cross-check  ====\n")
cat(sprintf("  %-10s  %-8s  %s\n", "Fixture", "Type", "Result"))
cat(sprintf("  %-10s  %-8s  %s\n", "-------", "----",  "------"))

results <- list()
for (fx in fixtures) {
  for (tp in fx$types) {
    key <- paste0(fx$label, ":", tp)
    out <- tryCatch(
      checkFixture(fx$label, fx$trees, tp),
      error = function(e) {
        cat(sprintf("  %-10s  %-8s  ERROR: %s\n", fx$label, tp, conditionMessage(e)))
        FALSE
      }
    )
    results[[key]] <- out
  }
}

cat("\n====  Summary  ====\n")
outcomes <- unlist(results[!sapply(results, is.null)])
nPass   <- sum(outcomes, na.rm = TRUE)
nFail   <- sum(!outcomes, na.rm = TRUE)
nSkip   <- sum(sapply(results, is.null))
cat(sprintf("PASS: %d   FAIL: %d   SKIP: %d\n", nPass, nFail, nSkip))
if (nFail > 0) {
  cat("FAILING fixtures:\n")
  for (nm in names(results)) {
    if (!is.null(results[[nm]]) && !results[[nm]]) cat(" ", nm, "\n")
  }
  quit(status = 1)
} else {
  cat("ALL non-NULL fixtures MATCH.\n")
}
