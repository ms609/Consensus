# Development-only bridge to the freqdiff oracle binary
# (dev/oracle/freqdiff/freqdiff.exe).
# NOT part of the package: used only to cross-validate ConsTree::Frequency()
# against the FDCT_new reference implementation by Jesper Jansson et al.
#
# Binary I/O (reverse-engineered from dev/reference/FDCT_new/main_new.cpp,
# nex_parser.h, and Tree.cpp):
#   Input:  inp.txt in the process CWD — one Newick per line, integer taxa
#           labels only (Tree.cpp asserts isdigit(*str) for leaf tokens).
#           No semicolons or other decoration required; bare "(1,(2,3))"
#           format works. Blank / short (<=2 char) lines terminate parsing.
#   Output: oup.txt in the process CWD — a single Newick line with the
#           original integer label strings. No trailing semicolon.
#   Stdout: two timing floats (from calc_w_knlogn + freqdiff); ignored here.
#
# IMPORTANT cwd dependency: both inp.txt and oup.txt are relative paths hard-
# coded in the binary (LOCAL_TEST block in main_new.cpp). The bridge sets
# working directory to a temporary scratch dir before invoking the binary,
# then restores the caller's cwd.

FREQDIFF_EXE <- normalizePath(
  "C:/Users/pjjg18/GitHub/Consensus/dev/oracle/freqdiff/freqdiff.exe",
  mustWork = TRUE
)

# Write trees as integer-labelled Newicks (one per line) to a file.
.WriteFreqdiffInput <- function(trees, path, labels) {
  n <- length(labels)
  treeLines <- vapply(trees, function(tr) {
    # Renumber tips to labels[1..n] order, then relabel with integers 1..n so
    # taxon i carries label "i".  Root every tree at taxon 1 ("1") BEFORE writing
    # -- frequency-difference weights ROOTED clusters, so both sides must feed
    # identically rooted trees (the C++ path roots at taxon 1 via .FactEdges).
    # Relabel before rooting so the integer<->label mapping survives any tip
    # reordering RootTree performs (the decode is labels[as.integer(label)]).
    tr <- TreeTools::RenumberTips(tr, labels)
    tr[["tip.label"]] <- as.character(seq_len(n))
    tr <- TreeTools::RootTree(tr, "1")
    tr[["edge.length"]] <- NULL
    ape::write.tree(tr)  # returns string; file=NULL is unreliable in some ape versions
  }, character(1))
  # Remove trailing semicolons — Tree.cpp doesn't need them and the assertion
  # on the first non-'(' character only tolerates digits.
  treeLines <- sub(";$", "", treeLines)
  writeLines(treeLines, path)
}

# Run the freqdiff binary for the frequency-difference consensus, returning
# the consensus as a `phylo` with the ORIGINAL tip labels restored.
FreqDiffOracle <- function(trees) {
  labels <- TreeTools::TipLabels(trees[[1]])

  # Use a per-call temp dir so concurrent calls don't stomp on inp/oup.txt
  scratch <- file.path(
    "C:/Users/pjjg18/GitHub/Consensus/dev/oracle/freqdiff",
    paste0(".scratch_", Sys.getpid())
  )
  dir.create(scratch, showWarnings = FALSE)
  on.exit({
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)

  inpPath <- file.path(scratch, "inp.txt")
  oupPath <- file.path(scratch, "oup.txt")

  .WriteFreqdiffInput(trees, inpPath, labels)

  # The binary reads/writes inp.txt and oup.txt relative to its cwd.
  oldwd <- setwd(scratch)
  on.exit(setwd(oldwd), add = TRUE)

  ret <- system2(FREQDIFF_EXE, stdout = TRUE, stderr = FALSE)

  if (!file.exists(oupPath)) {
    stop("freqdiff binary did not produce oup.txt.\nstdout:\n",
         paste(ret, collapse = "\n"))
  }

  nwk <- readLines(oupPath, warn = FALSE)
  nwk <- nwk[nzchar(nwk)]
  if (length(nwk) == 0L) {
    stop("oup.txt is empty — binary may have produced no valid consensus.")
  }
  nwk <- nwk[[1L]]
  nwk <- paste0(nwk, ";")  # ape::read.tree requires a trailing semicolon

  tr <- ape::read.tree(text = nwk)
  # tip.label values are the integer strings the binary echoed; map back.
  tr[["tip.label"]] <- labels[as.integer(tr[["tip.label"]])]
  # Return:
  tr
}

# Canonical, order-independent split-string set for comparison.
# Replicates the helper in dev/oracle/oracle.R.
SplitSetFD <- function(tree, labels) {
  if (TreeTools::NSplits(tree) == 0L) {
    character(0)
  } else {
    unname(as.character(TreeTools::as.Splits(tree, tipLabels = labels)))
  }
}
