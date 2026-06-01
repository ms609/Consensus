# Development-only bridge to the reference FACT binary (dev/oracle/fact.exe).
# NOT part of the package: it shells out to a compiled binary and is used only
# to cross-validate the R implementations against the original algorithms.
#
# FACT I/O (see dev/reference/FACT/src/{wrapper,tree}.cpp):
#  * Input: a NEXUS file whose Newick trees use integer taxon labels 1..n,
#    in a picky dialect (two header lines, a literal `translate` line, taxon
#    lines ending with `;`, `tree i = (...);` lines, then `END;`).
#  * Invocation: the filename, an algorithm bitmask `rule`, and a `rooted` flag
#    (0/1) are read from stdin.
#  * Output: the consensus printed by printNex() as an integer-labelled Newick
#    on its own line.

FACT_EXE <- normalizePath("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/fact.exe",
                          mustWork = TRUE)

# Algorithm bitmasks: bit i selects getConsensus(i) (main.cpp switch).
#  0 strict | 1 majSlow 2 majFast 3 majBest | 4 greedySlow 5 greedyFast |
#  6 looseSlow 7 looseFast | 8 majorityPlus | 9 adamsSlow 10 adamsFast
# Fast variants chosen, EXCEPT Adams: the fast (1024) variant mis-prints small
# trees (drops a leaf), so the reliable oracle is the classical slow Adams (512).
FACT_RULE <- c(strict = 1L, majority = 8L, greedy = 32L,
               loose = 128L, majorityPlus = 256L,
               adams = 512L, adamsFast = 1024L)

.WriteFactNexus <- function(trees, path, labels) {
  n <- length(labels)
  header <- c("#NEXUS", "[consensus oracle]", "translate")
  taxa <- c(paste0(seq_len(n - 1L), " t", seq_len(n - 1L), ","),
            paste0(n, " t", n, ";"))
  treeLines <- vapply(seq_along(trees), function(i) {
    tr <- TreeTools::RenumberTips(trees[[i]], labels)
    tr[["tip.label"]] <- as.character(seq_len(n))
    tr[["edge.length"]] <- NULL
    paste0("tree ", i, " = ", ape::write.tree(tr))
  }, character(1))
  writeLines(c(header, taxa, treeLines, "END;"), path)
}

# Run FACT for one method, returning the consensus as a `phylo` with the
# original leaf labels restored.
FactConsensus <- function(trees, method = "strict", rooted = 0L) {
  labels <- TreeTools::TipLabels(trees[[1]])
  nex <- file.path(dirname(FACT_EXE), "_oracle_input.nex")
  .WriteFactNexus(trees, nex, labels)
  out <- system2(FACT_EXE, input = paste(nex, FACT_RULE[[method]], rooted),
                 stdout = TRUE, stderr = FALSE)
  nwk <- grep("^\\(.*\\);$", out, value = TRUE)
  if (length(nwk) != 1L) {
    stop("Expected one Newick line from FACT; got ", length(nwk), ".\nOutput:\n",
         paste(out, collapse = "\n"))
  }
  tr <- ape::read.tree(text = nwk)
  tr[["tip.label"]] <- labels[as.integer(tr[["tip.label"]])]
  # Return:
  tr
}

# Canonical, order-independent split-string set for comparison (UNROOTED
# bipartitions) -- the right comparison for the split-based methods.
SplitSet <- function(tree, labels) {
  if (TreeTools::NSplits(tree) == 0L) {
    character(0)
  } else {
    unname(as.character(TreeTools::as.Splits(tree, tipLabels = labels)))
  }
}

# Canonical, order-independent set of ROOTED clades (sorted descendant
# tip-LABEL sets, size 2..n-1).  Necessary for rooted methods such as Adams,
# where split equality is insufficient (two trees can share all unrooted splits
# yet be rooted differently).  Uses each tree's OWN tip labels, so it is correct
# regardless of how the two trees order their tips.
CladeSet <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree)
  tl <- tree[["tip.label"]]
  clades <- vapply(seq_len(nrow(edge)), function(e) {
    below <- TreeTools::DescendantEdges(edge[, 1], edge[, 2], edge = e)
    tips <- edge[below, 2]
    paste(sort(tl[tips[tips <= nTip]]), collapse = ",")
  }, character(1))
  sizes <- lengths(strsplit(clades, ","))
  sort(unique(clades[sizes > 1 & sizes < nTip]))
}
