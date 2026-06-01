# Development-only bridge to the LOCAL-consensus oracle binary
# (dev/oracle/local/local.exe), a patched build of FDCT_new/main_new.cpp that
# takes the algorithm from argv[1] (minrlc_exact / minilc_exact / aho-build /
# freq) while keeping the inp.txt / oup.txt I/O, so this reuses the freqdiff
# bridge's machinery verbatim.  NOT part of the package.
#
# Reference: minRLC_exact / minILC_exact in
# dev/reference/FDCT_new/local_consensus.h (minRILC with minrs = TRUE / FALSE).
# Returns NULL when the binary reports "No valid consensus found." (minRILC
# returns NULL when the leaf set is a single inseparable Aho component).

LOCAL_EXE <- normalizePath(
  "C:/Users/pjjg18/GitHub/Consensus/dev/oracle/local/local.exe",
  mustWork = TRUE
)

.LocalAlgo <- c(rooted = "minrlc_exact", induced = "minilc_exact",
                aho = "aho-build")

.WriteLocalInput <- function(trees, path, labels) {
  n <- length(labels)
  treeLines <- vapply(trees, function(tr) {
    tr <- TreeTools::RenumberTips(tr, labels)
    tr[["tip.label"]] <- as.character(seq_len(n))
    tr[["edge.length"]] <- NULL
    sub(";$", "", ape::write.tree(tr))
  }, character(1))
  writeLines(treeLines, path)
}

# Run the local-consensus binary; returns a `phylo` with original labels, or
# NULL when the binary finds no valid consensus.
LocalOracle <- function(trees, type = c("rooted", "induced", "aho")) {
  type <- match.arg(type)
  labels <- TreeTools::TipLabels(trees[[1]])

  scratch <- file.path("C:/Users/pjjg18/GitHub/Consensus/dev/oracle/local",
                       paste0(".scratch_", Sys.getpid()))
  dir.create(scratch, showWarnings = FALSE)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  .WriteLocalInput(trees, file.path(scratch, "inp.txt"), labels)
  oupPath <- file.path(scratch, "oup.txt")

  oldwd <- setwd(scratch)
  on.exit(setwd(oldwd), add = TRUE)
  out <- system2(LOCAL_EXE, args = .LocalAlgo[[type]],
                 stdout = TRUE, stderr = FALSE)

  if (any(grepl("No valid consensus", out)) || !file.exists(oupPath)) {
    return(NULL)
  }
  nwk <- readLines(oupPath, warn = FALSE)
  nwk <- nwk[nzchar(nwk)]
  if (length(nwk) == 0L) {
    return(NULL)
  }
  tr <- ape::read.tree(text = paste0(nwk[[1L]], ";"))
  tr[["tip.label"]] <- labels[as.integer(tr[["tip.label"]])]
  # Return:
  tr
}
