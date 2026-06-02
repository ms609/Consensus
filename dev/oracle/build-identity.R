# Build-identity guard for oracle and profiling scripts.
# Source this file after library(ConsTree) and call assertConsTreeBuild().
#
# Three layers of defence against a stale install:
#   1. Version check  -- installed packageVersion must match DESCRIPTION.
#   2. Body check     -- each fast-path R wrapper must reference its C++ symbol
#                        (catches a stale R source, e.g. pre-port wrapper copied
#                        in place of the updated one).
#   3. Idempotence    -- Method(t, t, t) must reproduce t's split count
#                        (catches a stale DLL where version and body both pass
#                        but the compiled algorithm is wrong or a stub).
#
# Extend .FAST_PATH_CHECKS when a new C++ fast path is merged to main (i.e.
# when the R wrapper for a method is updated to call its *Cpp function).

.REBUILD_MSG <- paste0(
  "Force a clean rebuild:\n",
  "  # PowerShell (if DLL is locked, clear stale objects first):\n",
  "  Remove-Item src\\*.o, src\\*.dll -ErrorAction SilentlyContinue\n",
  "  R.exe CMD INSTALL --preclean --library=.agent-cons ."
)

# Methods with an active C++ fast path in the installed R wrapper.
.FAST_PATH_CHECKS <- c(
  Loose        = "looseConsensusCpp",
  Greedy       = "greedyConsensusCpp",
  MajorityPlus = "majorityPlusConsensusCpp",
  Frequency    = "frequencyConsensusCpp",
  RStar        = "rStarConsensus",
  Local        = "localConsensus"
)

# Fixed 12-tip rooted binary tree for idempotence probes.  n=12 is well inside
# Local's 20-leaf cap and RStar's 200-leaf cap.  Uses ape::read.tree (no RNG,
# no set.seed side-effects).  NSplits should equal 9 (fully binary, n-3).
.GUARD_TREE <- ape::read.tree(
  text = "(((1,2),(3,4)),((5,6),(7,8)),((9,10),(11,12)));"
)

assertConsTreeBuild <- function(pkgRoot = normalizePath(getwd())) {
  ns <- asNamespace("ConsTree")

  # --- 0. Print identity ---------------------------------------------------
  ver <- as.character(utils::packageVersion("ConsTree"))
  loc <- dirname(system.file(package = "ConsTree"))
  cat(sprintf("ConsTree %s  [%s]\n", ver, loc))

  # --- 1. Version must match DESCRIPTION ------------------------------------
  desc_path <- file.path(pkgRoot, "DESCRIPTION")
  if (file.exists(desc_path)) {
    lines    <- readLines(desc_path, warn = FALSE)
    m        <- regmatches(lines, regexpr("^Version:\\s*\\S+", lines))
    if (length(m)) {
      desc_ver <- trimws(sub("^Version:\\s*", "", m[[1L]]))
      if (ver != desc_ver) {
        stop(sprintf(
          "Stale build: installed ConsTree %s != DESCRIPTION %s.\n%s",
          ver, desc_ver, .REBUILD_MSG))
      }
    }
  }

  # --- 2. DLL is linked and callable ----------------------------------------
  sc <- tryCatch(
    get("consensus_rcpp_selfcheck", envir = ns, inherits = FALSE)(),
    error = function(e) NULL
  )
  if (!identical(sc, 42L)) {
    stop(sprintf(
      "DLL selfcheck failed (expected 42, got %s).\n%s",
      if (is.null(sc)) "ERROR (symbol not found)" else sc,
      .REBUILD_MSG))
  }

  # --- 3. Fast-path R wrappers reference their C++ symbol -------------------
  for (fn in names(.FAST_PATH_CHECKS)) {
    sym      <- .FAST_PATH_CHECKS[[fn]]
    fn_body  <- paste(deparse(body(get(fn, envir = ns, inherits = FALSE))),
                      collapse = "\n")
    if (!grepl(sym, fn_body, fixed = TRUE)) {
      stop(sprintf(
        "Stale build: %s() body does not call %s().\n%s",
        fn, sym, .REBUILD_MSG))
    }
  }

  # --- 4. Functional idempotence: Method(t, t, t) must reproduce t ----------
  # A stub or wrong algorithm returns a star (0 splits) or errors; the correct
  # algorithm returns t (9 splits for .GUARD_TREE).
  t      <- .GUARD_TREE
  trees3 <- structure(list(t, t, t), class = "multiPhylo")
  expected <- TreeTools::NSplits(t)
  for (fn in names(.FAST_PATH_CHECKS)) {
    result <- tryCatch(
      get(fn, envir = ns, inherits = FALSE)(trees3),
      error = function(e) NULL
    )
    if (is.null(result)) {
      stop(sprintf("Stale build: %s(t,t,t) errored.\n%s", fn, .REBUILD_MSG))
    }
    got <- TreeTools::NSplits(result)
    if (!identical(got, expected)) {
      stop(sprintf(
        "Stale build: %s(t,t,t) returned %d splits, expected %d.\n%s",
        fn, got, expected, .REBUILD_MSG))
    }
  }

  cat("Build identity OK.\n")
  invisible(NULL)
}
