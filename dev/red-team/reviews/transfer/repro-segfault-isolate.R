# Localise the Transfer() segfault on normal homogeneous input.
# The last line printed before the process dies pinpoints the crashing step.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))
f <- function(s) { cat(s, "\n"); flush.console() }
f("STEP loaded")
t6a <- BalancedTree(letters[1:6]); f("STEP t6a")
t6b <- PectinateTree(letters[1:6]); f("STEP t6b")
m <- structure(list(t6a, t6b), class = "multiPhylo"); f(paste("STEP multiPhylo len", length(m)))
# Inspect the R-side split pooling that feeds C++, before the C++ call.
sl <- lapply(m, function(tr) TreeTools::as.Splits(tr, tipLabels = TreeTools::TipLabels(t6a)))
f(paste("STEP as.Splits ok; ncols",
        paste(vapply(sl, ncol, integer(1)), collapse = ",")))
f("STEP calling Transfer ...")
res <- Transfer(m)
f("STEP Transfer returned")
f(paste("STEP class", paste(class(res), collapse = ",")))
f(paste("STEP NSplits", TreeTools::NSplits(res)))
print(res)
f("STEP DONE")
