# Verify Transfer() tip-set / duplicate-label validation against main source.
# TC-002: heterogeneous tip counts (subset) -> error or SILENT corruption?
# TC-003: duplicate tip labels in the first tree -> error or corruption?
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))

ok <- function(x) if (inherits(x, "try-error")) paste0("ERROR: ", conditionMessage(attr(x, "condition"))) else
  paste0("RETURNED phylo with ", NSplits(x), " splits, ", NTip(x), " tips")

t6a <- BalancedTree(letters[1:6]); t6b <- PectinateTree(letters[1:6])
t5  <- BalancedTree(letters[1:5])

cat("== sanity: homogeneous 6-tip ==\n")
cat(ok(try(Transfer(c(t6a, t6b)), silent = TRUE)), "\n\n")

cat("== TC-002a: trees[[1]] 6-tip, trees[[2]] 5-tip SUBSET ==\n")
het1 <- structure(list(t6a, t5), class = "multiPhylo")
cat(ok(try(Transfer(het1), silent = TRUE)), "\n\n")

cat("== TC-002b: trees[[1]] 5-tip, trees[[2]] 6-tip SUPERSET ==\n")
het2 <- structure(list(t5, t6a), class = "multiPhylo")
cat(ok(try(Transfer(het2), silent = TRUE)), "\n\n")

cat("== TC-003: duplicate label 'a' in trees[[1]] ==\n")
tdup <- t6a; tdup$tip.label <- c("a", "a", "b", "c", "d", "e")
cat(ok(try(Transfer(structure(list(tdup, t6a), class = "multiPhylo")), silent = TRUE)), "\n")
