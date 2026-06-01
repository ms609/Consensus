# Net concordant quartet information optimized by QuartetConsensus():
# summed over input trees, (same - penalty * different) quartet counts.
# The "abstain" convention means quartets unresolved in an input tree
# contribute to neither term, so this reads directly off QuartetStatus().
QcObjective <- function(consensus, trees, penalty = 1) {
  status <- Quartet::QuartetStatus(trees, cf = consensus)
  if (is.null(dim(status))) status <- t(status)
  sum(status[, "s"]) - penalty * sum(status[, "d"])
}


test_that("penalty must be a single non-negative number", {
  library("TreeTools")
  trees <- as.phylo(1:6, nTip = 8)
  expect_error(QuartetConsensus(trees, penalty = -1), "non-negative")
  expect_error(QuartetConsensus(trees, penalty = c(1, 2)), "single")
  expect_error(QuartetConsensus(trees, penalty = NA_real_), "non-negative")
  expect_error(QuartetConsensus(trees, penalty = "1"), "non-negative")
})


test_that("star tree scores zero; consensus beats star and majority rule", {
  skip_if_not_installed("Quartet")
  library("TreeTools")
  trees <- as.phylo(1:25, nTip = 8)

  qc <- QuartetConsensus(trees)
  star <- StarTree(TipLabels(trees[[1]]))
  mr <- UnrootTree(Consensus(trees, p = 0.5))

  # The star resolves no quartets, so the absolute objective is exactly zero:
  # this is the baseline that makes the score calibrated across leaf counts.
  expect_equal(QcObjective(star, trees), 0)

  # The greedy consensus starts from majority rule and only makes improving
  # moves, so it must beat the star and do no worse than majority rule.
  expect_gt(QcObjective(qc, trees), 0)
  expect_gte(QcObjective(qc, trees), QcObjective(mr, trees))
})


test_that("penalty = 1 on binary inputs is the quartet-distance median", {
  skip_if_not_installed("Quartet")
  library("TreeTools")
  # Maximizing (s - d) with binary inputs is equivalent to minimizing the
  # symmetric quartet distance 2d + r1 + r2, so the consensus is no worse
  # than majority rule under that distance.
  trees <- as.phylo(1:20, nTip = 8)

  qc <- QuartetConsensus(trees, penalty = 1)
  Qd <- function(tree) {
    sum(vapply(trees, function(tr) {
      s <- Quartet::QuartetStatus(tree, tr)
      2 * s[, "d"] + s[, "r1"] + s[, "r2"]
    }, numeric(1)))
  }
  expect_lte(Qd(qc), Qd(Consensus(trees, p = 0.5)))
})


test_that("penalty controls resolution: low resolves more than high", {
  library("TreeTools")
  set.seed(91)
  trees <- as.phylo(sample.int(945, 30), nTip = 9)

  # penalty 0 resolves every quartet it can; a large penalty resolves only
  # near-unanimous ones.
  qc_lenient <- QuartetConsensus(trees, penalty = 0)
  qc_strict <- QuartetConsensus(trees, penalty = 50)

  expect_gte(NSplits(qc_lenient), NSplits(qc_strict))
})


test_that("abstain: a rare but consistent resolution is recovered", {
  skip_if_not_installed("Quartet")
  library("TreeTools")
  resolved <- BalancedTree(8)
  star <- StarTree(TipLabels(resolved))
  # 7 uninformative (star) trees + 3 that agree on `resolved`: only 30% of
  # trees resolve anything, but they agree, so under the abstain convention
  # (input polytomies do not vote) the structure is recovered.
  trees <- structure(c(rep(list(star), 7), rep(list(resolved), 3)),
                     class = "multiPhylo")

  qc <- QuartetConsensus(trees, penalty = 1)

  expect_gt(NSplits(qc), 0)
  expect_equal(unname(Quartet::QuartetStatus(qc, resolved)[, "d"]), 0)
})


test_that("clean data is never pruned and stays fully resolved", {
  skip_if_not_installed("Quartet")
  library("TreeTools")
  tree <- BalancedTree(8)
  trees <- structure(rep(list(tree), 10), class = "multiPhylo")

  qc <- QuartetConsensus(trees, neverDrop = FALSE)

  expect_equal(NTip(qc), 8)
  expect_equal(length(attr(qc, "dropped")), 0)
  expect_equal(unname(Quartet::QuartetStatus(qc, tree)[, "d"]), 0)
})


test_that("drop_scores record the strictly increasing objective", {
  library("TreeTools")
  base <- BalancedTree(8)
  intermediate <- AddTipEverywhere(base, "rogue1")
  trees <- do.call(c, lapply(intermediate, function(tr) {
    AddTipEverywhere(tr, "rogue2")
  }))
  class(trees) <- "multiPhylo"

  qc <- QuartetConsensus(trees, neverDrop = FALSE)
  scores <- attr(qc, "drop_scores")
  dropped <- attr(qc, "dropped")

  expect_equal(length(scores), length(dropped))
  expect_true(all(is.finite(scores)))
  # Every greedy move (including a drop) strictly improves the objective.
  if (length(scores) > 1L) {
    expect_true(all(diff(scores) > 0))
  }
})
