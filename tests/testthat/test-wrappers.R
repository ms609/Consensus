test_that("Strict() matches TreeTools::Consensus(p = 1)", {
  trees <- ape::as.phylo(0:5, 8)
  expect_equal(Strict(trees), TreeTools::Consensus(trees, p = 1))
})

test_that("Majority() matches TreeTools::Consensus(p = 0.5)", {
  trees <- ape::as.phylo(0:5, 8)
  expect_equal(Majority(trees), TreeTools::Consensus(trees, p = 0.5))
  expect_equal(MajorityRule(trees), Majority(trees))
  expect_equal(Majority(trees, p = 1), Strict(trees))
})

test_that("A single tree is its own consensus", {
  tree <- ape::as.phylo(0, 8)
  expect_equal(Strict(list(tree)), tree)
  expect_equal(Majority(list(tree)), tree)
})
