# Consensus

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/ms609/Consensus/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/Consensus/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/Consensus/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/Consensus)
<!-- badges: end -->

'Consensus' is an R package providing a comprehensive, efficient suite of
methods for summarizing a collection of phylogenetic trees — for example a
bootstrap or Bayesian posterior sample — as a single **consensus tree**.

It builds on the tree and split infrastructure of
['TreeTools'](https://ms609.github.io/TreeTools/), and repackages a family of
asymptotically-efficient consensus algorithms developed by Jesper Jansson and
colleagues, several of which have not previously been available in R.

## Consensus methods

Each method takes a list of trees (or a `multiPhylo`) sharing the same leaves and
returns a rooted `phylo`. They differ in *which* groupings (splits) they retain:

| Function | Retains a grouping when… |
|----------|--------------------------|
| `Strict()` | it occurs in **every** tree |
| `Majority()` / `MajorityRule()` | it occurs in **> half** the trees (tunable via `p`) |
| `Loose()` | **no** tree contradicts it (semi-strict / combinable-component) |
| `MajorityPlus()` | **more** trees display it than contradict it |
| `Frequency()` | it is **more frequent than every** grouping that conflicts with it (frequency-difference) |
| `Greedy()` | added greedily, most frequent first, if compatible with those already kept (extended majority-rule) |
| `Average()` | (distance-based) the tree best fitting the mean path-length distances |

Forthcoming: `Adams()`, `Local()` (minimum rooted/induced local consensus), and
`RStar()`.

These fill genuine gaps in the R ecosystem: `ape::consensus()` offers only strict
and majority; `phangorn` partially overlaps on greedy (`allCompat`); loose,
majority-rule (+), and especially the **frequency-difference** consensus are not
otherwise readily available in R.

## Usage

```r
library("Consensus")

trees <- ape::as.phylo(1:100, 8)   # 100 eight-leaf trees

Strict(trees)        # most conservative
Majority(trees)      # the familiar 50% majority-rule tree
Loose(trees)         # everything not actively contradicted
Frequency(trees)     # frequency-difference: often more resolved than majority
Greedy(trees)        # most resolved of the split-based summaries
```

## Installation

Install the development version from GitHub:

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("ms609/Consensus")
```

'Consensus' is not yet on CRAN.

## Relationship to other packages

'Consensus' is the front-end consensus toolkit for the
[TreeTools](https://ms609.github.io/TreeTools/) ecosystem, which also includes
['TreeDist'](https://ms609.github.io/TreeDist/) (tree distances and
information-theoretic consensus) and
['TreeSearch'](https://ms609.github.io/TreeSearch/) (phylogenetic search).
TreeTools itself remains the fast engine for the strict and majority-rule
consensus, which 'Consensus' exposes through consistently-named wrappers.

## Citation and attribution

The algorithms repackaged here originate in the FACT, FACT2 and FDCT prototypes
of Jesper Jansson and colleagues, incorporated with permission; see
`inst/REFERENCES.bib` for the source papers, and `DESCRIPTION` for copyright
attribution. 'Consensus' is released under GPL (≥ 3).

Please note that this project is released with a
[Contributor Code of Conduct](https://ms609.github.io/TreeTools/CODE_OF_CONDUCT.html).
By contributing, you agree to abide by its terms.
