# ConsTree

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/ms609/ConsTree/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/ConsTree/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/ConsTree/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/ConsTree)
<!-- badges: end -->

'ConsTree' is an R package providing a comprehensive, efficient suite of
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
| `Adams()` | constructed from the finest root-level partition shared by **every** tree (may introduce novel groupings; rooted) |
| `Local()` | based on rooted triplets shared by **every** tree (minimum rooted/induced local consensus; ≤ 20 leaves) |
| `RStar()` | each rooted triplet grouping that wins a **plurality** against each alternative separately |

## Usage

```r
library("ConsTree")

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
remotes::install_github("ms609/ConsTree")
```

'ConsTree' is not yet on CRAN.

## Relationship to other packages

'ConsTree' is the front-end consensus toolkit for the
[TreeTools](https://ms609.github.io/TreeTools/) ecosystem, which also includes
['TreeDist'](https://ms609.github.io/TreeDist/) (tree distances and
information-theoretic consensus) and
['TreeSearch'](https://ms609.github.io/TreeSearch/) (phylogenetic search).
TreeTools itself remains the fast engine for the strict and majority-rule
consensus, which 'ConsTree' exposes through consistently-named wrappers.

The ['Rogue'](https://ms609.github.io/Rogue/) package identifies unstable
('rogue') leaves whose removal can improve the resolution and support of a
consensus tree; dropping rogues before summarising with 'ConsTree' often yields
a better-resolved result.

## Citation and attribution

The algorithms repackaged here originate in the FACT, FACT2 and FDCT prototypes
of Jesper Jansson and colleagues, incorporated with permission; see
`inst/REFERENCES.bib` for the source papers, and `DESCRIPTION` for copyright
attribution. 'ConsTree' is released under GPL (≥ 3).

Please note that this project is released with a
[Contributor Code of Conduct](https://ms609.github.io/TreeTools/CODE_OF_CONDUCT.html).
By contributing, you agree to abide by its terms.

