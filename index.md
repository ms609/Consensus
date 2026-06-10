# ConsTree

‘ConsTree’ is an R package providing a comprehensive, efficient suite of
methods for summarizing a collection of phylogenetic trees — for example
a bootstrap or Bayesian posterior sample — as a single **consensus
tree**.

It builds on the tree and split infrastructure of
[‘TreeTools’](https://ms609.github.io/TreeTools/), and repackages a
family of asymptotically-efficient consensus algorithms developed by
Jesper Jansson and colleagues, several of which have not previously been
available in R.

## Consensus methods

### Split-selection methods

Each method takes a list of trees (or a `multiPhylo`) sharing the same
leaves and returns a `phylo`. They differ in *which* groupings (splits
or clusters) they retain:

| Function | Retains a grouping when… |
|----|----|
| [`Strict()`](https://constree.github.io/reference/Strict.md) | it occurs in **every** tree |
| [`Majority()`](https://constree.github.io/reference/Majority.md) / [`MajorityRule()`](https://constree.github.io/reference/Majority.md) | it occurs in **\> half** the trees (tunable via `p`) |
| [`Loose()`](https://constree.github.io/reference/Loose.md) | **no** tree contradicts it (semi-strict / combinable-component) |
| [`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md) | **more** trees display it than contradict it |
| [`Frequency()`](https://constree.github.io/reference/Frequency.md) | it is **more frequent than every** grouping that conflicts with it (frequency-difference) |
| [`Greedy()`](https://constree.github.io/reference/Greedy.md) | added greedily, most frequent first, if compatible with those already kept (extended majority-rule) |
| [`Adams()`](https://constree.github.io/reference/Adams.md) | constructed from the finest root-level partition shared by **every** tree (may introduce novel groupings; rooted) |
| [`Local()`](https://constree.github.io/reference/Local.md) | based on rooted triplets shared by **every** tree (minimum rooted/induced local consensus; ≤ 20 leaves) |
| [`RStar()`](https://constree.github.io/reference/RStar.md) | each rooted triplet grouping that wins a **plurality** against each alternative separately |

### Distance and branch-length summaries

A second family summarizes the trees through a distance or treespace
criterion rather than by selecting groupings:

| Function | Summary |
|----|----|
| [`Average()`](https://constree.github.io/reference/Average.md) | the tree best fitting the **mean path-length** (patristic) distances of the inputs |
| [`Quartet()`](https://constree.github.io/reference/Quartet.md) | an approximate median minimizing the total **quartet distance** to the inputs; often more resolved than majority-rule |
| `Transfer()` | a greedy consensus minimizing total **transfer distance** to the inputs; often more resolved than majority-rule |
| [`BHVMean()`](https://constree.github.io/reference/BHVMean.md) | the Fréchet **mean tree** in Billera–Holmes–Vogtmann treespace, with branch lengths; [`BHVDistance()`](https://constree.github.io/reference/BHVDistance.md), `BHVPairwiseDistances()` and [`BHVVariance()`](https://constree.github.io/reference/BHVMean.md) provide the supporting geodesic distances and dispersion |

## Usage

``` r

library("ConsTree")

trees <- ape::as.phylo(1:100, 8)   # 100 eight-leaf trees

Strict(trees)        # most conservative
Majority(trees)      # the familiar 50% majority-rule tree
Loose(trees)         # everything not actively contradicted
Frequency(trees)     # frequency-difference: often more resolved than majority
Greedy(trees)        # most resolved of the split-based summaries
Transfer(trees)      # minimizes transfer distance; often more resolved than majority-rule
```

## Installation

Install the development version from GitHub:

``` r

if (!require("remotes")) install.packages("remotes")
remotes::install_github("ms609/ConsTree")
```

‘ConsTree’ is not yet on CRAN.

## Relationship to other packages

‘ConsTree’ is the front-end consensus toolkit for the
[TreeTools](https://ms609.github.io/TreeTools/) ecosystem, which also
includes [‘TreeDist’](https://ms609.github.io/TreeDist/) (tree distances
and information-theoretic consensus) and
[‘TreeSearch’](https://ms609.github.io/TreeSearch/) (phylogenetic
search). TreeTools itself remains the fast engine for the strict and
majority-rule consensus, which ‘ConsTree’ exposes through
consistently-named wrappers.

[‘TreeDist’](https://ms609.github.io/TreeDist/) also offers a
complementary summary that ‘ConsTree’ does not duplicate: the tree from
a sample that has the lowest median clustering-information distance
(CID) to the others — a single *representative* of the sample rather
than a constructed consensus.

The quartet machinery underlying
[`Quartet()`](https://constree.github.io/reference/Quartet.md) builds on
the [‘Quartet’](https://ms609.github.io/Quartet/) package, which counts
the resolved- and shared-quartet statistics between trees; and the BHV
summaries relate to
[‘distory’](https://cran.r-project.org/package=distory), which computes
geodesic distances in the same treespace.

The [‘Rogue’](https://ms609.github.io/Rogue/) package identifies
unstable (‘rogue’) leaves whose removal can improve the resolution and
support of a consensus tree; dropping rogues before summarizing with
‘ConsTree’ often yields a better-resolved result.

## Citation and attribution

The algorithms repackaged here originate in the FACT, FACT2 and FDCT
prototypes of Jesper Jansson and colleagues, incorporated with
permission; see `inst/REFERENCES.bib` for the source papers, and
`DESCRIPTION` for copyright attribution. ‘ConsTree’ is released under
GPL (≥ 3).

Please note that this project is released with a [Contributor Code of
Conduct](https://ms609.github.io/TreeTools/CODE_OF_CONDUCT.html). By
contributing, you agree to abide by its terms.
