# Geodesic (BHV) distance between trees

`BHVDistance()` returns the geodesic distance between phylogenetic trees
with edge lengths in the Billera-Holmes-Vogtmann (BHV) treespace
(Billera et al. 2001) , computed with the polynomial-time GTP algorithm
of Owen and Provan (2011) .

## Usage

``` r
BHVDistance(tree1, tree2 = NULL)

BHVDist(tree1, tree2 = NULL)

BHV(tree1, tree2 = NULL)
```

## Arguments

- tree1, tree2:

  A [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) tree, a
  `multiPhylo` object, or a list of `phylo` trees, each carrying
  `edge.length`. All trees across both arguments must share the same
  leaf labels. `tree2` may be omitted when `tree1` is a collection, in
  which case all pairwise distances within `tree1` are returned.

## Value

- Both `tree1` and `tree2` are single trees: a single non-negative
  number.

- One argument is a single tree, the other a collection: a named numeric
  vector with one entry per tree in the collection.

- `tree2` is omitted, or `tree1` and `tree2` are the same collection: a
  [`stats::dist`](https://rdrr.io/r/stats/dist.html) object of all
  pairwise distances.

- `tree1` and `tree2` are different collections: a numeric matrix with
  rows corresponding to `tree1` and columns to `tree2`.

## References

Billera LJ, Holmes SP, Vogtmann K (2001). “Geometry of the space of
phylogenetic trees.” *Advances in Applied Mathematics*, **27**(4),
733–767.
[doi:10.1006/aama.2001.0759](https://doi.org/10.1006/aama.2001.0759) .  
  
Owen M, Provan JS (2011). “A fast algorithm for computing geodesic
distances in tree space.” *IEEE/ACM Transactions on Computational
Biology and Bioinformatics*, **8**(1), 2–13.
[doi:10.1109/TCBB.2010.3](https://doi.org/10.1109/TCBB.2010.3) .

## See also

Other BHV summaries:
[`BHVMean()`](https://constree.github.io/reference/BHVMean.md)

## Examples

``` r
set.seed(2)
AddEdgeLengths <- function(tree) { tree$edge.length <- runif(nrow(tree$edge)); tree }
trees <- lapply(1:4, function(i) AddEdgeLengths(TreeTools::RandomTree(8, root = TRUE)))
t1 <- trees[[1]]; t2 <- trees[[2]]

BHVDistance(t1, t2)               # scalar
#> [1] 2.742226
BHVDistance(t1, trees)            # named vector
#> [1] 0.000000 2.742226 3.088498 3.158640
BHVDistance(trees)                # dist (pairwise)
#>          1        2        3
#> 2 2.742226                  
#> 3 3.088498 3.094400         
#> 4 3.158640 3.187792 3.344174
BHVDistance(trees[1:2], trees[3:4])  # matrix
#>          [,1]     [,2]
#> [1,] 3.088498 3.158640
#> [2,] 3.094400 3.187792
```
