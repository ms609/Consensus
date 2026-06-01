# Fréchet mean and variance in BHV treespace

`BHVMean()` returns the Fréchet (Karcher) mean of a set of trees in BHV
treespace: the tree that minimizes the sum of squared geodesic distances
to the sample (Brown and Owen 2020) . As there is no known closed form,
it is approximated by the iterative law-of-large-numbers algorithm of
Sturm (2003) and Miller et al. (2015) : starting from a sample tree,
each step walks a fraction \\1/(k+1)\\ of the way along the geodesic
towards a randomly chosen sample tree.

## Usage

``` r
BHVMean(trees, tolerance = 1e-04, maxIter = 100000L, cauchyLength = 10L)

BHVVariance(trees, mean = NULL, type = c("average", "sum"))
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels and carry `edge.length`.

- tolerance:

  Numeric convergence threshold, *relative* to the sample standard
  deviation: iteration stops once `cauchyLength` consecutive steps each
  move the estimate less than `tolerance` times the sample standard
  deviation. Smaller values give a more precise mean at the cost of more
  iterations.

- maxIter:

  Integer specifying the maximum number of iterations.

- cauchyLength:

  Integer specifying the number of consecutive small steps required to
  declare convergence.

- mean:

  Object of class `phylo` specifying a pre-computed mean tree; computed
  via `BHVMean()` if `NULL` (the default).

- type:

  Character specifying whether to return the mean squared distance
  (`"average"`, the default) or the total squared distance (`"sum"`).

## Value

`BHVMean()` returns the mean tree, an object of class `phylo`, with
attributes `iterations` (number of steps taken) and `converged`. Because
the step length shrinks as \\1/(k+1)\\, `converged = TRUE` indicates
that successive estimates have stopped moving appreciably (the stopping
rule was met before `maxIter`), not a guaranteed bound on the distance
to the exact mean; tighten `tolerance` for greater precision.

`BHVVariance()` returns a single non-negative number.

## Details

`BHVVariance()` returns the Fréchet variance: by default the mean
squared geodesic distance from the sample to its mean, \\(1/r)\sum_i
d(\bar T, T_i)^2\\; with `type = "sum"`, the total \\\sum_i d(\bar T,
T_i)^2\\.

The mean is "sticky": perturbing one sample tree need not move it, and
it is pulled towards lower-dimensional (less resolved) orthants, so it
may be unresolved even when the sample trees are binary (Brown and Owen
2020) .

## References

Brown DG, Owen M (2020). “Mean and variance of phylogenetic trees.”
*Systematic Biology*, **69**(1), 139–154.
[doi:10.1093/sysbio/syz041](https://doi.org/10.1093/sysbio/syz041) .  
  
Miller E, Owen M, Provan JS (2015). “Polyhedral computational geometry
for averaging metric phylogenetic trees.” *Advances in Applied
Mathematics*, **68**, 51–91.
[doi:10.1016/j.aam.2015.04.002](https://doi.org/10.1016/j.aam.2015.04.002)
.  
  
Sturm K (2003). “Probability measures on metric spaces of nonpositive
curvature.” In *Heat Kernels and Analysis on Manifolds, Graphs, and
Metric Spaces*, volume 338 of *Contemporary Mathematics*, 357–390.
[doi:10.1090/conm/338/06080](https://doi.org/10.1090/conm/338/06080) .

## See also

Other BHV summaries:
[`BHVDistance()`](https://constree.github.io/reference/BHVDistance.md)

## Examples

``` r
set.seed(0)
trees <- lapply(1:25, function(i) {
  tree <- TreeTools::RandomTree(6, root = FALSE)
  tree$edge.length <- runif(nrow(tree$edge))
  tree
})
meanTree <- BHVMean(trees)
BHVVariance(trees, mean = meanTree)
#> [1] 1.254778
```
