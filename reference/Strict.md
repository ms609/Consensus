# Strict consensus tree

`Strict()` returns the strict consensus of a set of trees: the tree that
contains exactly those splits (clades) present in *every* input tree
(Day 1985) .

## Usage

``` r
Strict(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`Strict()` returns the consensus tree, an object of class `phylo`,
rooted as in the first entry of `trees`.

## Details

This is a thin wrapper around
[`TreeTools::Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.html)
with `p = 1`, provided so that every consensus method in this package is
reachable through a common, consistently named interface.

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .

## See also

Less conservative summaries:
[`Majority()`](https://constree.github.io/reference/Majority.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Strict(trees)
#> 
#> Phylogenetic tree with 8 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
