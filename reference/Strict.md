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
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md).

Other consensus methods:
[`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md),
[`Average()`](https://ms609.github.io/ConsTree/reference/Average.md),
[`Frequency()`](https://ms609.github.io/ConsTree/reference/Frequency.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md),
[`Local()`](https://ms609.github.io/ConsTree/reference/Local.md),
[`Loose()`](https://ms609.github.io/ConsTree/reference/Loose.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`MajorityPlus()`](https://ms609.github.io/ConsTree/reference/MajorityPlus.md),
[`Quartet()`](https://ms609.github.io/ConsTree/reference/Quartet.md),
[`RStar()`](https://ms609.github.io/ConsTree/reference/RStar.md)

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
