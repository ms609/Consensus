# Majority-rule consensus tree

`Majority()` returns the majority-rule consensus (Margush and McMorris
1981) : the tree containing each split that occurs in more than half of
the input trees. Raising `p` retains only splits present in a greater
proportion of trees, up to the strict consensus at `p = 1`.

## Usage

``` r
Majority(trees, p = 0.5)

MajorityRule(trees, p = 0.5)

MR(trees, p = 0.5)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

- p:

  Numeric between 0.5 and 1: the minimum proportion of trees that must
  contain a split for it to be retained. `p = 0.5` (the default) gives
  the majority-rule consensus; `p = 1` gives the strict consensus.

## Value

`Majority()` returns the consensus tree, an object of class `phylo`,
rooted as in the first entry of `trees`.

## Details

A thin wrapper around
[`TreeTools::Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.html).

## References

Margush T, McMorris FR (1981). “Consensus n-trees.” *Bulletin of
Mathematical Biology*, **43**(2), 239–244.
[doi:10.1007/BF02459446](https://doi.org/10.1007/BF02459446) .

## See also

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Majority(trees)
#> 
#> Phylogenetic tree with 8 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
