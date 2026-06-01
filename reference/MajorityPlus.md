# Majority-rule (+) consensus tree

`MajorityPlus()` returns the majority-rule (+) consensus (Jansson et al.
2016) : a clade is retained when it occurs in more input trees than
contradict it – i.e. when the number of trees displaying the clade
exceeds the number of trees incompatible with it. A tree that is
compatible with a clade without displaying it counts neither for nor
against.

## Usage

``` r
MajorityPlus(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`MajorityPlus()` returns the consensus tree, an object of class `phylo`,
rooted as in the first entry of `trees`.

## Details

Every majority-rule split is retained (a split in more than half the
trees is contradicted by fewer than half), so `MajorityPlus()` contains
the majority-rule consensus
([`Majority()`](https://constree.github.io/reference/Majority.md)) and
may add further splits that are supported more often than they are
contradicted. The retained splits are necessarily mutually compatible,
so they define a valid tree.

## References

Jansson J, Shen C, Sung W (2016). “Improved algorithms for constructing
consensus trees.” *Journal of the ACM*, **63**(3), 1–24.
[doi:10.1145/2925985](https://doi.org/10.1145/2925985) .

## See also

Closely related:
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Loose()`](https://constree.github.io/reference/Loose.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
MajorityPlus(trees)
#> 
#> Phylogenetic tree with 8 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
