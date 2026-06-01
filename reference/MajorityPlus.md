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
([`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md))
and may add further splits that are supported more often than they are
contradicted. The retained splits are necessarily mutually compatible,
so they define a valid tree.

## References

Jansson J, Shen C, Sung W (2016). “Improved algorithms for constructing
consensus trees.” *Journal of the ACM*, **63**(3), 1–24.
[doi:10.1145/2925985](https://doi.org/10.1145/2925985) .

## See also

Closely related:
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md),
[`Loose()`](https://ms609.github.io/ConsTree/reference/Loose.md).

Other consensus methods:
[`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md),
[`Average()`](https://ms609.github.io/ConsTree/reference/Average.md),
[`Frequency()`](https://ms609.github.io/ConsTree/reference/Frequency.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md),
[`Local()`](https://ms609.github.io/ConsTree/reference/Local.md),
[`Loose()`](https://ms609.github.io/ConsTree/reference/Loose.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`Quartet()`](https://ms609.github.io/ConsTree/reference/Quartet.md),
[`RStar()`](https://ms609.github.io/ConsTree/reference/RStar.md),
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md)

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
