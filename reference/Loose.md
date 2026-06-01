# Loose consensus tree

`Loose()` returns the loose consensus, also known as the semi-strict or
combinable-component consensus (Bremer 1990) . It contains every split
that is contradicted by none of the input trees; equivalently, every
split that is compatible with each input tree.

## Usage

``` r
Loose(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`Loose()` returns the consensus tree, an object of class `phylo`, rooted
as in the first entry of `trees`.

## Details

The loose consensus is always at least as resolved as the strict
consensus
([`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md)) and
never includes a grouping that conflicts with any input tree. It is
incomparable with the majority-rule consensus
([`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md)):
a split present in most trees may still be contradicted by a minority,
and so excluded from the loose consensus, whereas a split occurring in a
single tree is retained if no other tree contradicts it.

## References

Bremer K (1990). “Combinable component consensus.” *Cladistics*,
**6**(4), 369–372.
[doi:10.1111/j.1096-0031.1990.tb00551.x](https://doi.org/10.1111/j.1096-0031.1990.tb00551.x)
.

## See also

Closely related:
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md).

Other consensus methods:
[`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md),
[`Average()`](https://ms609.github.io/ConsTree/reference/Average.md),
[`Frequency()`](https://ms609.github.io/ConsTree/reference/Frequency.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md),
[`Local()`](https://ms609.github.io/ConsTree/reference/Local.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`MajorityPlus()`](https://ms609.github.io/ConsTree/reference/MajorityPlus.md),
[`Quartet()`](https://ms609.github.io/ConsTree/reference/Quartet.md),
[`RStar()`](https://ms609.github.io/ConsTree/reference/RStar.md),
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Loose(trees)
#> 
#> Phylogenetic tree with 8 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
