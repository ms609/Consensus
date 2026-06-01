# Frequency-difference consensus tree

`Frequency()` returns the frequency-difference consensus: a split is
retained when it occurs strictly more often than every split that
conflicts with it. Equivalently, among each set of mutually incompatible
splits, the consensus keeps a split only if it is strictly more frequent
than all its rivals.

## Usage

``` r
Frequency(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`Frequency()` returns the consensus tree, an object of class `phylo`,
rooted as in the first entry of `trees`.

## Details

The frequency-difference consensus is at least as resolved as the
majority-rule consensus
([`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md))
– a split in more than half the trees is necessarily more frequent than
any conflicting split – and is contained within the greedy consensus
([`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md)).
The retained splits are mutually compatible, so they define a valid
tree.

An efficient algorithm for constructing this consensus was given by
(Jansson et al. 2024) ; the present implementation computes the same
tree directly from pooled split frequencies.

## References

Jansson J, Sung W, Tabatabaee SA, Yang Y (2024). “A Faster Algorithm for
Constructing the Frequency Difference Consensus Tree.” In Beyersdorff O,
Kanté MM, Kupferman O, Lokshtanov D (eds.), *41st International
Symposium on Theoretical Aspects of Computer Science (STACS 2024)*,
volume 289 of *Leibniz International Proceedings in Informatics
(LIPIcs)*, 43:1–43:17.
[doi:10.4230/LIPIcs.STACS.2024.43](https://doi.org/10.4230/LIPIcs.STACS.2024.43)
.

## See also

Closely related:
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`MajorityPlus()`](https://ms609.github.io/ConsTree/reference/MajorityPlus.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md).

Other consensus methods:
[`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md),
[`Average()`](https://ms609.github.io/ConsTree/reference/Average.md),
[`Greedy()`](https://ms609.github.io/ConsTree/reference/Greedy.md),
[`Local()`](https://ms609.github.io/ConsTree/reference/Local.md),
[`Loose()`](https://ms609.github.io/ConsTree/reference/Loose.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md),
[`MajorityPlus()`](https://ms609.github.io/ConsTree/reference/MajorityPlus.md),
[`Quartet()`](https://ms609.github.io/ConsTree/reference/Quartet.md),
[`RStar()`](https://ms609.github.io/ConsTree/reference/RStar.md),
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Frequency(trees)
#> 
#> Phylogenetic tree with 8 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
