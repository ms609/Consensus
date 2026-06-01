# Greedy (extended majority-rule) consensus tree

`Greedy()` returns the greedy consensus, also termed the extended
majority-rule consensus (Bryant 2003) . Distinct splits are considered
in decreasing order of their frequency across the input trees; each is
added to the growing consensus if it is compatible with every split
already accepted. The result is typically more resolved than the
majority-rule consensus
([`Majority()`](https://constree.github.io/reference/Majority.md)), and
contains it.

## Usage

``` r
Greedy(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`Greedy()` returns the consensus tree, an object of class `phylo`,
rooted as in the first entry of `trees`.

## Details

Splits that occur equally often are considered in a fixed, reproducible
order, so the result is deterministic. Where several mutually
incompatible splits are equally frequent, a different (but equally
valid) greedy resolution may be returned by other software.

## References

Bryant D (2003). “A classification of consensus methods for
phylogenetics.” In Janowitz MF, Lapointe F, McMorris FR, Mirkin B,
Roberts FS (eds.), *Bioconsensus*, volume 61 of *DIMACS Series in
Discrete Mathematics and Theoretical Computer Science*, 163–184.
American Mathematical Society.
[doi:10.1090/dimacs/061/11](https://doi.org/10.1090/dimacs/061/11) .

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Loose()`](https://constree.github.io/reference/Loose.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Greedy(trees)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
