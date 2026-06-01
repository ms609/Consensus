# Adams consensus tree

`Adams()` returns the Adams consensus (Adams 1972) , a summary of a set
of **rooted** trees.

## Usage

``` r
Adams(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`Adams()` returns the consensus tree, an object of class `phylo`, rooted
by construction.

## Details

The Adams consensus is defined recursively. At the top level, the leaves
are partitioned into the coarsest grouping that every input tree
refines: two leaves fall in the same block if and only if, in every
input tree, they sit below the same child of the most recent common
ancestor of the full leaf set. Each block becomes a child of the
consensus root, and the construction repeats within each block.

Unlike the split-based methods
([`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md)), the
Adams consensus can contain a cluster that appears in no individual
input tree, and it depends on how the input trees are rooted: it is a
statement about nesting, not about unrooted bipartitions. Input trees
are treated as rooted on their current root; root the trees as you
intend before calling `Adams()`.

## References

Adams EN (1972). “Consensus techniques and the comparison of taxonomic
trees.” *Systematic Zoology*, **21**(4), 390–397.
[doi:10.2307/2412432](https://doi.org/10.2307/2412432) .

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md).

Other consensus methods:
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
# Two rooted trees that disagree only on the position of one leaf
trees <- c(ape::read.tree(text = "(((a, b), c), d);"),
           ape::read.tree(text = "(((a, b), d), c);"))
Adams(trees) # keeps the clade (a, b); leaves c, d unresolved at the root
#> 
#> Phylogenetic tree with 4 tips and 2 internal nodes.
#> 
#> Tip labels:
#>   a, b, c, d
#> 
#> Unrooted; no branch length.
```
