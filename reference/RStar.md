# R\* consensus tree

`RStar()` returns the R\* consensus (Degnan et al. 2009) of a set of
**rooted** trees.

## Usage

``` r
RStar(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

## Value

`RStar()` returns the consensus tree, an object of class `phylo`. It is
rooted by construction, but when the resolved triplets leave the deepest
divergence unresolved the root is a polytomy (in the limit, a star).
Such a tree has a root node of degree greater than two, which
[`ape::is.rooted()`](https://rdrr.io/pkg/ape/man/root.html) reports as
unrooted; re-root or resolve downstream if a strictly binary root is
required. All input trees must share the same set of (unique) tip
labels.

## Details

The R\* consensus is a rooted-triplet method. For every set of three
leaves it tallies, across the input trees, the three possible resolved
rooted triplets (`ab|c`, `ac|b`, `bc|a`) and keeps the one that is
*uniquely favoured*: the resolution appearing in strictly more input
trees than **each** of the other two, considered separately. This is a
strict **plurality** rule, not a majority rule: a triplet can be
uniquely favoured with far fewer than half the votes (e.g. a 2–1–1 split
among four trees), and any tie leaves those taxa unresolved. The kept
triplets form the set of *majority resolved triplets*, \\R\_{maj}\\. The
R\* tree is then the unique tree whose clades are exactly the **strong
clusters** of \\R\_{maj}\\ (Degnan et al. 2009; Jansson et al. 2016) : a
leaf set `A` is a clade if and only if, for *every* pair of leaves in
`A` and *every* leaf `x` outside `A`, the triplet grouping that pair
against `x` is uniquely favoured. Equivalently, R\* is the most resolved
tree that displays no resolved triplet outside \\R\_{maj}\\.

R\* is always a *refinement* of the majority-rule consensus (every
majority clade also appears in `RStar()`) and is a statistically
consistent estimator of a species tree from gene trees. Unlike
[`Local()`](https://constree.github.io/reference/Local.md) it is
**polynomial** and imposes no hard leaf-count limit: it stores no dense
\\n^3\\ structure (memory is \\O(kn^2)\\ for `k` input trees), so the
practical bound is running time (the tally is \\O(kn^3)\\), not a memory
wall – quite unlike
[`Local()`](https://constree.github.io/reference/Local.md)'s
exact-exponential 20-leaf bound.

Like [`Adams()`](https://constree.github.io/reference/Adams.md), R\* is
a rooted method: triplet states depend on the rooting, so input trees
are treated as rooted on their current root. Root the trees as you
intend before calling `RStar()`.

## Construction and conventions

- *Fans.* When a non-binary input tree leaves three leaves unresolved (a
  fan), that tree does not count toward any resolution of that triplet;
  fans have no impact on \\R\_{maj}\\ (Jansson et al. 2016) .

- *Ties.* If no resolution of a triple is uniquely favoured, that triple
  contributes nothing, leaving the affected taxa unresolved.

- *Existence and uniqueness.* The strong-cluster construction always
  yields a single, well-defined tree (Lemma 1.1 of (Jansson et al. 2016)
  ); there is no incompatibility or "build-failure" case to resolve.

- *Algorithm.* Following (Jansson et al. 2016) , the tree is built from
  a leaf-pair similarity – for each pair `a`, `b`, the number of leaves
  `x` for which `ab|x` lies in \\R\_{maj}\\ – whose single-linkage
  (Apresjan) clusters form a laminar superset of the strong clusters;
  those are then filtered exactly and assembled. The tally runs in
  \\O(kn^3)\\ time using per-tree constant-time LCA queries and stores
  no \\n^3\\ tensor (memory \\O(kn^2)\\); the assembly is about
  \\O(n^3)\\. The asymptotically sub-cubic bounds of (Jansson and Sung
  2013; Jansson et al. 2016) rely on fast-matrix-multiplication and
  dynamic-connectivity machinery that does more work than this for every
  practical leaf count, so they are not used.

## References

Degnan JH, DeGiorgio M, Bryant D, Rosenberg NA (2009). “Properties of
consensus methods for inferring species trees from gene trees.”
*Systematic Biology*, **58**(1), 35–54.
[doi:10.1093/sysbio/syp008](https://doi.org/10.1093/sysbio/syp008) .  
  
Jansson J, Sung W (2013). “Constructing the R\* consensus tree of two
trees in subcubic time.” *Algorithmica*, **66**(2), 329–345.
[doi:10.1007/s00453-012-9639-1](https://doi.org/10.1007/s00453-012-9639-1)
.  
  
Jansson J, Sung W, Vu H, Yiu S (2016). “Faster algorithms for computing
the R\* consensus tree.” *Algorithmica*, **76**(4), 1224–1244.
[doi:10.1007/s00453-016-0122-2](https://doi.org/10.1007/s00453-016-0122-2)
.

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Local()`](https://constree.github.io/reference/Local.md).

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
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
# Five trees whose majority signal recovers the species tree (((a,b),c),d):
trees <- c(
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, c), b), d);"),
  ape::read.tree(text = "(((b, c), a), d);")
)
RStar(trees) # (a, b) wins {a,b,c} by plurality (3 vs 1 vs 1)
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   a, b, c, d
#> 
#> Rooted; no branch length.
```
