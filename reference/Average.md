# Average consensus tree

`Average()` returns the *average consensus* (Lapointe and Cucumel 1997)
: the tree whose path-length (patristic) distances most closely match
the average of the path-length distances of the input trees. Informally,
it places each leaf at its mean position across the input trees, making
it a natural distance-based summary of a posterior sample –
complementing the split-based
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md) and
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md)
methods.

## Usage

``` r
Average(
  trees,
  method = c("fastme.bal", "ls", "nj", "bionj", "fastme.ols"),
  weights = NULL,
  scale = c("none", "max"),
  edgeLengths = NA,
  outgroup = NULL,
  check.labels = TRUE,
  lsControl = list()
)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

- method:

  Character specifying how to build the tree from the average distance
  matrix:

  - `"fastme.bal"` (the default) returns the balanced minimum-evolution
    tree (Desper and Gascuel 2002) : a fast, accurate approximation of
    the least-squares tree;

  - `"ls"` searches for the least-squares tree itself – the criterion
    under which Lapointe & Cucumel's averaging guarantee holds – using
    [`TreeSearch::LeastSquaresTree()`](https://ms609.github.io/TreeSearch/reference/LeastSquaresTree.html),
    a compiled non-negative least-squares NNI/SPR search;

  - `"nj"`, `"bionj"` and `"fastme.ols"` return the corresponding
    distance tree (Saitou and Nei 1987; Gascuel 1997) .

- weights:

  Numeric vector specifying the weight of each tree in the average (e.g.
  posterior probabilities), with one entry per tree. Defaults to `NULL`,
  which weights every tree equally – appropriate for a posterior sample,
  in which a tree's frequency already encodes its probability.

- scale:

  Character specifying whether to rescale each tree's distance matrix
  before averaging. `"none"` (the default) leaves matrices unscaled,
  appropriate when the trees are already commensurable (e.g. a single
  posterior sample). `"max"` divides each matrix by its largest entry,
  the standardization recommended by Lapointe & Cucumel when combining
  trees from heterogeneous sources whose absolute distances are not
  comparable.

- edgeLengths:

  Logical specifying whether to use branch lengths when computing
  path-length distances. The default, `NA`, uses branch lengths when
  *every* tree has them and otherwise counts edges; `TRUE` requires
  branch lengths; `FALSE` always counts edges (a topology-only summary).

- outgroup:

  Character vector specifying tip label(s) on which to root the result.
  Defaults to `NULL`, which returns an unrooted tree: path-length
  distances are unaffected by rooting, so the method is intrinsically
  unrooted.

- check.labels:

  Logical specifying whether to confirm that every tree describes the
  same leaves. The default, `TRUE`, is safer; `FALSE` is faster when the
  trees are known to share an identical leaf set.

- lsControl:

  Named list of further arguments for the least-squares search
  (`method = "ls"`), passed to
  [`TreeSearch::LeastSquaresTree()`](https://ms609.github.io/TreeSearch/reference/LeastSquaresTree.html);
  for example `list(spr = FALSE, maxHits = 5L, weight = "fm")` to use
  Fitch-Margoliash weighting. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html); ignored by other
  methods.

## Value

`Average()` returns the average consensus tree, an object of class
`phylo` with fitted branch lengths, unrooted unless `outgroup` is given.

## Details

The procedure has two steps (Lapointe and Cucumel 1997) :

1.  Compute the path-length distance matrix of each input tree (using
    branch lengths where present, otherwise counting edges), optionally
    rescaling each matrix, and average the matrices.

2.  Find the tree whose own path-length distances best fit this average
    matrix, in the least-squares sense.

Because the average of several path-length matrices is usually not
itself realisable by any tree (it violates the four-point condition),
step 2 is a fit, not an inversion. By default `Average()` approximates
it with the fast balanced minimum-evolution tree; `method = "ls"`
instead performs the exact least-squares search, which – being NP-hard
(Day 1987) – uses tree rearrangements, as did the original FITCH
implementation. Branch lengths are fitted by non-negative least squares,
so that the fitted distances are realisable by a tree, as the criterion
requires.

A lone input tree is its own average: it is returned (unrooted unless
`outgroup` is given) without refitting, and `method`, `scale`, `weights`
and `edgeLengths` then have no effect.

## References

Day WHE (1987). “Computational complexity of inferring phylogenies from
dissimilarity matrices.” *Bulletin of Mathematical Biology*, **49**(4),
461–467. [doi:10.1007/BF02458863](https://doi.org/10.1007/BF02458863)
.  
  
Desper R, Gascuel O (2002). “Fast and accurate phylogeny reconstruction
algorithms based on the minimum-evolution principle.” *Journal of
Computational Biology*, **9**(5), 687–705.
[doi:10.1089/106652702761034136](https://doi.org/10.1089/106652702761034136)
.  
  
Gascuel O (1997). “BIONJ: an improved version of the NJ algorithm based
on a simple model of sequence data.” *Molecular Biology and Evolution*,
**14**(7), 685–695.
[doi:10.1093/oxfordjournals.molbev.a025808](https://doi.org/10.1093/oxfordjournals.molbev.a025808)
.  
  
Lapointe F, Cucumel G (1997). “The average consensus procedure:
combination of weighted trees containing identical or overlapping sets
of taxa.” *Systematic Biology*, **46**(2), 306–312.
[doi:10.1093/sysbio/46.2.306](https://doi.org/10.1093/sysbio/46.2.306)
.  
  
Saitou N, Nei M (1987). “The neighbor-joining method: a new method for
reconstructing phylogenetic trees.” *Molecular Biology and Evolution*,
**4**(4), 406–425.
[doi:10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
.

## See also

Split-based summaries:
[`Strict()`](https://ms609.github.io/ConsTree/reference/Strict.md),
[`Majority()`](https://ms609.github.io/ConsTree/reference/Majority.md).

Other consensus methods:
[`Adams()`](https://ms609.github.io/ConsTree/reference/Adams.md),
[`Frequency()`](https://ms609.github.io/ConsTree/reference/Frequency.md),
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
trees <- ape::rmtree(5, 8)         # five random eight-leaf trees
Average(trees)                     # fast (balanced minimum evolution) default
#> 
#> Phylogenetic tree with 8 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t5, t6, t2, t4, t7, t3, ...
#> 
#> Unrooted; includes branch length(s).
# \donttest{
if (requireNamespace("TreeSearch", quietly = TRUE) &&
    exists("LeastSquaresTree", where = asNamespace("TreeSearch"),
           mode = "function")) {
  Average(trees, method = "ls")    # faithful least-squares fit (slower)
}
#> 
#> Phylogenetic tree with 8 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t5, t6, t4, t1, t8, t7, ...
#> 
#> Unrooted; includes branch length(s).
# }
```
