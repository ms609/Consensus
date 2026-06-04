# Local consensus tree

`Local()` returns the local consensus (Jansson et al. 2018) of a set of
**rooted** trees.

## Usage

``` r
Local(trees, type = c("rooted", "induced"))
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object; all entries must share the
  same leaf labels.

- type:

  Character specifying whether to compute the minimum rooted local
  consensus (`"rooted"`, the default; MinRLC) or the minimum induced
  local consensus (`"induced"`; MinILC).

## Value

`Local()` returns the consensus tree, an object of class `phylo`, rooted
by construction.

## Details

The local consensus is the most conservative tree consistent with the
rooted triplets shared by *every* input tree. Two variants are offered:
the minimum rooted local consensus (MinRLC, `type = "rooted"`) and the
minimum induced local consensus (MinILC, `type = "induced"`), which
differ in how the resolution of the result is scored. The tree is
assembled from the Aho-graph decomposition of the common-triplet set by
dynamic programming over the subsets of leaves.

The MinRLC and MinILC variants, and the exact exponential-time algorithm
used here to construct them, are due to (Jansson et al. 2018) ; this is
a direct port of the reference C++ from their `FDCT_new` toolkit, used
with permission.

**Complexity note:** the algorithm is exponential, so `Local()` is
limited to `n <= 20` leaves. Running time depends not only on `n` but on
how *congruent* the input trees are: when few rooted triplets are shared
(highly incongruent trees, particularly with `type = "induced"`), the
dynamic programming can become intractable even within the 20-leaf
limit. A long-running call can be interrupted (e.g. with Ctrl-C).

**No-valid-consensus:** when the entire leaf set forms a single
inseparable Aho-graph component (no common triplets separate any pair),
the algorithm has no valid consensus; `Local()` then returns a star tree
(all leaves attached directly to the root, fully unresolved). The
reference binary reports "No valid consensus found." for this case.

Input trees are treated as rooted on their current root; root the trees
as you intend before calling `Local()`.

## References

Jansson J, Rajaby R, Sung W (2018). “Minimal phylogenetic supertrees and
local consensus trees.” *AIMS Medical Science*, **5**(2), 181–203.
[doi:10.3934/medsci.2018.2.181](https://doi.org/10.3934/medsci.2018.2.181)
.

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Adams()`](https://constree.github.io/reference/Adams.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
# Two trees that agree on one cherry but disagree on overall topology
t1 <- ape::read.tree(text = "(1,((2,3),4));")
t2 <- ape::read.tree(text = "(1,((2,4),3));")
Local(list(t1, t2), "rooted")  # keeps clade {2,3,4} only
#> 
#> Phylogenetic tree with 4 tips and 2 internal nodes.
#> 
#> Tip labels:
#>   1, 2, 3, 4
#> 
#> Rooted; no branch length.
```
