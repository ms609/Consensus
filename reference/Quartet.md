# Consensus tree minimizing quartet distance

Construct a consensus tree that minimizes the sum of symmetric quartet
distances to a set of input trees, using a greedy add-and-prune
heuristic.

## Usage

``` r
Quartet(
  trees,
  init = c("majority", "empty", "extended"),
  greedy = c("best", "first")
)
```

## Arguments

- trees:

  Object of class `multiPhylo` specifying the input trees. All trees
  must share the same tip labels. Trees may be non-binary (polytomies
  are handled correctly).

- init:

  Character string specifying the initial tree:

  - `"majority"` (default): start from the majority-rule consensus.

  - `"empty"`: start from a star tree (purely additive).

  - `"extended"`: start from the extended (greedy) majority-rule
    consensus.

- greedy:

  Character string specifying the greedy strategy:

  - `"best"` (default): evaluate all candidates and pick the single
    highest-benefit action at each step.

  - `"first"`: pick the first improving action encountered (faster but
    may give a slightly worse result).

## Value

`Quartet()` returns a consensus tree, an object of class `phylo`,
unrooted.

## Details

The majority-rule consensus minimizes the sum of Robinson-Foulds
distances to the input trees. Analogously, `Quartet()` finds an
approximate median tree under the symmetric quartet distance (Takazawa
et al. 2026) , which counts both false-positive and false-negative
resolved quartets equally.

Because the quartet distance gives greater weight to deep branches
(which resolve more quartets), quartet consensus trees tend to be more
resolved than majority-rule trees, especially when phylogenetic signal
is low.

The algorithm pools all splits observed across input trees and maintains
a quartet profile: for each of the \\\binom{n}{4}\\ quartets, a count of
how many input trees resolve it as each of the three possible
topologies. Splits are greedily added to (or removed from) the consensus
when doing so reduces the total symmetric quartet distance to the input
trees. Candidate splits must be compatible with all currently included
splits.

The function supports trees with up to 100 tips. For larger trees, the
explicit quartet enumeration becomes prohibitively expensive.

## References

Takazawa Y, Takeda A, Hayamizu M, Gascuel O (2026). “Outperforming the
majority-rule consensus tree using fine-grained dissimilarity measures.”
*bioRxiv*.
[doi:10.64898/2026.03.16.712085](https://doi.org/10.64898/2026.03.16.712085)
.

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md)

## Examples

``` r
library(TreeTools)
#> Loading required package: ape

# Generate bootstrap-like trees
trees <- as.phylo(1:30, nTip = 8)

# Quartet consensus
qc <- Quartet(trees)
plot(qc)


# Compare resolution with majority-rule
mr <- UnrootTree(Consensus(trees, p = 0.5))
cat("Majority-rule splits:", NSplits(mr), "\n")
#> Majority-rule splits: 2 
cat("Quartet consensus splits:", NSplits(qc), "\n")
#> Quartet consensus splits: 3 
```
