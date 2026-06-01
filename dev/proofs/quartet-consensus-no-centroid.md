# Quartet never resolves a "centroid" split

**Claim (no-centroid theorem).** Run with the default `neverDrop = TRUE`, every
non-trivial split in the tree returned by `Quartet(trees)` is present
in **at least one** input tree. Consequently `Quartet` can never assert
a relationship (bipartition) that *no* input tree contains — there is no
"centroid" / averaged placement that all input trees contradict.

This answers the question "can QC place a taxon where the data never puts it?":
**no, not as a single relationship — by construction.**

---

## Proof (constructive)

The candidate set of splits is fixed up front to the set of *observed* splits,
and the greedy search only ever toggles membership of candidates in that set; it
never synthesises a new bipartition.

1. **The candidate pool is the deduplicated union of input-tree splits.**
   `R/Quartet.R` converts each input tree to its split matrix
   (`splitsList <- lapply(trees, as.Splits)`, `Quartet.R:157-160`) and
   passes the list to `cpp_quartet_consensus`. In C++,
   `pool_splits()` (`src/Quartet.cpp:181-284`) iterates over every
   tree `t` and every split `s` of that tree, canonicalises it, and inserts it
   into the pool, incrementing `pool.count[idx]`. A pool entry is created
   *only* from some tree's split (`canon_buf` is copied from `mat(s, ·)`,
   `:228-232`). Hence every pooled split has `count >= 1`: it occurs in at least
   one input tree. The pool contains nothing else.

2. **The greedy only adds/removes pooled splits.** Both `greedy_best`
   (`:834-884`) and `greedy_first` (`:891-930`) loop over pooled indices
   `idx in 0..M-1` (`M = pool.n_splits`) and take exactly one of three actions:
   `do_add(idx)`, `do_remove(idx)` (both with `idx` a pool index), or
   `do_drop(tip)`. There is **no** code path that constructs a bipartition
   outside the pool. The included set is therefore always a subset of the pool.

3. **The returned splits are the included pooled splits.** The main routine
   returns `res$splits` = the matrix of currently-included pooled splits, which
   `R/Quartet.R:172-183` turns into the output tree via `as.phylo`.

Combining (1)–(3): output splits ⊆ pooled splits ⊆ {splits occurring in ≥1 input
tree}. ∎

**Drop caveat.** With `neverDrop = FALSE`, `do_drop` removes a tip and remaining
splits are remapped to the surviving tips (`:806-813`). A remapped split is the
restriction of an observed split to the active tip set, i.e. it is present in
the correspondingly pruned input tree — so the guarantee holds on the *reduced*
leaf set. The clean, full-leaf-set statement is for the default `neverDrop = TRUE`.

---

## Empirical confirmation

`dev/proofs/centroid-search.R` searches for a counterexample across two regimes
(seed 1; build of branch `consensus-rogues`, 2026-05-31):

| regime | what varies | trials | min split-freq | freq-0 (centroid) configs |
|--------|-------------|-------:|---------------:|--------------------------:|
| single rogue | one rogue's position on a fixed backbone | 1500 | **0.5000** | **0** |
| conflict | random-topology mixtures, large unstable clades, fully random trees | 1200 | **0.0500** (= 1/20) | **0** |

No centroid was ever produced. The minimum is `1/k` (a split present in exactly
one of `k` trees) in the conflict regime — confirming both that the floor is
strictly positive (the theorem) and that **QC does resolve minority splits**
(its "more resolved than majority-rule" behaviour: it resolves a split when its
*quartets* are net-supported, which can occur at low *bipartition* frequency).

### Methodological note (a real false-positive, caught)

An earlier search reported 74/1500 "centroids" (min freq 0.0000). These were
**false positives** from comparing only the "smaller side" of each split:
for a balanced `n/2 | n/2` split the smaller-side tiebreak is inconsistent
between trees, so a present split can look absent. Fix: canonicalise each split
to the side containing the lexicographically-first tip before comparing
(`canonKeys` in the search script). With the fix, the count drops to 0,
matching the theorem.

---

## Scope: what is and isn't proven

- **Proven:** QC never asserts a *single* relationship (split) absent from all
  inputs. The PI's literal concern — "park it somewhere ALL trees agree it does
  NOT [go]" — cannot happen for any one relationship.
- **The central/Adams placement is a *combination* of individually-observed
  splits**, which the per-split theorem does not forbid (keeping both deep
  half-clades `{1234}|rest` and `{5678}|rest` while a rogue sits between them —
  each present ~50% of the time, co-occurring in no tree). So this case must be
  checked empirically, including the realistic hard case of a rogue split
  ~50/50 between two distant homes.
- **Tested — QC does not centre, even at an exact tie** (`dev/`-adjacent probe,
  balanced-8, R = sister-`t1` vs sister-`t8`, opposite halves):
  | t1:t8 | QC places R at | R excluded from *both* halves (central)? |
  |------:|----------------|:--:|
  | 39:41 … 30:50 | `{R,t8}` (the larger mode) | no |
  | **40:40** | `{R,t1,t2}` (a 50%-supported clade) | **no** |
  At an exact tie QC still commits R to one 50%-mode clade rather than centring
  it — because *no* input tree places R centrally, so R's quartets never favour
  the central resolution. This matches the PI's stated preference ("somewhere it
  might go, >½ the time" over "where all trees agree it does not").

## What QC *does* do that is worth flagging: marginal minority splits

QC resolves a split when its **quartets** are net-supported (penalty = 1 ⇒ >½ of
resolving trees agree), which is **not** the same as the split being a >½
*bipartition*. In high-conflict data QC resolves splits with bipartition
frequency as low as **15–25%** on the strength of a **razor-thin ~51% quartet
majority** (measured: lowest-frequency QC splits across random-conflict configs
had quartet support 0.510–0.513). This is QC's "more resolved than majority-rule"
behaviour — the precision/accuracy trade-off of Smith (2019), governed by
`penalty`. It is *not* a centroid (the split is present in some trees and its
quartets are net-supported), but it *is* where a reader could be misled into
reading a barely-supported relationship as firm.

## Implication for the consensus contribution

Because QC's topology is already centroid-free and never centres a taxon, a
"dispersion correction" that *changes the topology* adds nothing over the
existing `Quartet`: the tree it would produce is the tree QC already
produces. The honest-uncertainty contribution therefore lives **not** in moving
taxa but in **calibrating / reporting confidence on the (unchanged) QC tree** —
controlling or flagging the marginal minority splits above (via the `penalty`
threshold, or a per-relationship confidence annotation), so a quartet-supported-
but-bipartition-rare split is not read as firm.
