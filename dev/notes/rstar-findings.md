# R\* — validation & correction record

How `RStar()` was validated, and the mid-course correction from a BUILD-based
assembly to the (correct) strong-cluster construction. Companion to
`rstar-spec.md` (literature) and `rstar-design.md` (design).

## TL;DR

`RStar()` is now an exact implementation of the R\* consensus tree as defined by
**Jansson, Sung, Vu & Yiu (2016)** (`Jansson2016a`), confirmed against the
definition three independent ways. An earlier BUILD-based assembly was found to
violate the definition on conflicting input and was replaced; the open questions
flagged in earlier drafts (OQ1/OQ2/OQ3) are **resolved by the primary source** and
are no longer open.

## The definition (Jansson et al. 2016, §1.1 + Lemma 1.1)

- `#ab|c` = number of input trees in which `ab|c` is a *consistent* resolved
  triplet. **Fan triplets have no impact** — they count toward nothing. (This
  resolves OQ1: fans abstain. Confirmed.)
- `R_maj = { ab|c : #ab|c > max(#ac|b, #bc|a) }` — strict **plurality**, no `k/2`
  threshold.
- The R\* tree is the unique tree `τ` with `r(τ) ⊆ R_maj` maximising internal
  nodes; **Lemma 1.1**: its clusters are exactly the **strong clusters** of
  `R_maj` (A is strong iff `aa'|x ∈ R_maj` for all pairs `a,a' ∈ A` and all
  `x ∉ A`). The tree always exists and is unique — there is **no
  build-failure/collapse case** (resolves OQ2) and **no over-resolution latitude**
  (resolves OQ3).

## The correction (what changed and why)

The first implementation assembled the tree with Aho **BUILD** on `R_maj`. BUILD
returns the least-resolved tree *consistent with* (displaying) `R_maj`, which is
**not** the R\* tree: it can display resolved triplets *outside* `R_maj`,
violating `r(τ) ⊆ R_maj`. Quantified over 400 random conflicting trials
(n = 5–12, k = 3–9), the BUILD version diverged from the definition both ways:

| Divergence (BUILD version) | Trials affected |
|---|---|
| Dropped a rooted **majority** clade (under-resolution) | 46 / 400 |
| Produced a clade violating the for-all-outgroup rule (over-resolution) | 51 / 400 |

Both are now understood as symptoms of using BUILD (an *exists*-outgroup grouping)
instead of strong clusters (the *for-all*-outgroup definition). `src/rstar.cpp`
was rewritten: keep the (correct) plurality tally, then assemble by extracting
strong clusters — candidate clusters are the single-linkage threshold components
of `s(a,b) = #{w : ab|w ∈ R_maj}` (a provable superset: strong ⊆ Apresjan ⊆
threshold-components), each tested against the strong-cluster rule, the laminar
survivors nested into the tree.

## Validation (current, strong-cluster version)

| Check | Result |
|---|---|
| Brute-force strong-cluster oracle: `RStar` clades == strong clusters of `R_maj`, over **all** 2ⁿ subsets, independent `ape::mrca` tally (n ≤ 12) | **60 / 60 trials exact** |
| Majority-rule refinement (every rooted majority clade ∈ R\*) | **0 / 400 violations** (was 46/400) |
| Over-resolution (clade violating for-all-outgroup rule) | **0 / 400** (was 51/400) |
| Strict-clade refinement (unanimous clades survive) | **0 / 150** |
| Identity: `RStar(k×T) == T`, n ≤ 40 | **0 / 60 fail** |
| Congruent oracle: `RStar` == FDCT `aho-build` (n = 4–16) | **6 / 6 MATCH** |
| Worked example / plurality-not-majority / tie → polytomy | exact |

The majority-refinement result is **provable**, not merely observed: a majority
clade C (in > k/2 trees) forces `#aa'|x > k/2` for every internal pair and outside
x, hence `aa'|x ∈ R_maj`, hence C is a strong cluster. The congruent-oracle match
is also provable (congruent input ⇒ unanimous triplets ⇒ strong clusters ≡ BUILD).

Harness: `dev/oracle/rstar/check-rstar.R` (must-pass: identity + congruent
oracle + strict + majority refinement), `check-strong-clusters.R` (brute-force
gate), `explore-rstar.R` / `diagnose-clade-rule.R` (broad sweeps, now 0/0 in both
directions — retained as before/after evidence of the BUILD→strong-cluster fix).
Package tests: `tests/testthat/test-rstar.R` (includes a self-contained
brute-force strong-cluster oracle + majority-refinement).

## Status of the former "open questions"

- **OQ1 (fans):** resolved — fans have no impact (Jansson et al. 2016, §1.1).
  Implementation already did this.
- **OQ2 (inconsistency/collapse):** resolved — does not arise; the strong-cluster
  tree always exists and is unique (Lemma 1.1).
- **OQ3 (assembly / over-resolution):** resolved — the definition mandates strong
  clusters (`r(τ) ⊆ R_maj`); BUILD was simply the wrong construction, now fixed.

## Deferred (performance only, not correctness)

The implementation is correctness-first: O(kn³) tally + O(n⁴) strong-cluster
assembly, capped at 200 leaves (a memory guard on the dense triplet tensor). The
sub-cubic (k = 2) and near-quadratic algorithms of `Jansson2013a` / `Jansson2016a`
(Apresjan-cluster hierarchy) are a future speed optimisation.
