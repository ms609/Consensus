# R\* consensus — implementation design

Companion to `rstar-spec.md` (the literature spec). This note records the
*implementation* decisions, especially how R\* relates to the Local C++ core and
the two questions that must go to the user / Jansson before the tree-collapse
behaviour can be finalised.

## Key insight: R\* is polynomial and reuses the Aho-BUILD core, NOT the DP

Local (MinRLC/MinILC) needed the exponential `minRILC` subset DP. **R\* does
not.** R\* is just:

1. **Tally** every rooted triplet's state (`ab|c`, `ac|b`, `bc|a`, fan) across
   the k input trees — O(k·n³), the same triplet machinery `get_common_triplets`
   already builds.
2. **Select** the "uniquely favoured" resolution per triplet — see rule below.
3. **Assemble** with Aho et al. BUILD — exactly `build_aho` + `ahoRec` (already
   in `local_consensus.h`, used by `ahoBuild`).

So R\* = (selectable triplet tensor) + `build_aho` + `ahoRec`. It scales to large
n, unlike Local. The shared C++ core to factor out is therefore the **triplet
tensor + Aho-BUILD**, parameterised by the triplet-selection rule:

| Method | Triplet selection rule | Assembly |
|--------|------------------------|----------|
| `aho-build` (FDCT) | resolution present in **all** k trees (unanimous) | BUILD |
| **R\*** | resolution **uniquely favoured** (plurality, see below) | BUILD |
| Local | unanimous (feeds the exponential DP, not BUILD) | minRILC DP |

## Selection rule (from Degnan et al. 2009, p. 36) — PLURALITY, not majority

A resolved triplet `ab|c` enters R\* iff it appears in **strictly more** input
trees than **each** of the other two resolutions `ac|b` and `bc|a`
*considered separately* — i.e. it is *uniquely favoured*. There is **no k/2
threshold**: in a 4-tree 2–1–1 split the 2-vote resolution wins with only 50%.
A tie (no unique winner) leaves those taxa unresolved (a polytomy).

R\* is always a refinement of the majority-rule consensus (every majority clade
appears in R\*), and is a statistically consistent species-tree estimator.

## Rooting

R\* is a **rooted** method (triplets are rooted). Treat input trees as rooted on
their current root, exactly as `Adams()` does. Document this.

## CRITICAL open questions — do NOT guess; route to user / Jansson

> **RESOLVED (see `rstar-findings.md`).** The user supplied the primary source
> `Jansson2016a` (Jansson, Sung, Vu & Yiu 2016). Its §1.1 + Lemma 1.1 settle all
> of these: fans have no impact (OQ1); the R\* tree is the unique tree whose
> clusters are the **strong clusters** of `R_maj`, which always exists — so there
> is no collapse case (OQ2) and no over-resolution latitude (OQ3). `RStar()` now
> uses the strong-cluster construction, not BUILD. The notes below are retained
> for historical context only.

- **OQ1 (fan triplets):** For non-binary input trees, is a fan triplet `(abc)`
  an *abstention* (ignored in the tally) or a *fourth competing state* that can
  itself win the plurality? Degnan et al. assume binary input and do not say.
  This changes the tally arithmetic.
- **OQ2 (BUILD failure / collapse):** When the selected R\* triplet set is
  inconsistent (the Aho graph is a single connected component and BUILD cannot
  split it), the paper says the taxa are "declared unresolved or partially
  unresolved" but gives no exact mechanism. We need the precise collapse rule
  before R\* output is well-defined on conflicting data.

(Full open-question list — OQ1–OQ6 — is in `rstar-spec.md`.)

## Validation strategy — note the gap

There is **no reference binary** for R\* (absent from FACT/FACT2/FDCT). So,
unlike every method so far, R\* cannot be oracle-validated directly. Plan:

- **Property:** R\* clades ⊇ majority-rule clades (`Majority()`); R\* is a
  refinement.
- **Worked example:** the concrete plurality case in `rstar-spec.md` (becomes a
  unit test).
- **Partial cross-check of the BUILD machinery:** the FDCT `aho-build` oracle
  (`LocalOracle(trees, "aho")`) exercises `build_aho`+`ahoRec` on *unanimous*
  triplets — so it validates the assembly half of R\*, just with a different
  selection rule.
- **Hand-built rooted cases** with known triplet tallies.

This validation gap (plus OQ1/OQ2) is why R\* is the riskier of the two heavy
methods and why the open questions should reach the authors.

## Implementation plan (after Local lands)

1. In C++, factor the triplet tally so per-state counts are available, and add a
   plurality selector producing the R\* triplet set; reuse `build_aho`+`ahoRec`.
2. `// [[Rcpp::export]] std::string rStarConsensus(List edgeList, int nTip)` →
   Newick (integer labels), or a signal for the inconsistent case (pending OQ2).
3. R wrapper `RStar(trees)` → rooted `phylo`; trivial guards; rooting per Adams.
4. Tests + property checks above; cite Bryant 2003, Degnan et al. 2009.
