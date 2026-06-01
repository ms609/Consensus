# Consensus — Agent Instructions

> **Cross-package coordination** (worktree discipline, commit attribution, GHA
> dispatch, CPU limits, multi-agent protocol) lives in
> [`../AGENTS.md`](../AGENTS.md). That file is authoritative for anything that
> spans packages. This file covers conventions specific to **ConsTree**.

## Project overview

**ConsTree** is an R package (GPL ≥ 3) providing a comprehensive suite of
phylogenetic consensus-tree methods, built on
[**TreeTools**](https://ms609.github.io/TreeTools/). It is the *front-end*
consensus toolkit; TreeTools remains the fast engine for the two methods it
already ships (strict, majority).

The package repackages a family of asymptotically-efficient algorithms
prototyped in C++ by Jesper Jansson and colleagues
([FACT](https://github.com/Mesh89/FACT),
[FACT2](https://github.com/Mesh89/FACT2),
[FDCT_new](https://github.com/tswddd2/FDCT_new)). We **have permission to
incorporate the source** — port directly, attribute in `DESCRIPTION`
(`Copyright:`) and cite the source papers (`inst/REFERENCES.bib`).

Not yet listed in the `../AGENTS.md` package registry and not (yet) part of the
multi-agent worktree rotation; develop on `main` or a feature branch as the user
directs.

## Two-home architecture (read before adding anything)

A decision taken with Martin (2026-05-31) governs *where* code lands:

- **Dependency-free upgrades to functionality that already exists in TreeTools**
  → go into **TreeTools** (keeps the foundation fast; everything downstream
  inherits the speed-up). Example: the O(kn) majority-consensus core upgrade,
  which runs as a *separate* TreeTools chip task — **do not** duplicate it here.
- **New methods** → go into **this package**. `Strict()` / `Majority()` are
  here only as thin, consistently-named wrappers around `TreeTools::Consensus()`.

**No `BH` (Boost) dependency.** Boost is a huge download and is unavailable on
webR. Audit confirmed the only Boost use in the reference repos is
`boost::dynamic_bitset` in `freqdiff.h`; it is replaceable by TreeTools' packed
splits / `DynamicBitset`. **Performance is king** — reimplement, don't sacrifice
speed, but ship nothing that pulls in Boost. (The dev *oracle* binaries below may
link Boost — they are never part of the package.)

## Repository layout

```
Consensus/
├── R/
│   ├── Consensus-package.R   # package doc + Rdpack glue
│   ├── wrappers.R            # Strict(), Majority(), MajorityRule() → TreeTools::Consensus
│   ├── selection.R           # split-selection family + shared .PoolSplits machinery
│   └── Average.R             # distance-based average consensus (Lapointe & Cucumel)
├── tests/testthat/           # property + idempotence + ecosystem cross-checks
├── inst/REFERENCES.bib       # Rdpack bibliography (cite via \insertCite{Key}{Consensus})
├── dev/                      # NOT shipped (.Rbuildignore: ^dev$); see "Dev oracle"
│   ├── oracle/               # FACT binary bridge + cross-validation driver
│   └── reference/            # cloned FACT / FACT2 / FDCT_new (read, don't ship)
└── man/                      # roxygen-generated; never edit by hand
```

`dev/` is `.Rbuildignore`d and untracked by git — so it is **absent from any
git worktree**. A subagent in an isolated worktree cannot reach the oracle
binary; run oracle validation from the main checkout (or pass results in).

## Core architecture: the split-pooling pipeline

Every split-based method shares one validated pipeline, implemented in
`R/selection.R`. New split-based methods should reuse it rather than re-deriving
it — only the *selection rule* differs.

1. **`.PoolSplits(trees)`** pools the bipartitions across all trees and returns a
   list with the contract:
   - `splits` — the **distinct** splits (a `Splits` object);
   - `members` — `as.logical(distinct)`, an (nSplit × nTip) logical matrix;
   - `counts` — occurrence count of each distinct split (`integer`);
   - `membership[[i]]` — indices (into `splits`) of the splits in tree *i*;
   - `labels`, `nTree`, `firstTree`;
   - `trivial` — a short-circuit tree for the degenerate cases (a single
     `phylo`, fewer than two trees, or fewer than four leaves). When
     `trivial` is non-`NULL`, return it immediately.
2. **`.CompatibilityMatrix(prep)`** — the pairwise split-compatibility matrix,
   `as.matrix(CompatibleSplits(prep$splits, prep$splits))`. Its logical negation
   is the *conflict* matrix used by `MajorityPlus()` / `Frequency()`.
3. **Selection** — a method-specific rule produces a logical `keep` over
   `prep$splits` (see each method in `selection.R`).
4. **`.SelectedConsensus(keep, prep)`** — rebuild via
   `as.phylo(as.Splits(members[keep, ], tipLabels = labels))` and root to match
   the first input tree (the same rooting `TreeTools::Consensus()` uses).

This pattern is **pure R orchestration over TreeTools' C++ primitives** — no
compiled code in this package (yet), hence webR-clean. The O(kn)/O(kn log n)
optima can follow later as compiled enhancements once the TreeTools core upgrade
lands; correctness-first R implementations come first and are the oracle
reference.

### Gotcha: `Splits` is an S4 class

Single-bracket `[` subsetting of a `Splits` object **drops the class**. To select
a subset of splits, round-trip through the logical matrix:
`as.Splits(prep$members[keep, , drop = FALSE], tipLabels = prep$labels)`. Never
`splits[keep]`.

### Gotcha: `as.phylo` is the **ape** generic

TreeTools registers S3 methods but does not export `as.phylo`. Use
`ape::as.phylo(...)` in code and examples (it dispatches to the TreeTools method
when a `Splits`/TreeTools object is passed and the TreeTools namespace is
loaded).

## Method roster & status

| `Consensus::` | Kind | Source | Status |
|---|---|---|---|
| `Strict()` | wrapper → `Consensus(p=1)` | TreeTools | ✅ FACT-exact |
| `Majority()` / `MajorityRule()` | wrapper → `Consensus(p)` | TreeTools | ✅ FACT-exact |
| `Loose()` | split-selection (compatible-with-all) | FACT `loose.cpp` | ✅ FACT-exact |
| `Greedy()` | split-selection (count-sorted accretion) | FACT `greedy.cpp` | ✅ FACT-match up to tie-break |
| `MajorityPlus()` | split-selection (present > incompatible) | FACT `majorityplus.cpp` | ✅ FACT-exact |
| `Frequency()` | split-selection (freq > every rival) | FACT2/FDCT `freqdiff.h` | ✅ property-validated; Boost oracle pending |
| `Average()` | distance-based (path-length LS / BME) | Lapointe & Cucumel | ✅ (user-authored) |
| `Adams()` | refinement to fixpoint (**not** split-selection) | FACT `adams.cpp` | ✅ FACT clade-oracle (exact) |
| `Local(type=)` | triplet + Aho-BUILD + 2^n DP (MinRLC/MinILC) | FACT2/FDCT `local_consensus.h` | ✅ compiled C++ (Rcpp); FDCT oracle exact (38); ≤20 leaves |
| `RStar()` | plurality triplets → strong clusters (Lemma 1.1; polynomial) | Jansson et al. 2016 (def.) | ✅ compiled C++ (Rcpp); definition-exact (brute-force strong-cluster oracle 60/60) + majority refinement + congruent `aho-build` oracle; ≤200 leaves (memory) |
| `Quartet()` | greedy add-and-prune; minimizes sum of symmetric quartet distances (Takazawa et al. 2026) | Quartet pkg (ported) | ✅ compiled C++ (Rcpp); brute-force oracle (n=5); ≤100 tips |

**Status (resume anchor).** Adams, Local, and **R\*** are all done and validated.
The package contains compiled C++ (Rcpp, no Boost — `src/local_consensus.cpp`,
`src/rstar.cpp`); Frequency is reference-validated (`dev/oracle/freqdiff/`).
**R\* is implemented and definition-exact** (own `src/rstar.cpp`: dense triplet
tally → strict plurality selection → **strong-cluster** assembly per Lemma 1.1 of
Jansson, Sung, Vu & Yiu 2016; polynomial, `std::vector`+DSU so it scales — ≤200
leaves as a *memory* guard, distinct from Local's *algorithmic* ≤20). The former
"open questions" are **resolved by the primary source** (`Jansson2016a`, supplied
by the user): fans have no impact (OQ1); the strong-cluster tree always exists and
is unique, so there is no collapse case (OQ2) and no over-resolution latitude
(OQ3). An earlier Aho-BUILD assembly was found to violate the definition
(`r(τ) ⊆ R_maj`) on conflicting input — under-resolving (~46/400 trials) and
over-resolving (~51/400) — and was replaced; see `dev/notes/rstar-findings.md`
for the before/after. Anchored by: a **brute-force strong-cluster oracle**
(`RStar` clades == strong clusters over all 2ⁿ subsets, independent `ape::mrca`
tally, 60/60), **majority-rule refinement** (every rooted majority clade present;
0/400, was 46/400), identity (`RStar(k×T)==T`, n≤40), the FDCT `aho-build` oracle
on congruent input (6/6), and strict-clade refinement. Harness in
`dev/oracle/rstar/` (`check-rstar.R`, `check-strong-clusters.R`); tests in
`tests/testthat/test-rstar.R`. The sub-cubic/Apresjan O(n²) algorithm of
`Jansson2013a`/`Jansson2016a` is a deferred speed optimisation.
Local is exact-exponential: the ≤20-leaf guard bounds memory, but runtime also
depends on input congruence (incongruent `type="induced"` can be intractable
even at n≤20); the call is interruptible. The phangorn `allCompat` cross-check is
gated behind `CONSENSUS_PHANGORN_TESTS=1` (R-devel/phangorn ABI can segfault).
`R CMD check` is clean bar one R-devel tooling NOTE (`read_symbols_from_dll`).
Open user decisions: the `Local()` citation (placeholder `JanssonShenSung2016` is
likely wrong for MinRLC/MinILC); R\* OQ1/OQ2/OQ3 (esp. the OQ2 collapse mechanism
above); whether Local needs a hard runtime guard. Next non-blocked work:
`Quartet()` migration + a vignette (Phase 4).

## The consensus lattice (invariants the tests enforce)

```
strict ⊆ loose  ⊆ greedy
strict ⊆ majority ⊆ greedy
majority ⊆ majorityPlus
majority ⊆ frequency ⊆ greedy
majority ∥ loose            (incomparable — neither always contains the other)
```

Property tests in `tests/testthat/test-selection.R` assert these via
order-independent split-string sets (`splitSet()` helper). Idempotence —
`Method(list(t, t, t))` returns `t` — holds for every method.

## Code conventions (house style)

- **Functions** `BigCamelCase`; **variables** `camelCase`; **private helpers**
  dot-prefixed (`.PoolSplits`). One exported method per concept; keep public API
  reachable as `Consensus::<Method>()`.
- **`@importFrom` in each function's roxygen block** — never a blanket `@import`.
- **`# Return:` comment** marks the final returned value instead of an explicit
  `return()`; early `return()` for guard clauses is fine (TreeTools precedent).
- **Prefer TreeTools over ape** for tree manipulation (`NTip`, `Preorder`,
  `RootTree`, `KeepTip`, `RenumberTips`, …); use `ape::` only where there is no
  TreeTools equivalent (`ape::as.phylo`, `ape::read.tree`, distance/NJ builders).
- **Full test coverage** of every new/changed line — happy path, error branches,
  edge cases (single tree, < 4 leaves, identical trees, conflicting splits).
  codecov gates the PR.
- **Cite** with `\insertCite{Key}{Consensus}` + `@references \insertAllCited{}`;
  add the BibTeX to `inst/REFERENCES.bib` (clean `AuthorYear` keys, matching the
  existing entries). Add legitimate technical terms flagged by
  `spelling::spell_check_package()` to `inst/WORDLIST`.
- After each user-visible change: update `NEWS.md` and bump the `.900X` dev
  suffix in `DESCRIPTION` (ecosystem rule).

## Build / test / docs

**Do not use `devtools::load_all()`** for validation here — install to an
isolated library and test against it (avoids DLL/namespace surprises and matches
how the oracle driver loads the package).

PowerShell aliases `R` to `Invoke-History`; **always call `R.exe`/`Rscript.exe`**
explicitly.

```bash
# Install to an isolated library (from the package root)
R.exe CMD INSTALL --library=.agent-cons .

# Regenerate Rd after any roxygen change (original package name required)
Rscript.exe -e "roxygen2::roxygenise()"

# Run the test suite against the isolated library
Rscript.exe -e ".libPaths(c('.agent-cons', .libPaths())); testthat::test_dir('tests/testthat')"

# Run a single test file
Rscript.exe -e ".libPaths(c('.agent-cons', .libPaths())); testthat::test_file('tests/testthat/test-selection.R')"
```

A task is complete only when the relevant tests pass and (for release) R CMD
check is clean. Prefer GHA for full R CMD check (see `../AGENTS.md`).

## Dev oracle (reference-grade validation)

`dev/oracle/` shells out to a patched **FACT** binary (`fact.exe`) to diff each
ported method against the original C++ on shared random fixtures — the strongest
check we have (we have permission to use the source).

- `dev/oracle/oracle.R` — `FactConsensus(trees, method, rooted)` writes a
  FACT-dialect NEXUS (integer taxon labels; literal `translate` line; picky
  header), invokes `fact.exe` over stdin (`file rule rooted`), and parses the
  integer-labelled Newick back to a `phylo` with original labels.
  `FACT_RULE` maps method → algorithm bitmask
  (`strict=1, majority=8, greedy=32, loose=128, majorityPlus=256, adams=1024`).
- `dev/oracle/check-oracle.R` — the cross-validation driver. Run with
  `Rscript.exe dev/oracle/check-oracle.R`.

**Confirmed:** strict / majority / loose / majorityPlus match FACT exactly;
greedy matches up to tie-break (equal-frequency conflicting splits — an
implementation-defined choice, documented in the roxygen).

### Building FACT oracle binaries (MinGW notes)

FACT's *generation*-only path uses POSIX/Unix bits absent on MinGW. The patched
build copy drops them: remove `#include <sys/resource.h>` and the `set_stack`
function, map `rand_r(&seed_)` → `rand()`, `#define M_PI` (and compile with
`-D_USE_MATH_DEFINES`). The consensus paths themselves are untouched. FACT2 /
FDCT_new additionally need `boost::dynamic_bitset` headers at build time — fine
for a dev oracle, never for the shipped package.

## Known environment gotchas

| Symptom | Cause | Fix |
|---|---|---|
| `R` runs `Invoke-History` | PowerShell alias | call `R.exe` / `Rscript.exe` |
| `as.phylo not exported from TreeTools` | it's the ape generic | use `ape::as.phylo` |
| `as.phylo.default(...)` on a subset | `[` dropped the `Splits` S4 class | logical-matrix round-trip |
| `INTEGER() ... not 'integer'` from phangorn | R-devel/phangorn ABI mismatch in this env | gate the cross-check behind an opt-in env var (`CONSENSUS_PHANGORN_TESTS=1`); **not** `tryCatch()` — the failure can be a SIGSEGV, which `tryCatch()` cannot catch |
| oracle binary missing | `dev/` is untracked → absent in worktrees | run oracle from the main checkout |

