# ConsTree — Agent Instructions

## Project overview

**ConsTree** (`ConsTree`, GPL ≥ 3) is the front-end consensus-tree toolkit built
on [**TreeTools**](https://ms609.github.io/TreeTools/), which remains the engine
for the methods it wraps (strict, majority). Several split- and triplet-based
methods are R/C++ ports of the asymptotically-efficient C++ in
[FACT](https://github.com/Mesh89/FACT),
[FACT2](https://github.com/Mesh89/FACT2) and
[FDCT_new](https://github.com/tswddd2/FDCT_new) (Jesper Jansson et al.; used with
permission). The committed FACT reference source lives in `dev/oracle/fact-src/`
— read it when porting or validating; never ship it. (FACT2 / FDCT_new are
cloned ad-hoc into `dev/reference/` when needed, e.g. for the Frequency port.)

## Repository layout

```
ConsTree/
├── R/
│   ├── ConsTree-package.R   # package doc + Rdpack glue
│   ├── wrappers.R           # Strict(), Majority(), MajorityRule() → TreeTools::Consensus
│   ├── selection.R          # split-selection family + shared .PoolSplits machinery
│   ├── adams.R              # Adams()  (refinement to fixpoint)
│   ├── local.R              # Local()  → src/local_consensus.cpp
│   ├── rstar.R              # RStar()  → src/rstar.cpp
│   ├── Quartet.R            # Quartet() → src/Quartet.cpp
│   ├── Average.R            # distance-based average consensus (Lapointe & Cucumel)
│   ├── BHV.R                # BHV geodesic distance / Fréchet mean → src/bhv.cpp
│   └── RcppExports.R        # generated; never edit
├── src/                     # Rcpp C++ (no Boost in shipped code)
├── tests/testthat/          # property + idempotence + ecosystem cross-checks
├── vignettes/ConsTree.Rmd   # introductory worked-examples vignette
├── inst/REFERENCES.bib      # Rdpack bibliography (cite via \insertCite{Key}{ConsTree})
├── inst/WORDLIST            # spelling::spell_check_package() allow-list
├── man/                     # roxygen-generated; never edit by hand
└── dev/                     # NOT shipped (.Rbuildignore: ^dev$); oracle + reference C++
```

## Core architecture: the C++ fast-path pipeline

Every split-based method in `R/selection.R` follows the same three-step
pattern — only the C++ selection logic differs:

1. **`.PrepareTrees(trees)`** — validates and normalises the input, returning a
   list with `trees`, `labels`, `nTree`, `firstTree`, and a `trivial` shortcut
   for degenerate cases (single `phylo`, < 2 trees, < 4 leaves).
2. **`.FactEdges(trees, labels)`** — renumbers each tree to the shared `labels`
   order, roots at `labels[[1]]`, and returns the `Preorder` edge matrices that
   the C++ entry points consume.
3. **C++ call** — `looseConsensusCpp` / `greedyConsensusCpp` /
   `majorityPlusConsensusCpp` / `frequencyConsensusCpp` — returns a Newick
   string; tip labels are decoded from integer indices, then
4. **`.RootLikeFirst(tree, firstTree)`** roots the result to match the first
   input tree, as `TreeTools::Consensus()` does.

## Method roster

| `ConsTree::` | File | Source | Status / limits |
|---|---|---|---|
| `Strict()` | wrappers.R | TreeTools | FACT-exact |
| `Majority()` / `MajorityRule()` | wrappers.R | TreeTools | FACT-exact |
| `Loose()` | selection.R | FACT `loose.cpp` | FACT-exact |
| `Greedy()` | selection.R | FACT `greedy.cpp` | FACT-match up to tie-break (documented) |
| `MajorityPlus()` | selection.R | FACT `majorityplus.cpp` | FACT-exact |
| `Frequency()` | selection.R | FACT2/FDCT `freqdiff.h` | property-validated (`dev/oracle/freqdiff/`) |
| `Adams()` | adams.R (Rcpp) | Jansson, Li & Sung 2017 (`src/cons_adams.cpp`) | slow-Adams clade-oracle exact |
| `Local(type=)` | local.R (Rcpp) | FACT2/FDCT `local_consensus.h` (Jansson, Rajaby & Sung 2018) | FDCT-oracle exact; **≤20 leaves**; see runtime caveat |
| `RStar()` | rstar.R (Rcpp) | Jansson et al. 2016 | definition-exact (strong-cluster oracle); **≤200 leaves** (memory) |
| `Quartet()` | Quartet.R (Rcpp) | Quartet pkg (ported) | brute-force oracle (n=5); ≤100 tips |
| `Average()` | Average.R | Lapointe & Cucumel | user-authored (path-length LS / BME) |
| `BHV…()` | BHV.R (Rcpp) | BHV geodesic | distance / Fréchet-mean utilities |

**`Local()` runtime caveat.** Exact-exponential. The ≤20-leaf guard bounds
*memory*, but runtime also depends on input congruence — incongruent
`type = "induced"` can be intractable even at n ≤ 20. The call is interruptible
(Ctrl-C).

## Open decisions (unresolved — surface to the user, don't silently pick)
- **`Local()` hard runtime guard**: none yet (only interruptible) — decide
  whether to add one.

## Resolved
- **`Local()` citation**: now `\insertCite{JanssonRajabySung2018}` (Jansson,
  Rajaby & Sung 2018, *AIMS Medical Science* 5(2):181–203,
  doi:10.3934/medsci.2018.2.181), the paper that defines the MinRLC/MinILC
  variants `Local()` implements.

## The consensus lattice (invariants the tests enforce)

```
strict ⊆ loose  ⊆ greedy
strict ⊆ majority ⊆ greedy
majority ⊆ majorityPlus
majority ⊆ frequency ⊆ greedy
majority ∥ loose            (incomparable)
```

`tests/testthat/test-selection.R` asserts these via order-independent split-string
sets (`splitSet()` helper). Idempotence — `Method(list(t, t, t))` returns `t` —
holds for every method.

## Code conventions (house style)

- **Functions** `BigCamelCase`; **variables** `camelCase`; **private helpers**
  dot-prefixed (`.PoolSplits`). One exported method per concept; keep the public
  API reachable as `ConsTree::<Method>()`.
- **`@importFrom` per roxygen block** — never a blanket `@import`.
- **`# Return:` comment** marks the returned value instead of explicit `return()`;
  early `return()` for guard clauses is fine (TreeTools precedent).
- **Prefer TreeTools over ape** for tree manipulation (`NTip`, `Preorder`,
  `RootTree`, `KeepTip`, `RenumberTips`, …); use `ape::` only where there is no
  equivalent (`ape::as.phylo`, `ape::read.tree`, distance/NJ builders).
- **Full test coverage** of every new/changed line — happy path, error branches,
  edge cases (single tree, < 4 leaves, identical trees, conflicting splits).
  codecov gates the PR.
- **Cite** with `\insertCite{Key}{ConsTree}` + `@references \insertAllCited{}`;
  add BibTeX to `inst/REFERENCES.bib` (clean `AuthorYear` keys). Add legitimate
  technical terms flagged by `spelling::spell_check_package()` to `inst/WORDLIST`.
- After each user-visible change: update `NEWS.md` and bump the `.900X` dev suffix
  in `DESCRIPTION`.

## Build / test / docs

- **Do not `devtools::load_all()`** for validation — install to an isolated
  library and test against it (avoids DLL/namespace surprises; matches how the
  oracle driver loads the package).
- PowerShell aliases `R` to `Invoke-History`; **always call `R.exe` / `Rscript.exe`**.

```bash
# Install to an isolated library (from the package root)
R.exe CMD INSTALL --library=.agent-cons .

# Regenerate Rd after any roxygen change
Rscript.exe -e "roxygen2::roxygenise()"

# Run the suite / a single file against the isolated library
Rscript.exe -e ".libPaths(c('.agent-cons', .libPaths())); testthat::test_dir('tests/testthat')"
Rscript.exe -e ".libPaths(c('.agent-cons', .libPaths())); testthat::test_file('tests/testthat/test-selection.R')"
```

A task is complete only when the relevant tests pass and (for release) `R CMD
check` is clean. Prefer GHA for full `R CMD check` (see `../AGENTS.md`). One
accepted NOTE: `read_symbols_from_dll` (R-devel tooling, not a package bug).

The phangorn `allCompat` cross-check is gated behind `CONSENSUS_PHANGORN_TESTS=1`
— R-devel/phangorn ABI mismatch can **SIGSEGV** (which `tryCatch()` cannot catch),
so it must be opt-in rather than wrapped.

## Dev directory protocol

`dev/` is excluded from the shipped R package (`.Rbuildignore: ^dev$`) and must
never be imported or loaded by package code.  Git tracks everything in `dev/`
**except** files whose names begin with `_` — those are ephemeral artefacts
(coverage logs, error dumps, scratch NEXUS fixtures, etc.) and are `.gitignore`d
via the pattern `dev/**/_*`.

**Convention:** if you create a file in `dev/` that you would not want to open in
a future session, prefix its name with `_`.  Scripts, notes, oracle sources, and
oracle binaries are tracked; build logs, intermediate outputs, and debug dumps are
not.

Agent worktrees (`git worktree add`) will contain everything git-tracked in `dev/`,
including the oracle binaries — the oracle can therefore be run from any worktree.

## Dev oracle (reference-grade validation)

`dev/oracle/` shells out to patched **FACT/FACT2/FDCT** binaries to diff each
ported method against the original C++ on shared random fixtures — the strongest
check available.

- `dev/oracle/oracle.R` — `FactConsensus(trees, method, rooted)` writes FACT-dialect
  NEXUS (integer labels, literal `translate` line), invokes `fact.exe` over stdin
  (`file rule rooted`), and parses integer-labelled Newick back to `phylo`.
  `FACT_RULE` maps method → algorithm bitmask:
  `strict=1, majority=8, greedy=32, loose=128, majorityPlus=256, adams=512`
  (the *slow*, classical Adams; the fast bit `1024` mis-prints small trees, so it
  is not used as an oracle).
- `dev/oracle/check-oracle.R` — cross-validation driver
  (`Rscript.exe dev/oracle/check-oracle.R`). Per-method harnesses also live under
  `dev/oracle/{freqdiff,local,rstar}/`.
- Oracle binaries (`fact.exe`, `freqdiff/freqdiff.exe`, `local/local.exe`) are
  git-tracked so any worktree can run the oracle without a rebuild.

### Rebuilding FACT oracle binaries (MinGW)
FACT's *generation*-only path uses POSIX bits absent on MinGW; the patched build
drops them: remove `#include <sys/resource.h>` and `set_stack`, map
`rand_r(&seed_)` → `rand()`, `#define M_PI` (compile `-D_USE_MATH_DEFINES`). The
consensus paths are untouched. FACT2/FDCT additionally need
`boost::dynamic_bitset` headers at build time — fine for the dev oracle, never for
the shipped package.
