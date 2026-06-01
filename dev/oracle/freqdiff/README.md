# dev/oracle/freqdiff — frequency-difference oracle

Dev-only tooling.  Not part of the package (covered by `.Rbuildignore`).
Cross-validates `ConsTree::Frequency()` against the FDCT_new reference
binary by Jesper Jansson et al. (<https://github.com/tswddd2/FDCT_new>).

## Files

| File | Purpose |
|------|---------|
| `freqdiff.exe` | Compiled oracle binary (Windows/MinGW) |
| `oracle_fd.R` | R bridge: `FreqDiffOracle(trees)` and `SplitSetFD()` helpers |
| `check-freqdiff.R` | Cross-validation driver; run with `Rscript.exe dev/oracle/freqdiff/check-freqdiff.R` |
| `README.md` | This file |

No patched source copies were needed; the binary is compiled directly from
`dev/reference/FDCT_new/` (see below).

## Binary I/O format (reverse-engineered)

The binary is compiled with `#define LOCAL_TEST` active in `main_new.cpp`,
which hard-wires the algorithm to `"freq"` (frequency-difference consensus).

- **Input**: `inp.txt` in the process CWD — one Newick per line, **integer
  leaf labels only** (`Tree.cpp` asserts `isdigit(*str)` for every leaf token;
  string labels cause an assertion failure at exit code 3).  Trailing
  semicolons are harmless (the parser skips to the first `(`), but are not
  required.  A blank or short (≤2 char) line terminates parsing.
- **Output**: `oup.txt` in the process CWD — a single Newick line using the
  same integer label strings from the input.  No trailing semicolon; append
  one before passing to `ape::read.tree()`.
- **Stdout**: two timing floats from the weight-calculation and filter passes;
  ignored by the R bridge.

The R bridge (`FreqDiffOracle`) uses `TreeTools::RenumberTips` to align each
tree with the reference label order, replaces tip labels with `1..n` integers,
writes `inp.txt`, invokes the binary, reads `oup.txt`, and maps the integer
labels back to the original names.  A per-PID scratch subdirectory is used so
concurrent calls don't stomp on each other.

## Rebuilding the binary

Requirements:
- g++ from Rtools 45 (`C:/rtools45/X86_64~1.POS/bin/g++.exe`, or on PATH as `g++`)
- BH R package installed (provides `boost::dynamic_bitset` headers)

```bash
# From the Consensus package root
BH_INC=$(Rscript.exe -e "cat(system.file('include', package='BH'))")

g++ -O2 -std=c++17 -D_USE_MATH_DEFINES \
    -I"$BH_INC" \
    dev/reference/FDCT_new/main_new.cpp \
    -o dev/oracle/freqdiff/freqdiff.exe
```

**No source patches were required.**  The only POSIX-specific code in FDCT_new
(`unistd.h` in `databuilder.cpp` and `databuilder_paper.cpp`) is not reachable
from `main_new.cpp`'s include chain.  `<chrono>` is used for timing (not
`gettimeofday`), and no `rand_r`/`sys/resource.h` calls appear in the compiled
translation unit.

## Running the cross-validation

```bash
Rscript.exe dev/oracle/freqdiff/check-freqdiff.R
```

Expected output (confirmed 2026-05-31):

```
-- random  n9  k21 --
  Frequency() mine= 5  freqdiff ref= 5  MATCH

-- random  n10 k31 --
  Frequency() mine= 5  freqdiff ref= 5  MATCH

-- conflict n8  k7 --
  Frequency() mine= 5  freqdiff ref= 5  MATCH

ALL MATCH
```

## Comparison methodology

Frequency-difference is an unrooted split-based method.  Both sides are
compared as unrooted bipartition sets via
`unname(as.character(TreeTools::as.Splits(tree, tipLabels = labels)))`,
using `setequal()` for order-independent equality — the same approach as
`dev/oracle/check-oracle.R`.

If discrepancies appear the script diagnoses whether they are tie-break
artifacts (equal-frequency conflicting splits, where the C++ and R
implementations may make different arbitrary choices) vs genuine algorithmic
differences.
