# dev/profiling — consensus-method benchmarks

Dev-only tooling (not shipped; covered by `.Rbuildignore: ^dev$`). Establishes
speed baselines for the six audited consensus methods and proves the speedups
from each fast-algorithm reimplementation.

## Files

| File | Purpose |
|------|---------|
| `bench-common.R` | shared helpers: `makeTrees()` generator (two regimes), `timeCall()` (median-of-reps, timeout-guarded), `benchGrid()` |
| `baseline.R` | capture pre-change timings of all six methods → `baseline-<date>.csv` |
| `compare.R` | re-run after a change and diff against a saved baseline (speedup + split-count check) |
| `baseline-<date>.csv` | saved baselines (git-tracked) |

## Workflow

```bash
# 1. BEFORE changing a method, capture the baseline (once per dev cycle):
Rscript.exe dev/profiling/baseline.R

# 2. Reimplement the method, then reinstall to the isolated library:
R.exe CMD INSTALL --library=.agent-cons .

# 3. Prove the win (all methods, or name some):
Rscript.exe dev/profiling/compare.R
Rscript.exe dev/profiling/compare.R baseline-2026-06-02.csv Greedy
```

## Notes

- Timing uses base `system.time` (median of 3 reps), not `bench` — avoids a
  dependency for a dev script; expected wins are order-of-magnitude.
- Trees are generated with `TreeTools::RandomTree` (project convention), not
  `ape::rtree`. Two regimes: `independent` (incongruent, large split pool —
  worst case for the `O(s²)` R pipeline) and `perturbed` (mostly congruent).
- Each call is guarded by a per-call elapsed `timeout`; cells that time out or
  exceed a method's bench cap (`BENCH_CAPS`) record `NA`, which is itself
  informative. (`RStar`'s former 200-leaf *memory* cap is gone — see below — so
  its `BENCH_CAPS` entry now only bounds grid runtime.)
- `compare.R` flags any cell where the output **split count** changed. For the
  unique-output methods that is a bug; for **Greedy** it may be the documented
  tie-break on equal-frequency incompatible splits — confirm against the FACT
  oracle and sign off, do not silently re-baseline.
- **RStar** (round 1): the dense `O(n^3)` triplet tensor and its hard 200-leaf
  cap were removed — memory is now `O(kn^2)` via per-tree constant-time LCA — and
  the strong-cluster assembly was tightened from `O(n^4)` to about `O(n^3)`. The
  R\* tree is unchanged (clade-exact vs the previous build on every grid cell to
  n = 200; `dev/oracle/rstar/check-vs-legacy.R`). At n ≤ 200 timing is at parity
  with the former code (the tally is still `O(kn^3)`); the practical win is the
  lifted cap (n = 300 ≈ 0.2 s, n = 500 ≈ 0.7 s at k = 10 — formerly an immediate
  error) and the far lower memory.
