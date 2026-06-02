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
  exceed a method's hard cap (e.g. `RStar` ≤ 200 tips) record `NA`, which is
  itself informative.
- `compare.R` flags any cell where the output **split count** changed. For the
  unique-output methods that is a bug; for **Greedy** it may be the documented
  tie-break on equal-frequency incompatible splits — confirm against the FACT
  oracle and sign off, do not silently re-baseline.

## Stale-build footgun and the build-identity guard

`R CMD INSTALL` copies `.R` sources and `DESCRIPTION` before compiling `src/`.
When `make` finds stale `.o`/`.dll` artefacts (timestamps newer than sources, or
a DLL lock that prevented the previous `INSTALL` from replacing the file), it
prints "Nothing to be done for 'all'" and the compiled library is not updated —
yet `INSTALL` still prints `* DONE`.  The result is a build where the R wrappers
and version are current but the compiled code is old: the oracle and profiling
scripts then silently validate the **wrong** implementation.

`dev/oracle/build-identity.R` addresses this with three layers checked at the top
of every `check-oracle.R`, `baseline.R`, and `compare.R` run:

1. **Version** — `packageVersion("ConsTree")` must equal `DESCRIPTION`'s
   `Version:` field.
2. **Body** — each fast-path R wrapper must reference its C++ symbol
   (e.g. `Greedy()` body must contain `greedyConsensusCpp`).  Catches a stale
   R source or wrong checkout.
3. **Idempotence** — `Method(t, t, t)` on a fixed 12-tip tree must return a
   tree with the same split count as `t`.  Catches a stale DLL where (1) and
   (2) both pass but the compiled algorithm is a stub or from a prior version.

### Force a clean rebuild (Windows / PowerShell)

```powershell
# Clear stale objects if the DLL is still locked by a previous R session:
Remove-Item src\*.o, src\*.dll -ErrorAction SilentlyContinue
R.exe CMD INSTALL --preclean --library=.agent-cons .
```

`--preclean` deletes all `src/*.o` before compiling, ensuring a full recompile
even when timestamps mislead `make`.  Do this whenever the guard aborts with a
"Stale build" message or when you suspect the DLL was not updated.
