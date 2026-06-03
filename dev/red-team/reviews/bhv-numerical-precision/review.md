# Red-team review: BHV geodesic numerical precision
**Date:** 2026-06-03
**Reviewer:** numerical-precision red team (focus: `src/bhv.cpp` GTP, ratio sequence, vertex cover)
**Verdict:** SHIP

The GTP geodesic, ratio-sequence handling, min-weight vertex cover, cone path and
Fréchet iteration are **numerically sound on valid inputs**. Every precision
hunt-target was checked against the source and confirmed with installed-library
repros (`.agent-cons`; not `load_all`). Distances reproduce closed forms to
machine precision (relerr 0–2e-16) across 500+ random many-orthant pairs, across
mixed edge-length scales spanning 20 decades, and the metric invariants
(symmetry, triangle) hold far tighter than the test tolerances claim. The single
real finding is a cosmetic-scale TOL effect that only bites at absurd downscaling
(< ~1e-10 absolute edge length); it is MINOR and documentable, not a correctness
bug. The prior review's F1–F3 fixes are intact and not re-litigated here.

Repro scripts: `dev/red-team/reviews/_bhv_repro.R`, `_bhv_repro2.R` (underscore =
ephemeral per dev protocol; delete-safe).

---

## F1 — Absolute drop-tolerance makes the metric scale-dependent at tiny edge lengths (MINOR)
Where: `src/bhv.cpp:24` (`static const double TOL = 1e-10`) and `:310`
(`if (len[s] <= TOL) continue;` in `from_r`), echoed at `:277,287,291` in
`tree_at`.

BHV distance is homogeneous of degree 1 — `d(sA, sB) = s·d(A,B)` — but the split
drop test is an **absolute** threshold. Any interior edge with `len <= 1e-10` is
deleted as "absent" before the geodesic is built, regardless of the tree's
overall length scale. So a uniformly *downscaled* tree silently loses real
topological structure.

Evidence (repro 1, `op_T`/`op_Tp`, interior lengths {3,4,10}, exact answer
`s·15√2`):

| scale s | distance | expected `s·15√2` | rel err |
|--------:|---------:|------------------:|--------:|
| 1e+00 … 1e-10 | tracks exactly | — | 0 |
| 1e-11   | **0.000e+00** | 2.121e-10 | **1.0** |
| 1e-13   | **0.000e+00** | 2.121e-12 | **1.0** |

At `s = 1e-11` the interior edges become {3e-11, 4e-11, 1e-10}, all `<= TOL`, so
both trees collapse to stars and the distance is reported as exactly 0 instead of
`2.12e-10`. The boundary is precisely the absolute `1e-10`.

Severity rationale: the failure scale (interior branch lengths below ~1e-10) is
pathological for phylogenetics — substitution-model branch lengths are O(1e-3 … 1e1).
Realistic downscaling (1e-3, 1e-6, even 1e-9) tracks the closed form to relerr 0
(repro 1). A genuinely short edge of 1e-8 is **kept** and handled correctly
(repro 1b: `d = 2.121e-07`, exact). So this is a documentable edge-of-domain
property, not a bug on plausible inputs — hence MINOR, not MAJOR.

Recommendation (optional, non-blocking): if scale-robustness is wanted, make the
drop relative — e.g. drop `len <= TOL * scale` where `scale` is the max edge
length across both trees (or the L2 norm of all lengths) — and document the
domain either way. Note this same TOL governs F1's fix from the prior review (a
length-0 interior edge is dropped, which is what makes `sA == 0` in
`min_weight_vc` unreachable); a relative TOL must preserve that exact-zero drop,
so keep an additional `len <= 0` / absolute floor guard if the threshold is made
relative.

---

## F2 — Fréchet `converged` flag tracks step size, not distance-to-mean (INFO)
Where: `src/bhv.cpp:417-421` (`stepDist = lambda * geo_dist(g)`, `lambda = 1/(i+1)`).

The convergence test fires when `cauchyLength` consecutive steps each move less
than `epsilon`. Because `lambda = 1/(i+1) → 0`, the step length shrinks toward 0
*independently of how close the estimate is to the true mean*. So `converged =
TRUE` means "the Sturm walk's steps got small", not "we are within epsilon of the
Fréchet mean". This is already correctly documented in `R/BHV.R:293-297`
("not a guaranteed bound on the distance to the exact mean; tighten `tolerance`").

Confirmed harmless / honest (repro2): at `tolerance = 1e-4` it converges
(`iterations = 10982`); at `tolerance = 1e-12, maxIter = 5000` it hits the cap and
honestly returns `converged = FALSE`. It never cycles or hangs — the step shrinks
monotonically in expectation and the loop is bounded by `maxIter`. No float
rounding can make it loop. INFO only; flagged so a future reader doesn't mistake
the flag for an error bound.

---

## Verified numerically stable

**Geodesic length accumulation — no cancellation (T1, T5).**
`geo_dist` (`:259-266`) accumulates `s += v*v` (ratio terms, `v = eNorm + fNorm >= 0`)
and `s += d*d` (common-edge and leaf terms), then a single `sqrt`. Every summand
is **non-negative**, so there is no subtraction of near-equal quantities and
therefore no catastrophic cancellation — only benign O(n·ε) accumulation rounding.
The only subtractions (`c.lenA - c.lenB`, `leafA[i] - leafB[i]`) occur *before*
squaring; by Sterbenz's lemma those differences are computed exactly when the
operands are within a factor of two, and otherwise contribute their honest
rounded difference. The cone/shared-orthant Euclidean formula `sqrt(Σ(la-lb)²)`
is the same pattern and is likewise clean. Empirically the same-topology test
(`sqrt(nEdge)`) and mixed-scale cones are exact to 0 ulp (repro 1, 2).

**Ratio "sequence" cannot loop or be reordered into error (T2).**
The `std::sort(g.ratios …)` at `:254-255` is **cosmetic**: nothing consumes the
order. `geo_dist` sums ratios order-independently; `tree_at` (`:279-297`) processes
each ratio independently via its own crossing time `time = eN/(eN+fN)`, which is in
[0,1] by construction (both norms ≥ 0). FP rounding that reorders two near-equal
ratios therefore cannot cause an orthant revisit or loop — there is no stateful
walk over the sorted list. (The task's "ratios in [0,1]" conflates the *crossing
time*, guaranteed in [0,1], with the *ratio value* `e/f ∈ [0,∞]`, which is fine
and used only as a sort key.) The GTP partition loop in `gtp_no_common`
(`:165-199`) terminates structurally for positive lengths: a split happens only
when `cov.weight < 1 - TOL`; an empty side would force cover weight exactly 1
(caught as `finalRatio`), so both `r1` and `r2` are non-empty and every split
strictly shrinks the subproblem. Empirical additivity `d(A,M)+d(M,B)=d(A,B)` over
500 random **12-tip** pairs (deep vertex-cover recursion): max rel err
**2.45e-16** (repro2).

**Min-weight vertex cover is well-conditioned across extreme scale spread (T3).**
Capacities are *normalized*: `cap[S→a_i] = wA[i]/sA`, `cap[b_j→T] = wB[j]/sB`, so
each side sums to 1 and the augmenting-path threshold `eps = 1e-13` (`:77`) is a
**relative** tolerance — scale-invariant by design. `sA = Σ wA[i]` is a sum of
non-negative squared lengths (`:73-74`), no cancellation. A cone pair with
interior lengths `{1e-10, 1e10}` returns the exact cone distance `1e10` (relerr 0,
repro 1, both orderings). An 8-tip incompatible pair with **every** interior
length multiplied by `10^U(-10,10)` (a 20-decade spread within single vertex-cover
groups) gives additivity rel err **1.24e-16** (repro2). The only theoretical
exposure — a sub-`1e-13` *normalized* weight becoming invisible to the max-flow —
costs at most an error of that same negligible relative order, and was not
observed. Not a finding.

**Near-zero legitimate edges are preserved (T4 happy path).**
An edge of `1e-8` is well above `TOL = 1e-10` and is kept and handled exactly
(repro 1b). Only the pathological sub-1e-10 regime in F1 drops real structure.

**Fréchet iteration converges and never cycles (T6).**
Covered in F2: bounded by `maxIter`, step `1/(i+1) → 0`, honest `converged` flag,
no float-induced cycling. The pre-pass epsilon scaling (`:396-404`) guards the
identical-trees degenerate case (`epsilon = tol` when stdev is 0).

**Metric invariants hold far inside the test tolerances (T7).**
- Symmetry: max `|d(a,b) − d(b,a)|` over 2000 random 8-tip pairs = **4.44e-16**
  (1 ulp). The `1e-9`/`1e-2` test slacks are orders of magnitude looser than the
  true precision, so they are **conservative, not concealing** — a real precision
  defect producing violations ≫ 1e-9 would still fail the test.
- Triangle: max `d(t1,t3) − d(t1,t2) − d(t2,t3)` over 2000 random 8-tip triples =
  **−1.08** (always strongly negative; never approaches 0⁺). The `+1e-9` slack in
  `test-BHV.R:54` is never exercised.
- Recommendation (optional): the triangle/symmetry tolerances could be tightened
  to ~1e-12 to make the suite a tighter precision tripwire, but this is polish,
  not a requirement.

---

## Scope notes
- Did **not** re-report F1 (zero-length-edge loop), F2 (singleton split), F3
  (missing `edge.length`) from `claude-cranky-kare` — confirmed fixed and intact.
  The prior F1 fix (the `len <= TOL` drop in `from_r`) is what makes
  `min_weight_vc`'s `sA == 0` / `sB == 0` division unreachable; verified by the
  zero-length-edge test still passing and by the TOL repro.
- Build identity: ran against the pre-installed `.agent-cons` library in the
  shared checkout; git-tracked `src/`/`R/` source is unmodified, so the installed
  binary matches the source under review (no edits made).
