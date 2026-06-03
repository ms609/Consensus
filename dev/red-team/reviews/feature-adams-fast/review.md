# Review: feature/adams-fast @ 7b42a9c (base 9c471a5)

Independent review of the Adams O(kn log n) centroid-path patch (PR #7).
Begun by the external-reviewer subagent; completed in the main session after that
agent hung mid-run (during the oracle execution). All claims below were
re-verified here, not carried over on trust.

## Resolution (F1 fixed, post-review)
F1 was resolved by making the bound *true* rather than restating it. Both
comparison sorts were replaced with linear radix passes in `src/cons_adams.cpp`:
- `buildNub` (was `:190`): stable LSD radix on `firstVisit` (key < eulerLen) — O(m).
- meet-key grouping (was `:402`): LSD radix over the k coordinates with
  epoch-stamped dense coordinate compression — O(k·sz) per spine step.

Re-validated on the rebuilt package (R-devel): test suite 60/60; oracle
clade-exact at scale (counts 7/11/27/38/12 unchanged, "All oracle comparisons
MATCH"); independent brute-force 0/320 binary + 0/350 multifurcating + 10/10
adversarial; timing 0.02 s at n=200/k=50 (no regression vs the prior C++ build),
clean scaling to n=1600. The `O(kn log n)` claim in NEWS/`@details`/header is now
accurate, so no documentation change was required.

Minors F2–F5 were then folded in:
- **F2** — a reachability probe (guards temporarily `Rcpp::stop`, 16 042 fuzz
  cases) showed `kids.size()==1` **is reachable and load-bearing** (degree-1
  suppression at a deepest single-block step; e.g. the n=7/k=2 pair now in
  `test-adams.R`), *not* "defensive only" as the comment claimed — corrected the
  comment and added that regression. `kids.empty()` is genuinely unreachable and
  is now wrapped in `// # nocov` with an `Rcpp::stop`, per house style.
- **F3** — `R/adams.R` now validates a shared leaf set up front with a clear
  message (was an opaque `RenumberTips` failure); covered by a new test.
- **F4** — `@return` notes the `is.rooted()==FALSE` ape quirk on multifurcating
  roots; `Adams.Rd` regenerated.
- **F5** — version paper-trail: no code change (end state 9007 is correct).

## Verdict
**Do not block on correctness — block the *claim*.** The implementation is
correct (1000+ independent agreements; oracle clade-exact at n=137) and ~45x
faster than the R recursion. But its central advertised property is false as
written: the code is **O(kn log² n)**, not the documented **O(kn log n)**. That
one finding (F1) must be resolved before merge — either make the bound true
(radix-sort the two comparison sorts) or restate it honestly. Everything else is
minor/nit and does not gate.

## F1. Documented complexity O(kn log n) is actually O(kn log² n) — MAJOR (claim accuracy, not a correctness bug)
Where: `src/cons_adams.cpp:190` and `:402`; the claim appears in
`R/adams.R:31` (`@details \eqn{O(kn \log n)}`), `NEWS.md:8`, the
`src/cons_adams.cpp:10` header, and the commit message.

Two **comparison** sorts, each fresh at every recursion level, each adding an
independent log factor:
- `:190` — `buildNub` does `std::sort(nodes, by firstVisit)` over the ≤2m-1
  auxiliary nodes. O(m log m) per tree per call.
- `:402` — the meet-key grouping does `std::sort(idx, k-element comparator)`.
  O(k·sz·log sz) per spine step; Σ over steps ≤ O(k·m·log m) per call.

Grouped by recursion depth: side blocks are disjoint subsets of the leaves, so
Σ m ≤ n per depth level ⇒ Σ m log m ≤ n log n per level; depth is O(log n) (each
side block ≤ m/2) ⇒ **O(kn log² n)** total. Either sort alone gives the extra
log, so *both* must be removed to reach the paper's bound — fixing one is not
enough.

The paper (Jansson, Li & Sung 2017, Thm 2) achieves O(kn log n) precisely by
**radix-sorting** the meet vectors and building each T'_j once with incremental
leaf deletion (Lemma 5.2: T_j|B in linear time), never a comparison sort per
level. Both keys here are bounded integers (`firstVisit` < eulerLen ≤ 2·nNode;
nub node ids < 2m), so a 2-pass LSD radix / counting sort restores O(kn log n)
with no comparator. The `:345-351` counting sort by branchLevel is already
linear and correct — not part of the problem.

This was foreseen: the plan said "Acceptable simpler start: per-call std::sort →
O(kn log² n) … Be honest about which is shipped," and the scheduled
`adams-profile-optimize` job named both sort sites as its targets. That job
**fired** (2026-06-02 21:29Z, per `list_scheduled_tasks`) **but landed nothing** —
HEAD is still `53fe6be` from 19:37, the working tree is clean, and there is no new
`dev/profiling/` round log. So the optimisation that would make the claim true
never happened, and the overshoot shipped silently. **Decision required (the
user's):** make it true (radix) vs. restate the bound. Do not merge with the
claim as-is.

## F2. Defensive degree-1 / empty guards may be dead code, unmarked — MINOR (house style)
Where: `src/cons_adams.cpp:465-467` (`if (kids.empty()) continue;` and
`if (kids.size() == 1) deeper = kids[0];`), commented "the degree-1 guard is
defensive only."

The spine-bottom leaf always becomes `lastLeaf` (it has the max branchLevel and
is peeled last), so `deeper` is initialised to a real node before the chain
assembly in every case I traced; each step also yields ≥1 block, so `kids` is
never empty. If these branches are genuinely unreachable, house style (cf. the
Frequency / MajorityPlus commits that `// # nocov` their unreachable LCA/guard
arms) is to mark them so C++ coverage stays clean. If covr shows them *reached*,
then they are load-bearing and the "defensive only" comment is wrong. **Action:**
run covr on the new lines; mark `// # nocov` or add a reaching test accordingly.

## F3. Wrapper does not validate a shared leaf set — MINOR
Where: `R/adams.R:72,80-82`. `labels <- TipLabels(trees[[1]])` then
`RenumberTips(tr, labels)` on every tree. Adams requires Λ(T₁)=…=Λ(Tₖ); a tree
with a different leaf set yields a cryptic `RenumberTips` error rather than a
clear "all trees must share the same leaves" message. Sibling split methods route
through the same TreeTools machinery, so this is consistent with the package, but
a guarded `stop()` would be friendlier. Low priority.

## F4. `@return … rooted` vs `is.rooted() == FALSE` — NIT (doc)
`R/adams.R:43-44` says the result is "rooted by construction" (true: it has a
definite root and can carry no-input clades). But a multifurcating-root Adams
tree returns `is.rooted() == FALSE` under ape's strict 2-children heuristic. Not a
bug — matches the old behaviour and the maths — but a user testing `is.rooted()`
may be surprised. Optional one-line note.

## F5. Version paper-trail — NIT (resolved, no action)
Commit `7b42a9c` bumped DESCRIPTION to 0.0.0.9008; merge `53fe6be` reverted to
9007. End state (DESCRIPTION 9007, NEWS heading 9007) is correct and matches the
user's explicit "keep 9007." Flagged only so the git-history blip isn't mistaken
for a regression.

## Coverage / test-discrimination notes
- **v1 over-resolution is discriminated.** The at-scale oracle datasets produce
  non-trivial clade counts on incongruent inputs (indep n137/k10 → 27; perturb
  n50/k20 → 11; caterpillar distinct n40/k4 → 12) and match `fact.exe` rule 512
  exactly. v1's fixed-global-level bug inflated `mine` above `fact` on exactly
  these; the exact match (+ 680+ independent brute-force agreements incl. deep
  multi-level-compression cases) means a regression to v1 would be caught. No
  coverage gap here.
- The in-suite `refAdamsCladeSet()` reference uses `KeepTip` per node — a
  *different* derivation from the C++ nub/centroid machinery, so it is not
  tautological. The reviewer additionally built a fully independent brute force
  (`brute_adams.R`, ancestor/tipset signatures, no KeepTip) — see below.

## Verified fine (re-confirmed here)
- **Correctness, independently.** `brute_adams.R` derives "which child of the
  LCA" from ancestor/tipset membership on the original trees — no shared bug mode
  with the C++. Agreement: 320/320 binary, 350/350 multifurcating (root degree up
  to 15 / full star), 10/10 targeted adversarial (`fuzz3.R`, reproduced here:
  dup+conflict, spine-bottom-vs-early-branch, reversed caterpillars, star-vs-
  resolved, nested polytomies, balanced 8-leaf 3-way, idempotent polytomy,
  refinement pair, 4-way meet). On a provably-unique tree this is strong.
- **Oracle clade-exact** (`check-oracle.R`, exit 0, re-run here): Adams vs FACT
  rule 512 rooted=1 MATCH on n9/n10/n8 and at scale n50(×2)/n137/n40-caterpillars.
- **`heavyChild[j][curPos[j]]` is in-bounds.** `curPos[j]` is the *min*
  branchLevel over the remaining block; the `heavyChild` read fires only when
  `branchLevel[j][t] != curPos[j]` (i.e. strictly deeper), so the index is
  ≤ len-1, the filled range — never -1, never OOB.
- **`curPos[j]` never stuck at INT_MAX in the loop.** `removed` is global across
  trees, so `remaining ≥ 2` ⇒ every tree's pointer finds a non-removed leaf
  (`p < m`). The INT_MAX init is dead-defensive.
- **Scratch reuse is reentrancy-safe.** Each call consumes
  `branchLevel/sideChild/blockOf/removed/stamp` into local
  `stepBlocks`/`blockLeaves`/`blockNode` *before* any recursive call; nothing
  reads them post-recursion. `std::vector<char>` (not `vector<bool>`) is the
  correct choice. No file-scope mutable state; recursion depth O(log n).
- **Nested chain ≡ flat meet.** Leaves sharing a root-LCA child all carry the
  same branchLevel and are peeled in one step → assembled as one flat node
  (star → `(1,2,3,4)`, no spurious degree-1 node); nesting only across genuinely
  distinct spine levels (identical caterpillar → `(((a,b),c),d)`). Verified by
  trace and by the star/caterpillar oracle cases.
- **Rooting.** Each input marshalled on its OWN root (`Preorder(RenumberTips(...))`,
  not `.FactEdges`/taxon-1); output not re-rooted. Correct for a rooted method.
- **NAMESPACE / RcppExports.** No bare `KeepTip` left in `R/`;
  `adamsConsensusCpp(Rcpp::List, int)` registered with 2 args, signature matches
  the pre-registered stub exactly.
