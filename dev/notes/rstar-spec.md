# R* Consensus Tree — Implementation Specification

**Status:** Research spec — do NOT write package code from this file alone; resolve the OPEN QUESTIONS first.

**Package target:** `Consensus` (C:/Users/pjjg18/GitHub/Consensus)

**Primary sources read directly:**

- Degnan, J.H., DeGiorgio, M., Bryant, D., Rosenberg, N.A. (2009). "Properties of consensus methods for inferring species trees from gene trees." *Systematic Biology* 58(1):35–54. doi:10.1093/sysbio/syp008. [arXiv preprint 0802.2355 read in full.]
- Bansal, M.S., Dong, J., Fernández-Baca, D. (2009 preprint). "Comparing and Aggregating Partially Resolved Trees." arXiv:0906.5089. [For formal triplet/fan definitions.]
- Bryant, D., Berry, V. (2001). "A structured family of clustering and tree construction methods." *Advances in Applied Mathematics* 27:705–732. [The tree-construction algorithm referenced by Degnan et al.; not read directly — see Open Questions.]
- Bryant, D. (2003). "A classification of consensus methods for phylogenetics." In *Bioconsensus*, DIMACS Series 61:163–184. [Referenced throughout; not publicly accessible — see Open Questions.]
- Jansson, J., Sung, W.-K., Vu, H., Yiu, S.-M. (2016). "Faster algorithms for computing the R* consensus tree." *Algorithmica* 76:1224–1244. [Not read directly; complexity figures from secondary sources.]

---

## 1. Precise Definition of R* Consensus

### 1.1 Terminology

A **rooted triplet** (also called a rooted triple or 3-taxon statement) on leaf set {a, b, c} is a rooted binary phylogeny with exactly three leaves. Using the notation from Degnan et al. (2009, p. 37, arXiv p. 3):

> "We use the notation (AB)C for the three-taxon statement (rooted triple) that the most recent common ancestor (MRCA) of gene lineages A and B is not an ancestor of lineage C."

For any three taxa {a, b, c} in a rooted tree there are **exactly three possible resolved rooted triplet states**:
- ab|c  — MRCA(a, b) is a proper descendant of MRCA(a, c) = MRCA(b, c)
- ac|b  — MRCA(a, c) is a proper descendant of MRCA(a, b) = MRCA(b, c)
- bc|a  — MRCA(b, c) is a proper descendant of MRCA(a, b) = MRCA(a, c)

A **fan** (unresolved triplet / polytomy) is the 4th possible state: MRCA(a, b) = MRCA(a, c) = MRCA(b, c) (a single internal node with all three leaves as children). A fully binary (bifurcating) input tree induces only resolved triplets; partially resolved input trees can produce fan states.

**Critical note on fan in R* tallying:** Degnan et al. (2009) define "uniquely favored" by comparing only among the **three resolved** alternatives. The fan is **not listed as a competing state** in the selection rule. However, the paper's setting assumes fully binary input gene trees, so fans never arise in the tally in their presentation. See Section 5 (Open Questions) for whether fans should be treated as a 4th competing alternative when input trees are non-binary.

### 1.2 The Uniquely-Favored Rule — Exact Definition

The following is quoted or closely paraphrased from Degnan et al. (2009), p. 36 (journal), p. 2 (arXiv):

> "A rooted triple (AB)C on 3 taxa is said to be **uniquely favored** if it appears in more trees than either of the other 2 rooted triples, (AC)B or (BC)A, on the same set of 3 taxa."

**This is a strict plurality rule, not a strict majority rule:**

- ab|c is uniquely favored if and only if count(ab|c) > count(ac|b) AND count(ab|c) > count(bc|a).
- The threshold is NOT count > k/2; a triplet can be uniquely favored with far fewer than half the votes, as long as it beats each other alternative separately.
- In the case of any tie (count(ab|c) = count(ac|b), or count(ab|c) = count(bc|a), or a three-way tie), **no triplet is uniquely favored** for that set of taxa.

**Degnan et al. (2009, p. 36, journal) state explicitly:**

> "There might not be a uniquely favored triple for some sets of taxa."

**Consequence of ties:**

> "...if the set of rooted triples is incompatible or if there is a tie for the most frequently occurring rooted triple, the R* tree is declared unresolved or partially unresolved for those taxa causing the incompatibility." (Degnan et al. 2009, p. 36, arXiv p. 2.)

### 1.3 Step-by-Step R* Construction (Two Steps)

**Step 1 — Select uniquely-favored triplets.**  
For each of the C(n,3) = n!/(3!(n-3)!) 3-element subsets {a,b,c} of the n taxa:
1. Count how many of the k input trees have each resolved state: n_ab|c, n_ac|b, n_bc|a.  
   If input trees may be non-binary, also count n_fan (see Open Questions).
2. Determine the maximum: m = max(n_ab|c, n_ac|b, n_bc|a).
3. If exactly one state achieves m, that state is uniquely favored and is added to R*.
4. Otherwise (tie or fan-only), no triplet from {a,b,c} enters R*.

**Step 2 — Construct the most-resolved tree consistent with R*.**  
Build the most resolved rooted tree that contains (displays) exactly the triplets in R*, using the algorithm of Bryant and Berry (2001) (a generalization of Aho et al. BUILD). Degnan et al. (2009, p. 36, journal) state this constructs the tree "for example, using the algorithm of Bryant and Berry (2001), Corollary 2.2."

The resulting tree is the **Rooted Ancestral Consensus Tree (RACT)** (Degnan et al. 2009, p. 38, arXiv p. 4).

### 1.4 Clade Membership Rules (Explicit)

Degnan et al. (2009, p. 36, journal) give explicit rules for the 4-taxon case that generalize to n taxa:

> "1. Clades of sizes 1 and 4 [i.e., 1 and n] are included automatically.
> 2. The set {XY} is a clade exactly when (XY)Z and (XY)W are uniquely favored.
> 3. The set {XYZ} is a clade exactly when (XY)W, (XZ)W, and (YZ)W are uniquely favored."

For general n, with leaf set S, subset A ⊆ S (|A| ≥ 2, |A| ≤ n-1) is a clade in the R* tree if and only if for each pair A_i, A_j ∈ A (i ≠ j), and every Z ∈ S \ A, the triplet (A_i A_j)Z is uniquely favored. (Degnan et al. 2009, p. 36, journal, rules 1' and 2'.)

---

## 2. Construction Algorithm (BUILD)

### 2.1 The Aho et al. (1981) BUILD Algorithm

The BUILD algorithm (Aho, Sagiv, Szymanski, Ullman, 1981; also Bryant and Berry 2001, Corollary 2.2) takes as input a set R of rooted triplets on leaf set S and either constructs a rooted tree T that *displays* every triplet in R (i.e., is consistent with R), or reports that R is inconsistent (no such tree exists).

**The algorithm (recursive):**

```
BUILD(S, R):
  If |S| = 1: return the single-leaf tree.
  Construct the Aho graph G = (S, E) where edge {a,b} ∈ E iff
    there exists ab|c ∈ R  (a and b appear as the "close pair").
  Find connected components C_1, C_2, ..., C_m of G.
  If m = 1: FAIL — R is inconsistent (no rooted tree can display all triplets).
  For each component C_i:
    R_i = restriction of R to leaf sets within C_i
    T_i = BUILD(C_i, R_i)   [recurse]
    If T_i = FAIL: return FAIL.
  Create a new root node r; connect r to the root of each T_i.
  Return the resulting tree.
```

The tree returned by BUILD is consistent with R. When R is a complete and consistent set of triplets (i.e., the triplets uniquely specify the topology — which holds whenever the R* triplet set is fully informative), BUILD returns a unique fully-resolved tree. When R is not fully informative, BUILD may return a partially resolved tree. Bryant and Berry (2001), Corollary 2.2, is cited by Degnan et al. as the "most resolved" construction — the precise sense in which this extends BUILD is addressed in Open Question OQ3.

**Key property:** A set R of rooted triplets is consistent if and only if BUILD does not return FAIL (Steel 1992; Aho et al. 1981).

**Time complexity:**
- Naive BUILD: O(|S|^3) or equivalently O(n^3) for n = |S| taxa (because the Aho graph has O(n^2) edges and the algorithm makes O(n) recursive calls).
- With faster graph algorithms (Henzinger et al. 1999): O(|R| + n^2) for a single BUILD call.
- For R* with k input trees: Step 1 requires O(k * C(n,3)) = O(kn^3) time to tally all triplets. Step 2 (BUILD) runs in O(n^3). Total naive: **O(kn^3)**.
- Faster algorithms (Jansson, Sung, Vu, Yiu 2016): O(n^2 sqrt(log n)) for k=2; O(n^2 log^(4/3) n) for k=3; O(kn^2 log^(k+2) n) for unbounded k.

### 2.2 Consistency Guarantee and Partial Resolution Fallback

**Is R* always consistent?** Degnan et al. (2009, p. 36–37, journal; pp. 2, 12 arXiv) state explicitly:

> "We use the convention that if the set of rooted triples is incompatible or if there is a tie for the most frequently occurring rooted triple, the R* tree is declared unresolved or partially unresolved for those taxa causing the incompatibility."

This means R* is **not** guaranteed to be consistent. When BUILD returns FAIL on the full R*, the algorithm must **localise the incompatibility** and produce a partially-resolved tree. The mechanism Degnan et al. use is: when BUILD detects an inconsistency (single connected component), those taxa are left unresolved (their common ancestor in the R* tree becomes a multifurcation). This is analogous to how the majority-rule tree handles incompatible clades by omitting them.

**Worked example of incompatibility (Degnan et al. 2009, arXiv p. 12):**

With four input gene trees: (((AB)C)D), (((AD)C)B), (((BC)A)D), ((((CD)A)B):

- (AD)B and (AB)D each appear twice — so each is uniquely favored for {A,B,D}.
- But (AB)C and (BC)D are also each uniquely favored (appear in 2 of 3 input trees).
- These triplets (AB)C and (BC)D imply any tree containing both must also contain (AC)D, but the input trees have (AD)C as uniquely favored (Ranwez et al. 2007 logic). The triplet set is inconsistent.
- Result: the R* tree is partially unresolved — specifically for the taxa causing the incompatibility.

**Practical fallback for implementation:** When BUILD returns FAIL on a subset of taxa during recursion, collapse those taxa to a multifurcation at their lowest common ancestor in the partial tree. This ensures BUILD always terminates with a (possibly star-like) tree rather than an error.

### 2.3 Rooting Assumption

R* is defined for **rooted input trees**. The method requires rooted trees because the state of a rooted triplet (which taxon is the outgroup) depends on the root. The "most probable rooted triple in the gene tree distribution" is the one matching the species tree's rooting (Degnan et al. 2009, Lemma 2, p. 38 arXiv/p. 38 journal).

If input trees are unrooted, an outgroup must be specified to root them before applying R*. There is no published version of R* for unrooted input (see Open Questions).

---

## 3. Worked Example

**Setup:** 4 taxa (A, B, C, D), 5 input trees, species tree = (((AB)C)D).

**Input trees:**

| Tree | Topology            |
|------|---------------------|
| T1   | (((AB)C)D)          |
| T2   | (((AB)C)D)          |
| T3   | (((AB)C)D)          |
| T4   | (((AC)B)D)          |
| T5   | (((BC)A)D)          |

k = 5 input trees.

**Per-triplet tally** (C(4,3) = 4 triplets):

**{A, B, C}:**
- T1,T2,T3 each have (AB)C (A and B closer together, C outgroup)
- T4 has (AC)B
- T5 has (BC)A
- Counts: (AB)C = 3, (AC)B = 1, (BC)A = 1
- Maximum = 3; unique winner = (AB)C. **Uniquely favored: (AB)C.**

**{A, B, D}:**
- All five trees: in T1,T2,T3 the tree is (((AB)C)D) — reading {A,B,D}: A and B are sister, so (AB)D.
- T4: (((AC)B)D) — reading {A,B,D}: the MRCA(A,B) is deeper than MRCA(A,D)? No — the structure is (((AC)B)D), so MRCA(A,C) < MRCA(A,B) < MRCA(A,D). For {A,B,D}: MRCA(A,B) < MRCA(A,D) = MRCA(B,D), so triplet is (AB)D.
- T5: (((BC)A)D) — for {A,B,D}: MRCA(B,C) < MRCA(B,A) < MRCA(B,D). For {A,B,D}: MRCA(A,B) < MRCA(A,D) = MRCA(B,D), so triplet is (AB)D.
- Counts: (AB)D = 5, (AD)B = 0, (BD)A = 0.
- **Uniquely favored: (AB)D.**

**{A, C, D}:**
- T1,T2,T3: (((AB)C)D) — for {A,C,D}: MRCA(A,B) < MRCA(A,C) < MRCA(C,D). So MRCA(A,C) < MRCA(A,D) = MRCA(C,D). Triplet: (AC)D.
- T4: (((AC)B)D) — for {A,C,D}: MRCA(A,C) < MRCA(A,D) = MRCA(C,D). Triplet: (AC)D.
- T5: (((BC)A)D) — for {A,C,D}: MRCA(B,C) < MRCA(A,C) (since A joins the BC clade), then MRCA(A,BC) < MRCA(A,D)... Rereading: the tree is ((BC)(A))D overall, so MRCA(A,C) = MRCA(A,B) at the node joining A to (BC). For {A,C,D}: MRCA(A,C) < MRCA(C,D). Triplet: (AC)D.
- Counts: (AC)D = 5, (AD)C = 0, (CD)A = 0.
- **Uniquely favored: (AC)D.**

**{B, C, D}:**
- T1,T2,T3: For {B,C,D}: MRCA(B,C) at node (AB)C < MRCA(B,D). Triplet: (BC)D.
- T4: (((AC)B)D): For {B,C,D}: MRCA(A,C) < MRCA(B,(AC)) < MRCA(B,D)... The deepest node grouping B with C is the node joining B to (AC): MRCA(B,C) = that node. MRCA(B,D) is root. MRCA(C,D) = MRCA(B,D) = root. So MRCA(B,C) < MRCA(B,D) = MRCA(C,D). Triplet: (BC)D.
- T5: (((BC)A)D): MRCA(B,C) < MRCA(B,D) = MRCA(C,D). Triplet: (BC)D.
- Counts: (BC)D = 5, (BD)C = 0, (CD)B = 0.
- **Uniquely favored: (BC)D.**

**R* triplet set:** { (AB)C, (AB)D, (AC)D, (BC)D }

**Consistency check via Aho graph:**

BUILD on {A,B,C,D} with R* = {(AB)C, (AB)D, (AC)D, (BC)D}:
- Aho graph edges (pair appears as "close" pair): A-B (from (AB)C, (AB)D), A-C (from (AC)D), B-C (from (BC)D)
- Aho graph: A-B-C all connected, D isolated.
- Components: {A,B,C} and {D}.
- Root joins T({A,B,C}) and T({D}).
- Recurse on {A,B,C} with R restricted = {(AB)C}:
  - Aho graph: A-B edge (from (AB)C); C isolated.
  - Components: {A,B} and {C}. Root joins ({A,B}) and ({C}).
  - Recurse on {A,B}: single edge, return cherry (A,B).
  - Recurse on {C}: return leaf C.
  - T({A,B,C}) = ((AB)C).
- T({D}) = leaf D.
- **R* tree = (((AB)C)D)** ✓ — matches the species tree.

**Unit-test assertions:**
```
trees <- list(
  read.tree(text = "(((A,B),C),D);"),  # T1
  read.tree(text = "(((A,B),C),D);"),  # T2
  read.tree(text = "(((A,B),C),D);"),  # T3
  read.tree(text = "(((A,C),B),D);"),  # T4
  read.tree(text = "(((B,C),A),D);")   # T5
)
result <- RStar(trees)
# Expected: (((A,B),C),D)
# Specifically:
#   {A,B,C} → (AB)C (3 vs 1 vs 1)   — plurality, not majority
#   {A,B,D} → (AB)D (5 vs 0 vs 0)
#   {A,C,D} → (AC)D (5 vs 0 vs 0)
#   {B,C,D} → (BC)D (5 vs 0 vs 0)
```

**Tie-producing example** — to test that ties produce unresolved output, use k=4 trees where {A,B,C} splits 2-2-0:

```
# T1,T2 = (((A,B),C),D)  → {A,B,C}: (AB)C=2
# T3,T4 = (((A,C),B),D)  → {A,B,C}: (AC)B=2
# For {A,B,C}: tie 2-2-0 → no uniquely-favored triplet → {A,B,C} node unresolved
# Expected R* tree: ((A,B,C),D) — the {A,B,C} clade exists but is an unresolved polytomy
```

**Plurality-vs-majority discriminating example** — k=4 trees where the plurality winner has only 2 out of 4 votes (below strict majority):

```
# T1,T2 = (((A,B),C),D) → {A,B,C}: (AB)C
# T3    = (((A,C),B),D) → {A,B,C}: (AC)B
# T4    = (((B,C),A),D) → {A,B,C}: (BC)A
# Counts: (AB)C=2, (AC)B=1, (BC)A=1 → max=2, unique → (AB)C uniquely favored
# Under strict majority (>k/2=2): (AB)C=2 is NOT strictly greater than 2 → rejected
# Under plurality: (AB)C=2 > max(1,1) → accepted
# A correct R* implementation must select (AB)C here.
```

---

## 4. Implementation Notes for R + TreeTools

### 4.1 Enumerating Rooted Triplets and Their States

**Approach:** For each input tree T (a `phylo` object), iterate over all C(n,3) 3-subsets {a,b,c}. For each subset, determine the rooted-triplet state by comparing MRCA depths.

**Key observation:** For a rooted tree, the triplet state of {a,b,c} is determined by which pair has the deepest (most recent) MRCA. Given the MRCA matrix M (where M[i,j] is the node index of the MRCA of tips i and j, in a preorder-numbered tree where ancestors have smaller node numbers than their descendants), the identity rule is:

For triplet {i, j, l}, compare M[i,j], M[i,l], M[j,l]:
- Exactly two of these three values will be equal (the triplet's root — the LCA of all three), and exactly one will be unique and strictly larger (deeper).
- If M[i,j] is the unique maximum (M[i,j] > M[i,l] = M[j,l]): state is ij|l.
- If M[i,l] is the unique maximum (M[i,l] > M[i,j] = M[j,l]): state is il|j.
- If M[j,l] is the unique maximum (M[j,l] > M[i,j] = M[i,l]): state is jl|i.
- If all three values are equal (all the same node): the triplet is a fan for this tree.

This identity rule works because, in any rooted tree, for any three leaves there is always one pair whose MRCA is a proper descendant of the other pairs' MRCA — unless the tree is unresolved at their common ancestor (polytomy), in which case all three MRCAs coincide. The rule does **not** require converting node indices to depths; in preorder numbering, within any triplet, the deeper MRCA always has a strictly larger node index than the shallower one, so the comparison M[i,j] vs M[i,l] vs M[j,l] correctly identifies the close pair.

**TreeTools functions to use:**
- `TreeTools::Preorder()` — put tree in preorder before ancestry queries (required so that ancestor node numbers are strictly smaller than descendant node numbers).
- `TreeTools::ListAncestors()` (internal via `AllAncestors()`) — get the full ancestor list for each node.
- `TreeTools::MRCA(x1, x2, ancestors)` — compute the MRCA node of two tips, returning its node index. In a preorder-numbered tree, MRCA node index increases monotonically as you go deeper: MRCA(a,b) < MRCA(a,b's child, ...) because the ancestor has a smaller preorder index than its descendants.
- `TreeTools::DescendantTips()` or the C-level ancestry arrays — for fast batch MRCA computation across all triplets.
- **Quartet package analogy:** The `Quartet::QuartetStates()` function provides a model for how to enumerate all quartet states efficiently. A `RootedTripletStates()` function would follow the same pattern but over C(n,3) subsets and using MRCA node-identity comparison.

**Efficient batch computation:** Precompute an n × n MRCA matrix (using `MRCA()` for each pair, with a shared `ancestors` list built once via `ListAncestors()`) in O(n^2 * depth) per tree (O(n^3) worst case for balanced trees), then tally triplet states from this matrix in O(n^3) total. This is the standard approach for an R-level implementation.

**R-level pseudo-code:**

```r
# For one tree: build MRCA matrix using node indices (preorder numbering).
# In preorder, MRCA(a, b) has a SMALLER node index than either a or b,
# and a LARGER node index than any ancestor of MRCA(a,b).
# Thus comparing M[i,j], M[i,l], M[j,l] by node index correctly identifies
# the close pair (the unique maximum is the deepest MRCA = the close pair's MRCA).
mrca_matrix <- function(tree) {
  n <- length(tree$tip.label)
  tree <- Preorder(tree)
  anc <- ListAncestors(tree$edge[, 1], tree$edge[, 2])
  mat <- matrix(0L, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- mat[j, i] <- MRCA(i, j, anc)
    }
  }
  mat
}

# Triplet state from MRCA matrix:
triplet_state <- function(M, i, j, l) {
  mij <- M[i, j]; mil <- M[i, l]; mjl <- M[j, l]
  if (mij == mil && mij == mjl) return("fan")  # all same → polytomy
  if (mij > mil && mij > mjl)   return("ij|l")
  if (mil > mij && mil > mjl)   return("il|j")
  if (mjl > mij && mjl > mil)   return("jl|i")
  stop("unexpected tie in MRCA values")  # should not occur in a bifurcating tree
}
```

**Performance note:** For large n, the O(n^3) loop over triplets should be implemented in C++ (as in the Quartet package). For the initial implementation, an R loop over C(n,3) subsets with precomputed MRCA matrix is acceptable.

### 4.2 Tally Array

Maintain a 3D tally array of dimension C(n,3) × 3 (one row per 3-subset, three columns for the three states), accumulating counts across k input trees. Optionally a 4th column for fan counts if needed (see Open Questions).

Indexing: for taxa ordered 1 … n, a canonical index for {i, j, l} (i < j < l) can be used to address the tally array. The canonical assignment of states is:
- state 1 = ij|l  (i and j are the close pair)
- state 2 = il|j  (i and l are the close pair)
- state 3 = jl|i  (j and l are the close pair)

### 4.3 BUILD Implementation

The BUILD algorithm operates on a set of rooted triplets. In R:

```r
build <- function(taxa, triplets) {
  # triplets: data.frame or matrix with columns close1, close2, out
  # (taxa a,b closer than c → row c(a, b, c))
  if (length(taxa) <= 2) return(taxa)  # leaf or cherry

  # Construct Aho graph: connect i--j if any triplet has {i,j} as close pair
  # Find connected components via union-find or igraph
  # If only 1 component: R is inconsistent for this subset → collapse to polytomy
  # Otherwise recurse on each component
  ...
}
```

The R implementation can use `igraph::components()` for the connected-components step, keeping it readable. The C++ implementation should use union-find for speed.

### 4.4 Handling Non-Binary Input Trees

If any input tree has polytomies (multifurcations), some triplets will be fans (MRCA depths are tied). Implementors must decide (see Open Questions) whether fans count as a 4th competing state or are treated as abstentions.

**Recommended conservative approach (pending author clarification):** Treat fans as abstentions — do not count them in any state's tally. A triplet in a fan node of an input tree contributes 0 to all three resolved states. This is consistent with the spirit of Degnan et al.'s formulation, where "more trees than either of the other 2 rooted triples" is compared purely among resolved states.

### 4.5 R* Requires Rooted Input

All input trees must be rooted. The `phylo` root is the first edge's parent. Trees should be checked with `ape::is.rooted()`. If unrooted trees are provided, the user must supply an outgroup taxon for rooting; the package should not silently root by the first tip.

### 4.6 Edge Cases

| Situation | Expected behaviour |
|-----------|-------------------|
| k = 1 tree | R* returns that tree itself (every triplet is uniquely favored 1-0-0) |
| k = 2 trees, same topology | R* returns that topology |
| k = 2 trees, different topologies | All triplets that differ are tied 1-1-0; only shared triplets enter R*; result is usually a star or partially unresolved tree |
| All input trees are star trees | All triplets are fans; R* = star tree |
| n = 3 taxa | One triplet; if uniquely favored, fully resolved 3-taxon tree; otherwise star |
| Incompatible R* triplet set | BUILD fails locally → collapse those taxa to polytomy |
| Identical taxa labels required | All trees must share exactly the same leaf set (unlike supertree methods) |

---

## 5. Open Questions

> **RESOLVED — see `rstar-findings.md`.** The user supplied the primary source
> `Jansson2016a` (Jansson, Sung, Vu & Yiu 2016), whose §1.1 + Lemma 1.1 settle the
> implementation-critical questions: fans have no impact (OQ1); the R* tree is the
> unique tree whose clusters are the **strong clusters** of `R_maj`, which always
> exists (so OQ2's collapse case does not arise and OQ3's assembly is fixed — it
> is the strong-cluster construction, not BUILD). `RStar()` implements this and is
> validated by a brute-force strong-cluster oracle. The questions below are
> retained for historical context.

These questions could not be definitively resolved from the literature and should be directed to the authors (Bryant, Degnan) before coding.

**OQ1 — Fan treatment in non-binary input [CRITICAL].**  
Degnan et al. (2009) define R* exclusively for fully-resolved (binary) input gene trees. When input trees have polytomies, a triplet {a,b,c} may be unresolved (fan state) in some input trees. It is not stated whether:
(a) Fans are excluded from the tally (abstentions), so "uniquely favored" is computed over the resolved-only subset of input trees; or  
(b) Fans are counted as a 4th state, and ab|c is uniquely favored only if it beats all three alternatives including the fan count.  
The wording "more trees than either of the other 2 rooted triples" suggests (a), but this is based on a binary-input context. **Implementation choice (a) is the safer default** until confirmed.

**OQ2 — Incompatibility resolution mechanism [CRITICAL].**  
Degnan et al. (2009) state that incompatible R* triplet sets produce a "partially unresolved" tree, but do not specify the exact algorithm for localising the incompatibility. Bryant and Berry (2001, Corollary 2.2) is cited as the construction algorithm. The precise mechanism — which taxa are collapsed to a polytomy when BUILD fails — needs to be confirmed. Possible mechanisms:
(a) Collapse all taxa in the failing connected component to a polytomy.  
(b) Iteratively remove the least-supported conflicting triplets until BUILD succeeds.  
(c) Apply BUILD recursively, collapsing each failed sub-problem to a polytomy.  
The paper's example (arXiv p. 12) implies (a) or (c), but needs confirmation.

**OQ3 — Bryant and Berry (2001) "most resolved" construction.**  
The R* construction calls for "the most resolved tree containing only uniquely-favored triplets" (Degnan et al. 2009, p. 36). Standard BUILD (Aho et al. 1981) returns the *least* resolved consistent tree. Bryant and Berry (2001) Corollary 2.2 is cited as giving the construction; the precise additional step (if any) to ensure maximal resolution from the consistent triplet set is not confirmed from secondary sources. **Access to Bryant and Berry (2001) is needed to verify this.** The paper is at: https://www.sciencedirect.com/science/article/pii/S0196885801907584.

**OQ4 — Behaviour when all triplets for a 3-subset are equally tied (three-way tie).**  
If n_ab|c = n_ac|b = n_bc|a (including the case where all three counts are 0, which occurs when all input trees have a fan for that 3-subset), what is the correct output? The paper says "unresolved" but does not distinguish between a pairwise tie and a three-way tie. Both presumably produce the same outcome (no triplet favored, those taxa left unresolved), but confirmation is useful.

**OQ5 — Consistency with the majority-rule tree — which direction.**  
Degnan et al. (2009, p. 36) state "the greedy and R* consensus trees are always resolutions of the majority-rule tree (Bryant 2003), meaning that every clade on the majority-rule consensus tree is also on the greedy and R* consensus trees." This is referenced to Bryant (2003) but requires access to that paper to confirm the mechanism. Implementation consequence: if the R* tree is always a refinement of majority-rule, one optimisation could build majority-rule first and only run R* within each polytomy — but this should be verified.

**OQ6 — Handling taxa subsets (missing taxa) in some input trees.**  
R* as defined assumes all input trees have the same leaf set. The Degnan et al. paper only considers the identical-leaf-set case. If some input trees are missing taxa (a common practical scenario), it is unclear how to tally triplets that include missing taxa. This would require a supertree-style generalisation not covered in the 2009 paper.

---

## 6. References

- Aho, A.V., Sagiv, Y., Szymanski, T.G., Ullman, J.D. (1981). "Inferring a tree from lowest common ancestors with an application to the optimization of relational expressions." *SIAM Journal on Computing* 10(3):405–421.
- Bansal, M.S., Dong, J., Fernández-Baca, D. (2010). "Comparing and aggregating partially resolved trees." *Theoretical Computer Science* 412(48):6634–6652. [arXiv:0906.5089]
- Bryant, D. (2003). "A classification of consensus methods for phylogenetics." In Janowitz, M.F. *et al.* (eds), *Bioconsensus*, DIMACS Series in Discrete Mathematics and Theoretical Computer Science 61:163–184. American Mathematical Society.
- Bryant, D., Berry, V. (2001). "A structured family of clustering and tree construction methods." *Advances in Applied Mathematics* 27:705–732. doi:10.1006/aama.2001.0758.
- Degnan, J.H., DeGiorgio, M., Bryant, D., Rosenberg, N.A. (2009). "Properties of consensus methods for inferring species trees from gene trees." *Systematic Biology* 58(1):35–54. doi:10.1093/sysbio/syp008.
- Henzinger, M.R., King, V., Warnow, T. (1999). "Constructing a tree from homeomorphic subtrees, with applications to computational evolutionary biology." *Algorithmica* 24:1–13.
- Jansson, J., Sung, W.-K., Vu, H., Yiu, S.-M. (2016). "Faster algorithms for computing the R* consensus tree." *Algorithmica* 76:1224–1244. doi:10.1007/s00453-016-0122-2.
- Jansson, J., Sung, W.-K. (2012). "Constructing the R* consensus tree of two trees in subcubic time." *Algorithmica* 67:329–351. doi:10.1007/s00453-012-9639-1.
- Steel, M. (1992). "The complexity of reconstructing trees from qualitative characters and subtrees." *Journal of Classification* 9:91–116.
