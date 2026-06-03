source('C:/Users/pjjg18/GitHub/worktrees/Consensus/adams/dev/red-team/reviews/feature-adams-fast/brute_adams.R')

check <- function(label, trees) {
  trees <- structure(trees, class = "multiPhylo")
  cpp <- tryCatch(cladeSetLab(Adams(trees)),
                  error = function(e) paste("ERR:", conditionMessage(e)))
  bru <- bruteAdamsClades(trees)
  ok <- identical(cpp, bru)
  cat(sprintf("[%s] %s\n", if (ok) "OK " else "FAIL", label))
  if (!ok) {
    cat("   cpp  :", paste(cpp, collapse=" ; "), "\n")
    cat("   brute:", paste(bru, collapse=" ; "), "\n")
  }
  ok
}

rd <- function(x) ape::read.tree(text = x)

# ----- targeted adversarial cases -----

# (1) mixed duplicates + conflict: 3 copies of T1, 1 conflicting T2
T1 <- rd("(((a,b),(c,d)),(e,f));")
T2 <- rd("(((a,c),(b,e)),(d,f));")
check("3xT1 + 1xT2 dup+conflict", list(T1,T1,T1,T2))

# (2) leaf is spine-bottom in one tree, branches early in another
#   tree A: caterpillar a,b,c,d,e (e deepest)
#   tree B: e branches at the root
A <- rd("((((a,b),c),d),e);")
B <- rd("(e,(((a,b),c),d));")
check("spine-bottom vs early-branch leaf", list(A,B))

# (3) k=2 deep asymmetry, conflicting deep nodes
A2 <- rd("(a,(b,(c,(d,(e,(f,g))))));")
B2 <- rd("(g,(f,(e,(d,(c,(b,a))))));")
check("k=2 reversed caterpillars", list(A2,B2))

# (4) high-degree root in some, resolved in others
S <- rd("(a,b,c,d,e,f);")        # full star
R <- rd("(((a,b),(c,d)),(e,f));") # fully resolved
check("star vs resolved", list(S,R))
check("star vs resolved vs resolved2", list(S,R,rd("(((a,c),(b,d)),(e,f));")))

# (5) deeply nested polytomies that conflict on nesting depth
P1 <- rd("((a,b,c),(d,e,f));")
P2 <- rd("((a,b),(c,d,e,f));")
P3 <- rd("(a,(b,c),(d,e),f);")
check("nested polytomies x3", list(P1,P2,P3))

# (6) one tree where a whole side-block is itself the heavy child elsewhere
H1 <- rd("(((a,b),(c,d)),((e,f),(g,h)));")
H2 <- rd("(((a,b),(e,f)),((c,d),(g,h)));")
H3 <- rd("(((a,e),(b,f)),((c,g),(d,h)));")
check("balanced 8-leaf 3-way", list(H1,H2,H3))

# (7) repeated identical multifurcating tree (idempotence on polytomy)
M <- rd("((a,b,c),(d,e),(f,g,h));")
check("idempotent polytomy", list(M,M,M))

# (8) two trees, one a refinement of the other
C1 <- rd("(((a,b),c),(d,e));")
C2 <- rd("((a,b,c),(d,e));")
check("refinement pair", list(C1,C2))

# (9) leaves swap which child of root across trees in a complex meet
W1 <- rd("((a,b,c,d),(e,f,g,h));")
W2 <- rd("((a,b,e,f),(c,d,g,h));")
W3 <- rd("((a,c,e,g),(b,d,f,h));")
W4 <- rd("((a,d,f,g),(b,c,e,h));")
check("4-way meet 8 leaves", list(W1,W2,W3,W4))
