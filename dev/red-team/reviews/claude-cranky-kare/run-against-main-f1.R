# Run cranky-kare F1 (zero-length incompatible interior edge -> infinite loop)
# against CURRENT main source. If this script is killed by the OS timeout, the
# hang is LIVE on main. If it prints a finite distance, F1 is fixed.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
cat("=== F1: zero-length incompatible interior edge ===\n"); flush.console()
A <- ape::read.tree(text="(0:1,(1:0,2:0):0,3:1,4:1,5:1);")  # interior split {1,2}, length 0
B <- ape::read.tree(text="(0:1,(2:1,3:1):1,1:1,4:1,5:1);")  # {2,3}:1 incompatible with {1,2}
cat("calling BHVDistance (may hang) ...\n"); flush.console()
d <- BHVDistance(A, B)
cat("d =", d, "\n")
cat("=== F1 VERDICT: FIXED (returned", d, "without hanging) ===\n")
