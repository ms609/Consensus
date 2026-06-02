library("ConsTree")

example_db <- tools::Rd_db("ConsTree")

cat("Running examples for", length(example_db), "topics\n")

failures <- character(0)

for (topic in names(example_db)) {
  cat("\n>>> Example:", topic, "\n")
  ex <- tools::Rd2ex(example_db[[topic]])

  if (length(ex) == 0L) {
    cat("No example found for topic:", topic, "\n")
    next
  }

  ex_file <- tempfile(fileext = ".R")
  writeLines(ex, ex_file)

  tryCatch(
    {
      sys.source(ex_file, envir = globalenv())
      cat("✓ Success:", topic, "\n")
    },
    error = function(e) {
      cat("✘ Error in topic:", topic, "\n", conditionMessage(e), "\n")
      failures <<- c(failures, topic)
    }
  )
}
cat("\nFinished running examples.\n")

if (length(failures)) {
  cat("❌ Failures in", length(failures), "topics:\n")
  print(failures)
  quit(status = 1)
} else {
  cat("✅ All examples ran successfully.\n")
}
