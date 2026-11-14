#!/usr/bin/env Rscript
# Compare two core-species result trees.
#
# Usage:
#   Rscript core_species/verification/compare_core_results.R REFERENCE CANDIDATE
#
# Tables are compared by shape, column names and values.  Numeric columns use a
# tolerance of 1e-9, which absorbs last-digit differences in unseeded
# eigenvector centrality.  run_manifest.tsv and run_summary.tsv record paths
# and timestamps, so they are reported as skipped unless --strict is given.

args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args
args <- args[args != "--strict"]
if (length(args) != 2) {
  stop("Usage: compare_core_results.R REFERENCE CANDIDATE [--strict]")
}
reference <- args[1]
candidate <- args[2]

volatile <- c("run_manifest.tsv", "run_summary.tsv")

files_in <- function(root) {
  list.files(root, recursive = TRUE, full.names = FALSE)
}

ref_files <- files_in(reference)
cand_files <- files_in(candidate)
missing <- setdiff(ref_files, cand_files)
extra <- setdiff(cand_files, ref_files)

problems <- character()
skipped <- character()

compare_table <- function(path, rel) {
  reader <- if (grepl("\\.csv$", rel)) read.csv else read.delim
  a <- reader(path, stringsAsFactors = FALSE, check.names = FALSE)
  b <- reader(file.path(candidate, rel),
    stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(dim(a), dim(b))) {
    problems <<- c(problems, paste0(rel, ": shape ", paste(dim(a), collapse = "x"),
      " vs ", paste(dim(b), collapse = "x")))
    return(invisible(NULL))
  }
  if (!identical(names(a), names(b))) {
    problems <<- c(problems, paste0(rel, ": column names differ"))
    return(invisible(NULL))
  }
  for (col in names(a)) {
    x <- a[[col]]
    y <- b[[col]]
    if (is.numeric(x) && is.numeric(y)) {
      ok <- isTRUE(all.equal(x, y, tolerance = 1e-9, check.attributes = FALSE))
    } else {
      ok <- identical(as.character(x), as.character(y))
    }
    if (!ok) problems <<- c(problems, paste0(rel, ": column ", col, " differs"))
  }
}

for (rel in intersect(ref_files, cand_files)) {
  if (grepl("\\.(tsv|csv)$", rel)) {
    if (basename(rel) %in% volatile && !strict) {
      skipped <- c(skipped, rel)
      next
    }
    compare_table(file.path(reference, rel), rel)
  }
}

cat("Reference:", normalizePath(reference), "\n")
cat("Candidate:", normalizePath(candidate), "\n")
cat("Common files:", length(intersect(ref_files, cand_files)), "\n")
if (length(missing)) cat("Missing in candidate:\n  ",
  paste(missing, collapse = "\n  "), "\n", sep = "")
if (length(extra)) cat("Extra in candidate:\n  ",
  paste(extra, collapse = "\n  "), "\n", sep = "")
if (length(skipped)) cat("Skipped volatile files:\n  ",
  paste(skipped, collapse = "\n  "), "\n", sep = "")

if (length(problems)) {
  cat("Differences:\n  ", paste(problems, collapse = "\n  "), "\n", sep = "")
  quit(status = 1)
}
if (length(missing) || length(extra)) quit(status = 1)
cat("OK: all compared tables match.\n")
