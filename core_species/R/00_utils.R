# Shared helpers for the core-species network analysis.
#
# The module reads existing network RDS files and writes tables under
# `output$results`; it never modifies the upstream network construction.

need_file <- function(path, label) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path) ||
      !file.exists(path)) {
    stop(label, " not found: ", path)
  }
  normalizePath(path)
}

read_species <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"species" %in% names(x)) stop("Missing 'species' column: ", path)
  unique(as.character(x$species))
}

safe_fraction <- function(n, d) ifelse(d == 0, NA_real_, n / d)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_core_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

bind_rows_safe <- function(x) {
  x <- purrr::compact(x)
  if (!length(x)) data.frame() else dplyr::bind_rows(x)
}

abundance_class <- function(x) {
  dplyr::case_when(
    x == "high" ~ "abundant",
    x %in% c("low", "very_low") ~ "rare",
    TRUE ~ "other"
  )
}

is_abundance_eligible <- function(x) {
  x %in% c("abundant", "rare")
}

# Resolve a path from the config against the directory that holds the corpus
# (the repository root by default), so recipes do not depend on the caller's
# working directory.
resolve_core_path <- function(path, root) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) {
    return(path)
  }
  if (grepl("^/", path)) path else file.path(root, path)
}

prepare_output_dirs <- function(outdir) {
  dirs <- file.path(outdir, c("qc", "tables", "plot_data"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}
