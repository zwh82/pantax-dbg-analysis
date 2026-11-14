#!/usr/bin/env Rscript
# Compare two analysis result trees produced by the legacy scripts and by
# analysis/run_analysis.R.
#
#   Rscript analysis/verification/compare_analysis.R \
#     <old_results> <new_results> <old_plots> <new_plots>
#
# Tables, RDS objects and text files are compared by value.  PNG/PDF outputs
# are skipped; the plot objects captured by verification/prelude.R are rendered
# with a fixed seed and compared pixel-for-pixel instead.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: compare_analysis.R <old_results> <new_results> ",
       "<old_plots> <new_plots>")
}
old_dir <- normalizePath(args[1])
new_dir <- normalizePath(args[2])
old_plots <- normalizePath(args[3])
new_plots <- normalizePath(args[4])

suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
})

list_files <- function(dir) {
  files <- list.files(dir, recursive = TRUE, full.names = FALSE, all.files = FALSE)
  sort(files[!dir.exists(file.path(dir, files))])
}

issues <- character()
checks <- 0

add_issue <- function(msg) {
  issues <<- c(issues, msg)
}

compare_vector <- function(a, b, path) {
  if (is.factor(a) || is.factor(b)) {
    if (!is.factor(a) || !is.factor(b)) {
      add_issue(paste0(path, ": factor vs non-factor"))
      return(invisible())
    }
    if (!identical(levels(a), levels(b))) {
      add_issue(paste0(path, ": factor levels differ"))
    }
    if (!identical(is.ordered(a), is.ordered(b))) {
      add_issue(paste0(path, ": ordered flag differs"))
    }
    a <- as.character(a)
    b <- as.character(b)
  }
  if (is.numeric(a) || is.numeric(b)) {
    if (!isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = TRUE))) {
      add_issue(paste0(path, ": numeric values differ"))
    }
  } else if (!identical(a, b)) {
    add_issue(paste0(path, ": values differ"))
  }
}

compare_object <- function(a, b, path) {
  if (is.null(a) || is.null(b)) {
    if (!is.null(a) || !is.null(b)) add_issue(paste0(path, ": NULL mismatch"))
    return(invisible())
  }
  if (is.data.frame(a) || is.data.frame(b)) {
    if (!is.data.frame(a) || !is.data.frame(b)) {
      add_issue(paste0(path, ": data.frame vs other type"))
      return(invisible())
    }
    if (!identical(names(a), names(b))) {
      add_issue(paste0(path, ": column names differ"))
    }
    if (nrow(a) != nrow(b)) {
      add_issue(paste0(path, ": row counts differ (", nrow(a), " vs ", nrow(b), ")"))
    }
    for (nm in intersect(names(a), names(b))) {
      compare_vector(a[[nm]], b[[nm]], paste0(path, "$", nm))
    }
    return(invisible())
  }
  if (is.list(a) && is.list(b)) {
    if (!identical(names(a), names(b))) {
      add_issue(paste0(path, ": list names differ"))
    }
    for (nm in intersect(names(a), names(b))) {
      compare_object(a[[nm]], b[[nm]], paste0(path, "$", nm))
    }
    return(invisible())
  }
  if (is.environment(a) || is.environment(b)) return(invisible())
  compare_vector(a, b, path)
}

compare_file <- function(rel) {
  old_file <- file.path(old_dir, rel)
  new_file <- file.path(new_dir, rel)
  checks <<- checks + 1

  if (!file.exists(new_file)) {
    add_issue(paste0(rel, ": missing in new tree"))
    return(invisible())
  }
  if (identical(unname(tools::md5sum(old_file)), unname(tools::md5sum(new_file)))) {
    return(invisible())
  }

  ext <- tolower(tools::file_ext(rel))
  if (ext %in% c("png", "pdf")) return(invisible())

  if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    a <- read.table(old_file, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    b <- read.table(new_file, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    compare_object(a, b, rel)
  } else if (ext == "rds") {
    a <- readRDS(old_file)
    b <- readRDS(new_file)
    compare_object(a, b, rel)
  } else if (ext == "txt") {
    if (!identical(readLines(old_file, warn = FALSE),
                   readLines(new_file, warn = FALSE))) {
      add_issue(paste0(rel, ": text differs"))
    }
  } else {
    add_issue(paste0(rel, ": binary content differs"))
  }
}

render_plot <- function(plot, path) {
  set.seed(1)
  ggplot2::ggsave(
    filename = path, plot = plot, width = 6, height = 6,
    units = "in", dpi = 100, limitsize = FALSE, bg = "white"
  )
}

compare_plot <- function(key) {
  old_file <- file.path(old_plots, key)
  new_file <- file.path(new_plots, key)
  checks <<- checks + 1

  if (!file.exists(new_file)) {
    add_issue(paste0("plot ", key, ": missing in new run"))
    return(invisible())
  }

  tmp_old <- tempfile(fileext = ".png")
  tmp_new <- tempfile(fileext = ".png")
  on.exit(unlink(c(tmp_old, tmp_new)), add = TRUE)

  result <- tryCatch({
    render_plot(readRDS(old_file), tmp_old)
    render_plot(readRDS(new_file), tmp_new)
    TRUE
  }, error = function(e) {
    add_issue(paste0("plot ", key, ": render failed: ", conditionMessage(e)))
    FALSE
  })
  if (!result) return(invisible())

  if (identical(unname(tools::md5sum(tmp_old)), unname(tools::md5sum(tmp_new)))) {
    return(invisible())
  }

  a <- png::readPNG(tmp_old)
  b <- png::readPNG(tmp_new)
  if (!identical(dim(a), dim(b))) {
    add_issue(paste0("plot ", key, ": rendered image sizes differ"))
    return(invisible())
  }
  diff <- max(abs(a - b))
  if (diff > 0) {
    add_issue(paste0("plot ", key, ": pixels differ (max channel diff ",
                     signif(diff, 3), ")"))
  }
}

old_files <- list_files(old_dir)
new_files <- list_files(new_dir)

missing_new <- setdiff(old_files, new_files)
extra_new <- setdiff(new_files, old_files)
if (length(missing_new) > 0) {
  add_issue(paste0("files missing in new tree: ", paste(missing_new, collapse = ", ")))
}
if (length(extra_new) > 0) {
  add_issue(paste0("extra files in new tree: ", paste(extra_new, collapse = ", ")))
}

for (rel in intersect(old_files, new_files)) {
  compare_file(rel)
}

old_plot_files <- list.files(old_plots, pattern = "\\.plot\\.rds$")
new_plot_files <- list.files(new_plots, pattern = "\\.plot\\.rds$")
missing_plots <- setdiff(old_plot_files, new_plot_files)
extra_plots <- setdiff(new_plot_files, old_plot_files)
if (length(missing_plots) > 0) {
  add_issue(paste0("plots missing in new run: ", paste(missing_plots, collapse = ", ")))
}
if (length(extra_plots) > 0) {
  add_issue(paste0("extra plots in new run: ", paste(extra_plots, collapse = ", ")))
}
for (key in intersect(old_plot_files, new_plot_files)) {
  compare_plot(key)
}

cat("Checked", checks, "files/plots,", length(issues), "failures\n")
if (length(issues) > 0) {
  cat(paste0("  - ", issues, collapse = "\n"), "\n")
  quit(status = 1)
}
