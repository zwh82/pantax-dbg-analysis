#!/usr/bin/env Rscript
# Core-species network analysis entry point.
#
# Reads existing per-cohort network RDS files (g_*_full_network.rds) and writes
# enrichment, node-role, disease-loss and ERAS recovery tables, figures and a
# manuscript summary.  It never rebuilds or changes upstream networks.

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
})
repo_root <- dirname(script_dir)

parse_args <- function(args) {
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i < length(args) && !grepl("^--", args[i + 1])) {
        value <- args[i + 1]
        i <- i + 1
      } else {
        value <- TRUE
      }
    } else {
      i <- i + 1
      next
    }
    if (identical(key, "set")) {
      opts$set <- c(opts$set, value)
    } else {
      opts[[key]] <- value
    }
    i <- i + 1
  }
  opts
}

as_int <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop("--", name, " expects an integer, got: ", value)
  out
}

as_number <- function(value, name) {
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) stop("--", name, " expects a number, got: ", value)
  out
}

parse_cli_value <- function(value) {
  if (is.logical(value)) return(value)
  v <- as.character(value)
  if (v %in% c("NULL", "null", "NA", "na", "none")) return(NULL)
  if (v %in% c("TRUE", "true", "True")) return(TRUE)
  if (v %in% c("FALSE", "false", "False")) return(FALSE)
  number <- suppressWarnings(as.numeric(v))
  if (!is.na(number)) return(number)
  v
}

set_config_field <- function(cfg, path, value) {
  keys <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (any(!nzchar(keys))) stop("Invalid --set path: ", path)
  assign_field <- function(node, keys) {
    if (!is.list(node) || !keys[1] %in% names(node)) {
      stop("Unknown --set path: ", path)
    }
    if (length(keys) == 1) {
      node[[keys[1]]] <- value
    } else {
      node[[keys[1]]] <- assign_field(node[[keys[1]]], keys[-1])
    }
    node
  }
  assign_field(cfg, keys)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(opts$help) || isTRUE(opts$h)) {
  cat(
    "Usage: Rscript core_species/run_core_species.R [options]\n\n",
    "  --config FILE            config file (default core_species/config/core_species.yaml)\n",
    "  --step LIST              analysis, figures, manuscript or all (default: analysis)\n",
    "  --base-dir DIR           directory holding PR*/res_rc networks\n",
    "  --core-rare FILE         core-rare species table\n",
    "  --core-abundant FILE     core-abundant species table\n",
    "  --outdir DIR             analysis result root\n",
    "  --figure-dir DIR         figure output directory\n",
    "  --manuscript-dir DIR     manuscript-statistics output directory\n",
    "  --n-perm 1000            label permutations per network\n",
    "  --n-attack-perm 100      targeted-removal permutations\n",
    "  --seed 123\n",
    "  --set path/to/field=value\n",
    "  --print-config           print the effective config and exit\n",
    "  --save-config FILE       save the effective config, then run\n",
    sep = ""
  )
  quit(save = "no")
}

source(file.path(repo_root, "config", "read_config.R"))
source(file.path(script_dir, "R", "00_utils.R"))

config_path <- if (!is.null(opts$config)) {
  resolve_core_path(opts$config, repo_root)
} else {
  file.path(script_dir, "config", "core_species.yaml")
}
if (!file.exists(config_path)) stop("Config not found: ", config_path)
cfg <- read_config_file(config_path)

if (!is.null(opts$step)) cfg$step <- opts$step
if (!is.null(opts$`base-dir`)) cfg$network_dir <- opts$`base-dir`
if (!is.null(opts$`core-rare`)) cfg$core_rare <- opts$`core-rare`
if (!is.null(opts$`core-abundant`)) cfg$core_abundant <- opts$`core-abundant`
if (!is.null(opts$outdir)) cfg$output$results <- opts$outdir
if (!is.null(opts$`figure-dir`)) cfg$output$figures <- opts$`figure-dir`
if (!is.null(opts$`manuscript-dir`)) cfg$output$manuscript <- opts$`manuscript-dir`
if (!is.null(opts$`n-perm`)) cfg$analysis$n_perm <- as_int(opts$`n-perm`, "n-perm")
if (!is.null(opts$`n-attack-perm`)) {
  cfg$analysis$n_attack_perm <- as_int(opts$`n-attack-perm`, "n-attack-perm")
}
if (!is.null(opts$seed)) cfg$analysis$seed <- as_int(opts$seed, "seed")
if (!is.null(opts$`top-percent`)) {
  cfg$analysis$consensus_hub$top_percent <-
    as_number(opts$`top-percent`, "top-percent")
}
if (!is.null(opts$`min-degree`)) {
  cfg$analysis$consensus_hub$min_degree <-
    as_number(opts$`min-degree`, "min-degree")
}
if (!is.null(opts$`vote-thresh`)) {
  cfg$analysis$consensus_hub$vote_thresh <-
    as_number(opts$`vote-thresh`, "vote-thresh")
}

if (!is.null(opts$set)) {
  for (assignment in opts$set) {
    parts <- strsplit(as.character(assignment), "=", fixed = TRUE)[[1]]
    if (length(parts) < 2) stop("--set expects path=value, got: ", assignment)
    path <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    cfg <- set_config_field(cfg, path, parse_cli_value(value))
  }
}

steps <- trimws(unlist(strsplit(as.character(cfg$step), ",")))
allowed <- c("analysis", "figures", "manuscript", "all")
if (!all(steps %in% allowed)) {
  stop("Unknown --step: ", paste(setdiff(steps, allowed), collapse = ", "),
       ". Use ", paste(allowed, collapse = ", "), ".")
}
if ("all" %in% steps) steps <- c("analysis", "figures", "manuscript")
steps <- unique(steps)

# Resolve every path against the repository root so the recipe does not depend
# on the caller's working directory.
cfg$network_dir <- resolve_core_path(cfg$network_dir, repo_root)
cfg$core_rare <- resolve_core_path(cfg$core_rare, repo_root)
cfg$core_abundant <- resolve_core_path(cfg$core_abundant, repo_root)
cfg$output$results <- resolve_core_path(cfg$output$results, repo_root)
cfg$output$figures <- resolve_core_path(cfg$output$figures, repo_root)
cfg$output$manuscript <- resolve_core_path(cfg$output$manuscript, repo_root)

if (!is.null(opts$`save-config`)) {
  save_path <- opts$`save-config`
  if (isTRUE(save_path)) stop("--save-config expects a file path")
  write_config_file(cfg, save_path)
  message("Saved effective config: ", save_path)
}

if (isTRUE(opts$`print-config`)) {
  if (requireNamespace("yaml", quietly = TRUE)) {
    cat(yaml::as.yaml(config_for_yaml(cfg)))
  } else {
    cat("config <- "); dput(cfg)
  }
  quit(save = "no")
}

if (!is.null(cfg$runtime$library_paths)) {
  .libPaths(c(unlist(cfg$runtime$library_paths), .libPaths()))
}

core_rare_file <- need_file(cfg$core_rare, "Core-rare file")
core_abundant_file <- need_file(cfg$core_abundant, "Core-abundant file")
core_rare <- read_species(core_rare_file)
core_abundant <- read_species(core_abundant_file)
if (length(intersect(core_rare, core_abundant))) stop("Core lists overlap.")

source(file.path(script_dir, "R", "01_networks.R"))
source(file.path(script_dir, "R", "02_edges.R"))
source(file.path(script_dir, "R", "03_permutation.R"))
source(file.path(script_dir, "R", "04_roles.R"))
source(file.path(script_dir, "R", "05_disease.R"))
source(file.path(script_dir, "R", "06_recovery.R"))
source(file.path(script_dir, "R", "07_analysis.R"))

if ("analysis" %in% steps) {
  suppressPackageStartupMessages({
    library(igraph)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(survival)
  })
  run_core_species_analysis(
    base_dir = cfg$network_dir,
    core_rare = core_rare,
    core_abundant = core_abundant,
    outdir = cfg$output$results,
    n_perm = cfg$analysis$n_perm,
    n_attack_perm = cfg$analysis$n_attack_perm,
    seed = cfg$analysis$seed,
    core_rare_file = core_rare_file,
    core_abundant_file = core_abundant_file
  )
}

if ("figures" %in% steps) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(patchwork)
    library(scales)
    library(igraph)
  })
  source(file.path(script_dir, "R", "08_figures.R"))
  run_core_species_figures(cfg$output$results, cfg$output$figures)
  hub <- cfg$analysis$consensus_hub
  run_core_species_consensus_figures(
    cfg$output$results, cfg$output$figures,
    top_percent = hub$top_percent,
    min_degree = hub$min_degree,
    vote_thresh = hub$vote_thresh
  )
}

if ("manuscript" %in% steps) {
  source(file.path(script_dir, "R", "09_manuscript.R"))
  run_core_species_manuscript(cfg$output$results, cfg$output$manuscript)
}

message("Core-species steps complete: ", paste(steps, collapse = ", "))
