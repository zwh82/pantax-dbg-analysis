#!/usr/bin/env Rscript
# Entry point for the high_rare_net4 network pipeline.

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
})

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

as_bool <- function(value, name) {
  if (is.logical(value)) return(value)
  v <- tolower(as.character(value))
  if (v %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (v %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("--", name, " expects TRUE or FALSE, got: ", value)
}

as_number <- function(value, name) {
  number <- suppressWarnings(as.numeric(value))
  if (is.na(number)) stop("--", name, " expects a number, got: ", value)
  number
}

as_nullable_string <- function(value) {
  if (is.logical(value) && isTRUE(value)) return(NULL)
  if (is.null(value) || is.na(value)) return(NULL)
  if (tolower(as.character(value)) %in% c("", "auto", "none", "null", "na")) return(NULL)
  as.character(value)
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

# Set a nested config field. Paths use "/" so that field names containing a
# dot (rep.num, subsample.ratio, ...) stay unambiguous.
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
    "Usage: Rscript run_pipeline.R [options]\n\n",
    "  --config FILE                 config file (default config/pipeline.yaml)\n",
    "  --cohorts FILE                cohort table (overrides config)\n",
    "  --outdir DIR                  output root (overrides config)\n",
    "  --only COHORT[,COHORT]        process only these cohorts\n",
    "  --data-type read_count|relative_abundance\n",
    "  --tool spiec-easi|flashweave|sparcc\n",
    "  --prevalence 0.1              prevalence threshold and file-name value\n",
    "  --present 0                   present threshold\n",
    "  --floor 0                     values below this become zero\n",
    "  --fill-zero-counts TRUE|FALSE\n",
    "  --skip-prevalence-filter TRUE|FALSE\n",
    "  --filter-weight TRUE|FALSE    |weight| >= 0.1 edge filter\n",
    "  --filter-confidence TRUE|FALSE\n",
    "  --graph-file FILE|none        read pre-built GML instead of fitting\n",
    "  --merged TRUE|FALSE           merged-cohort table layout\n",
    "  --eras TRUE|FALSE             ERAS time ordering\n",
    "  --network-only TRUE|FALSE     build Control networks only\n",
    "  --core-network-only TRUE|FALSE\n",
    "  --prjna398089-ibd-not TRUE|FALSE\n",
    "  --spiec-easi-ra TRUE|FALSE    relative-abundance SpiecEasi\n",
    "  --suffix NAME|auto            output suffix (auto = derive)\n",
    "  --flag TEXT                   output flag\n",
    "  --hub-top-percent 0.1\n",
    "  --hub-min-degree 5\n",
    "  --seed 123\n",
    "  --set path/to/field=value     any config field, e.g.\n",
    "                                --set network/spiec_easi/pulsar/rep.num=50\n",
    "  --print-config                print the effective config and exit\n",
    "  --save-config FILE            save the effective config, then run\n",
    sep = ""
  )
  quit(save = "no")
}

resolve_input_path <- function(path) {
  if (is.null(path) || grepl("^/", path)) path else file.path(script_dir, path)
}

config_path <- resolve_input_path(
  if (!is.null(opts$config)) opts$config else "config/pipeline.yaml"
)
if (!file.exists(config_path)) stop("Config not found: ", config_path)

source(file.path(script_dir, "config", "read_config.R"))
cfg <- read_config_file(config_path)

# ---- command-line overrides -------------------------------------------------
# Every switch below corresponds to a field in config/pipeline.yaml (or the
# config passed with --config).  CLI values take precedence over the config
# file, so no config file editing is needed for a one-off run.
if (!is.null(opts$`data-type`)) cfg$input$data_type <- opts$`data-type`
if (!is.null(opts$merged)) cfg$input$merged <- as_bool(opts$merged, "merged")
if (!is.null(opts$eras)) cfg$runs$eras <- as_bool(opts$eras, "eras")
if (!is.null(opts$outdir)) cfg$output$dir <- opts$outdir
if (!is.null(opts$cohorts)) cfg$cohorts_file <- opts$cohorts

if (!is.null(opts$tool)) {
  tool <- tolower(opts$tool)
  if (tool %in% c("spieceasi", "spiec_easi")) tool <- "spiec-easi"
  cfg$network$tool <- tool
}
if (!is.null(opts$prevalence)) cfg$input$prevalence <- as_number(opts$prevalence, "prevalence")
if (!is.null(opts$present)) cfg$input$present <- as_number(opts$present, "present")
if (!is.null(opts$floor)) cfg$input$floor <- as_number(opts$floor, "floor")
if (!is.null(opts$`fill-zero-counts`)) {
  cfg$input$fill_zero_counts <- as_bool(opts$`fill-zero-counts`, "fill-zero-counts")
}
if (!is.null(opts$`skip-prevalence-filter`)) {
  cfg$input$skip_prevalence_filter <- as_bool(opts$`skip-prevalence-filter`, "skip-prevalence-filter")
}
if (!is.null(opts$`filter-weight`)) {
  cfg$network$filter_by_weight <- as_bool(opts$`filter-weight`, "filter-weight")
}
if (!is.null(opts$`filter-confidence`)) {
  cfg$network$filter_by_confidence <- as_bool(opts$`filter-confidence`, "filter-confidence")
}
if (!is.null(opts$`network-only`)) {
  cfg$runs$network_only <- as_bool(opts$`network-only`, "network-only")
}
if (!is.null(opts$`core-network-only`)) {
  cfg$runs$core_network_only <- as_bool(opts$`core-network-only`, "core-network-only")
}
if (!is.null(opts$`prjna398089-ibd-not`)) {
  cfg$runs$prjna398089_ibd_not <- as_bool(opts$`prjna398089-ibd-not`, "prjna398089-ibd-not")
}
if (!is.null(opts$`spiec-easi-ra`)) {
  cfg$network$spiec_easi$relative_abundance <-
    as_bool(opts$`spiec-easi-ra`, "spiec-easi-ra")
}
if (!is.null(opts$`graph-file`)) cfg$network$graph_file <- as_nullable_string(opts$`graph-file`)
if (!is.null(opts$suffix)) cfg$output$suffix <- as_nullable_string(opts$suffix)
if (!is.null(opts$flag)) cfg$output$flag <- opts$flag
if (!is.null(opts$`hub-top-percent`)) {
  cfg$hub$top_percent <- as_number(opts$`hub-top-percent`, "hub-top-percent")
}
if (!is.null(opts$`hub-min-degree`)) {
  cfg$hub$min_degree <- as_number(opts$`hub-min-degree`, "hub-min-degree")
}
if (!is.null(opts$seed)) cfg$runtime$seed <- as_number(opts$seed, "seed")

# Generic escape hatch for any other field, e.g.
#   --set abundance/thresholds/high=2e-3
#   --set network/spiec_easi/pulsar/rep.num=50
if (!is.null(opts$set)) {
  for (assignment in opts$set) {
    parts <- strsplit(as.character(assignment), "=", fixed = TRUE)[[1]]
    if (length(parts) < 2) stop("--set expects path=value, got: ", assignment)
    path <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    cfg <- set_config_field(cfg, path, parse_cli_value(value))
  }
}

# --save-config persists the effective config and then continues with the run,
# so a recipe can always leave a reproducible config next to its results.
# --print-config only inspects the effective config and exits.
if (!is.null(opts$`save-config`)) {
  save_path <- opts$`save-config`
  if (isTRUE(save_path)) stop("--save-config expects a file path")
  write_config_file(cfg, save_path)
  message("Saved effective config: ", save_path)
}

if (isTRUE(opts$`print-config`)) {
  source(file.path(script_dir, "R/00_config.R"))
  validate_config(cfg)
  cat("# resolved output dir: ", resolve_output_dir(cfg), "\n", sep = "")
  cat("# resolved suffix: ", resolve_output_suffix(cfg), "\n", sep = "")
  if (requireNamespace("yaml", quietly = TRUE)) {
    cat(yaml::as.yaml(config_for_yaml(cfg)))
  } else {
    cat("config <- ")
    dput(cfg)
  }
  quit(save = "no")
}

if (!is.null(cfg$runtime$library_paths)) {
  .libPaths(c(cfg$runtime$library_paths, .libPaths()))
}

se_lib <- if (isTRUE(cfg$network$spiec_easi$relative_abundance)) {
  cfg$runtime$spiec_easi_library$relative_abundance
} else {
  cfg$runtime$spiec_easi_library$default
}

suppressPackageStartupMessages({
  library(SpiecEasi, lib.loc = se_lib)
  library(igraph)
  library(RColorBrewer)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  library(tidyr)
  library(cowplot)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(Matrix)
  library(ggsci)
})

source(file.path(script_dir, "R/00_config.R"))
source(file.path(script_dir, "R/01_abundance.R"))
source(file.path(script_dir, "R/02_network.R"))
source(file.path(script_dir, "R/03_hub_analysis.R"))
source(file.path(script_dir, "R/04_edge_plots.R"))
source(file.path(script_dir, "R/05_centrality_plots.R"))
source(file.path(script_dir, "R/06_network_plot.R"))
source(file.path(script_dir, "R/07_process_cohort.R"))

validate_config(cfg)

cohorts_path <- resolve_input_path(cfg$cohorts_file)
config_data <- read_cohort_table(cohorts_path, merged = cfg$input$merged)

if (!is.null(opts$only)) {
  keep <- trimws(unlist(strsplit(opts$only, ",")))
  key_col <- if (cfg$input$merged) "project_name" else "cohort"
  config_data <- config_data[config_data[[key_col]] %in% keep, , drop = FALSE]
  if (nrow(config_data) == 0) stop("No cohorts matched --only: ", opts$only)
}

cat("Total diseases to process:", nrow(config_data), "\n")

for (i in 1:nrow(config_data)) {
  process_one_cohort(i, config_data, cfg)
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("All done!\n")
cat(rep("=", 60), "\n", sep = "")
