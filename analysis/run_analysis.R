#!/usr/bin/env Rscript
# Entry point for the downstream analysis and figure scripts.
#
# Steps read the network result tree produced by run_pipeline.R and write
# tables and figures back into the same tree:
#
#   Rscript analysis/run_analysis.R --step all
#   Rscript analysis/run_analysis.R --step net_compare,rho,rho_edge
#   Rscript analysis/run_analysis.R --step hub_summary \
#     --input /path/ConsensusHub_master_summary.csv --outdir /path/out

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
      opts[[key]] <- sub("^--[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i < length(args) && !grepl("^--", args[i + 1])) {
        opts[[key]] <- args[i + 1]
        i <- i + 1
      } else {
        opts[[key]] <- TRUE
      }
    }
    i <- i + 1
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

resolve_input_path <- function(path) {
  if (is.null(path) || grepl("^/", path)) path else file.path(script_dir, path)
}

config_path <- resolve_input_path(
  if (!is.null(opts$config)) opts$config else "config/analysis.yaml"
)
if (!file.exists(config_path)) stop("Config not found: ", config_path)

source(file.path(script_dir, "..", "config", "read_config.R"))
cfg <- read_config_file(config_path)

if (!is.null(opts$results)) cfg$results_dir <- opts$results
if (!is.null(opts$`network-dir`)) cfg$network_dir <- opts$`network-dir`
if (is.null(cfg$network_dir) || !nzchar(cfg$network_dir)) {
  cfg$network_dir <- cfg$results_dir
}
if (!is.null(opts$suffix)) cfg$suffix <- trimws(unlist(strsplit(opts$suffix, ",")))
if (!is.null(opts$prevalence)) cfg$prevalence <- as.numeric(opts$prevalence)
if (!is.null(opts$`disease-map`)) cfg$disease_map_file <- opts$`disease-map`

if (!is.null(cfg$runtime$library_paths)) {
  .libPaths(c(cfg$runtime$library_paths, .libPaths()))
}

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(stringr)
  library(ggsci)
  library(scales)
  library(gridExtra)
  library(clinfun)
  library(RColorBrewer)
})

source(file.path(script_dir, "R/00_utils.R"))
source(file.path(script_dir, "R/01_stability.R"))
source(file.path(script_dir, "R/02_node_overlap.R"))
source(file.path(script_dir, "R/03_degree_change.R"))
source(file.path(script_dir, "R/04_net_compare.R"))
source(file.path(script_dir, "R/05_rho.R"))
source(file.path(script_dir, "R/06_consensus_hub.R"))
source(file.path(script_dir, "R/07_final_figure.R"))

ensure_dir(cfg$results_dir)
validate_analysis_config(cfg)
ensure_dir(file.path(cfg$results_dir, cfg$paths$final))

steps <- c(
  "stability", "node", "degree", "net_compare",
  "rho", "rho_edge", "hub", "hub_summary", "final"
)

requested <- if (is.null(opts$step)) "all" else opts$step
requested <- trimws(unlist(strsplit(requested, ",")))
if (identical(requested, "all")) {
  requested <- setdiff(steps, "hub_summary")
}
unknown <- setdiff(requested, steps)
if (length(unknown) > 0) {
  stop("Unknown step(s): ", paste(unknown, collapse = ", "),
       "\nAvailable: ", paste(steps, collapse = ", "))
}

run_step <- function(step) {
  message("\n=== ", step, " ===")
  switch(
    step,
    stability   = run_stability(cfg),
    node        = run_node_overlap(cfg),
    degree      = run_degree_change(cfg),
    net_compare = run_net_compare(cfg),
    rho         = run_rho_global(cfg),
    rho_edge    = run_rho_by_edge_type(cfg),
    hub         = run_meta_consensus_hub(cfg),
    hub_summary = run_hub_summary_plot(
      cfg,
      in_file = if (!is.null(opts$input)) opts$input else NULL,
      out_dir = if (!is.null(opts$outdir)) opts$outdir else NULL
    ),
    final       = run_final_figure(cfg)
  )
}

for (step in requested) {
  run_step(step)
}

message("\nDone: ", paste(requested, collapse = ", "))
