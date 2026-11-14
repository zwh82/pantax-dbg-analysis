# Analysis / plotting configuration.
# NOTE: the default config is config/analysis.yaml.  This file is kept as an
# R-format example and is only read when passed explicitly:
#   Rscript analysis/run_analysis.R --config analysis/config/analysis.R
#
# All downstream analysis scripts write their tables and figures under
# `results_dir` and read cohort networks from `network_dir` (NULL means
# results_dir, which reproduces the legacy layout).
#   Rscript analysis/run_analysis.R --config analysis/config/analysis.R

analysis_config <- list(

  # Root that contains one directory per cohort, plus the derived
  # net_compare2_V3/, consensus_hub/, stab/ and scripts/final_figures/ trees.
  results_dir = "/home/work/wenhai/meta/high_rare_net4",

  # Cohort network root; NULL keeps the legacy single-root layout.
  network_dir = NULL,

  # One or more network result suffixes ("res_rc", "res_flashweave_rc", ...).
  # The first suffix is used by the steps that historically took a single suffix.
  suffix = "res_rc",

  # Prevalence used in the species-prevalence file names written by the
  # pipeline.  It is part of the file name even when the filter was skipped.
  prevalence = 0.1,

  # Disease labels that each step skips.  These reproduce the per-script
  # `next` conditions in the legacy code.
  skip_diseases = list(
    stability   = c("Health"),
    node        = c("Health"),
    degree      = character(0),
    net_compare = c("Health"),
    rho         = c("ERAS", "Health"),
    rho_edge    = c("ERAS", "Health"),
    hub         = c("ERAS")
  ),

  # Consensus-hub definition shared by both hub scripts.
  hub = list(
    top_percent    = 0.10,
    min_degree     = 3,
    vote_threshold = 3,
    seed           = 123,
    levels         = c("rare", "abundant"),
    colors         = c(abundant = "#e34a33", rare = "#2b8cbe"),
    step_within    = 0.65,
    step_between   = 1.25,
    bar_width      = 0.55,
    figure_height  = 11
  ),

  # Output sub-directories, relative to results_dir.  The names are part of
  # the published layout and must not be changed without re-running the
  # figure assembly.
  paths = list(
    stability    = "stab",
    node_overlap = "net_compare/node_overlap",
    degree       = "mean_degree_change",
    net_compare  = "net_compare2_V3",
    rho          = "rho_compare",
    hub          = "consensus_hub",
    final        = "scripts/final_figures",
    final_split  = "scripts/final_figures/main_split"
  ),

  # Disease-name mapping used by the analysis tables and plot labels.
  # NULL means <results_dir>/scripts/disease_name_mapping.tsv.
  disease_map_file = NULL,

  # The optional hub-summary plot accepts an externally produced summary
  # table.  These defaults reproduce the legacy CLI defaults; pass --input
  # and --outdir to point at another table.
  hub_summary = list(
    input   = "/home/work/meishi/Themis/module_hub/outputs_net4/ConsensusHub_master_summary.csv",
    out_dir = "/home/work/meishi/Themis/module_hub/outputs_net4/Meta_ConsensusHub_Abundance_res_rc_net4"
  ),

  # ERAS cohort used by the final main figure.  Its network result directory
  # supplies the centrality and edge-count inputs.
  eras_cohort = "PRJEB20800_57_ERAS",

  runtime = list(
    library_paths = c(
      "/home/wenhai/miniconda3/envs/r_4.3/lib/R/library",
      "/home/wenhai/application/Qmatey/tools/R"
    )
  )
)
