# Refactor review

Scope: the code reachable from `scripts/00_main.R` only. The original
`high_rare_net4` tree was not modified.

## 1. File mapping

| Old | New |
|---|---|
| `00_main.R` | `run_pipeline.R` |
| `01_config.R` | `config/pipeline.yaml` (R format still accepted via `--config`) |
| `02_utils.R` | `R/01_abundance.R` |
| `03_functions_network.R` | `R/02_network.R` |
| `04_functions_analysis.R` | `R/03_hub_analysis.R` |
| `05_functions_plot.R` | `R/04_edge_plots.R`, `R/05_centrality_plots.R`, `R/06_network_plot.R` |
| `06_process_one.R` | `R/07_process_cohort.R` |
| — | `R/00_config.R` (new: cohort table, output paths, validation) |

Entry-reachable code shrank from 3,513 lines to 2,618 lines. The reduction
comes from deleting unreachable functions, commented-out alternatives, and
duplicated inline blocks, not from dropping a mode.

## 2. Behavior-preservation proof

Each case ran the old code and the new code on the same synthetic input, then
compared every output file:

- graph `.rds`: vertex table, edge table, and graph attributes, sorted by
  name; raw RDS bytes are not compared because igraph serializes internal
  caches differently across identical runs.
- `.csv`/`.tsv`: parsed table equality with `all.equal`.
- `.txt`: exact text equality.
- `fig_rds/*.rds`: ggplot labels, mapping, and `ggplot_build()$data`.
- `.png`: skipped (device-dependent bytes).

A determinism control (old code run twice) passes with 0 failures.

| Case | Files | Result |
|---|---|---|
| read count, SPIEC-EASI, 10% prevalence | 36 | 0 failures |
| fixed GML input, multi-group plots | 40 | 0 failures |
| skip prevalence filter | 36 | 0 failures |
| relative-abundance input | 40 | 0 failures |
| network-only mode | 6 | 0 failures |
| merged-cohort mode | 36 | 0 failures |
| edge-weight filter | 4 | 0 failures |
| edge-confidence filter | 4 | 0 failures |
| SparCC | 43 | 0 failures |
| FlashWeave | 41 | 0 failures |

The comparator is `verification/compare_outputs.R`. The synthetic fixtures are
not shipped; `verification/README.md` records how they were built.

## 3. Behavior changes

These are not pure refactors and are called out explicitly.

1. **Merged-cohort mode repaired.** The legacy branch could not run: it used
   `project_name`, `all_groups`, `group_abundance_categories`, and `cohort`
   without defining them, and `metadata$Disease <- ...` crashed when the
   metadata had no `Disease` column. The refactor defines the four values and
   drops the unused `Disease` assignment. Verified against a locally patched
   legacy copy with the same repairs (36/36 outputs match). The original
   merged configs are commented out in `01_config.R`, so this branch was not
   in active use.
2. **Cohort table validation.** Empty tables now stop with an error; a
   6-column relative-abundance table is accepted (the old entry only accepted
   5 columns in that mode). Valid legacy inputs are unchanged.
3. **No `setwd()`.** Config, module, and cohort paths resolve from the entry
   script's directory. File contents are unchanged.
4. **Package loading.** `library(tidyverse)` and the unused `library(readr)`
   were replaced by explicit package loads. The used functions come from the
   same packages. This also avoids a local environment where the `tidyverse`
   metapackage fails while its components load.
5. **Merged mode `metadata$Disease` assignment removed.** It was unused and
   could crash; see item 1.
6. **Console output.** Debug dumps (`positive_edges`, `edge_count_summary`,
   `sorted_names`, `isnull`, `Group2`, `print(g$Stability)`) and threshold
   numbers in two log lines were removed. Output files are unchanged.
7. **Abundance legend labels moved to config.** Defaults reproduce the old
   strings exactly.

## 4. code-humanizer scan (Mode A)

Catalog references are from `code-humanizer` SKILL.md.

### Fixed during the refactor

| Pattern | Location | Severity | Action |
|---|---|---|---|
| 1. Reimplemented helper | Rich/Rare mapping in edge and centrality modules | 2 | Extracted `abundance_type_from_category()` |
| 1. Reimplemented helper | Three inline `fread` + rowname blocks | 2 | Extracted `read_matrix_file()` |
| 1. Reimplemented helper | Prevalence filtering repeated in both input branches | 2 | Extracted `filter_by_prevalence()` |
| 1. Reimplemented helper | Directory setup, metric table, stability table | 2 | Extracted helpers |
| 5. Dead code | `plot_edges_connect` (called nonexistent `plot_edges_connect_multi`) | 3 | Removed; unreachable |
| 5. Dead code | `plot_edge_type_heatmap` (unreachable; needs unloaded `pheatmap`) | 3 | Removed |
| 5. Dead code | `analyze_network_collapse`, `assess_robustness` | 3 | Removed from the entry copy; still in the original repo |
| 5. Dead code | `read_matrix`, `get_disease_cohorts`, `replace_disease_name*` | 2 | Removed from the entry copy; downstream scripts still use the original |
| 5. Dead code | `abundance_matrix_full` parameter, `n_high`/`n_low`/`n_very_low`, `base_group`, duplicate `all_edges_list` | 2 | Removed |
| 12/14. Noise | Debug `print()`/`cat()` and a decorative banner | 2 | Removed |
| 14. Dead imports | `tidyverse`, `readr` | 2 | Removed |

### Reported, not fixed

| Pattern | Location | Severity | Why left alone |
|---|---|---|---|
| 7. Config/API sprawl | ACA group-order branches in `R/04_edge_plots.R` and `R/05_centrality_plots.R` | 3 | Three places use inconsistent orders (`Control,ACA,CRC` vs `Control,CRC,ACA`). Changing them changes legend order. Recommend a `group_levels` config field in a follow-up. |
| 6. Wrapper that adds nothing | `read_metadata_file()` | 1 | One-line wrapper, used twice; kept as a named input seam. |
| 1. Duplication | `plot_multi_edges()` and `plot_multi_edges_connect()` share color/Rho setup | 2 | They produce different outputs; merging risks behavior drift. |
| 1. Duplication | `R/04_edge_plots.R::classify_edges()` and `net_compare2_V3.R::get_edge_df()` | 2 | Cross-scope duplication with a downstream, out-of-scope script. |
| 8. Defensive code | `tryCatch()` around `mean_distance()` | 1 | Justified for disconnected networks; only wraps one call. |
| 14. Abrupt failure | `quit(1)` when no graph source is configured | 2 | Pre-existing failure behavior; changing it changes the exit path. Recommend `stop()` later. |
| 14. Hidden constant | Edge-strength breaks `0, 0.1, 0.3, 0.5, 1` | 2 | Pre-existing plot definition; config candidate. |
| 15/16. Test gap | No in-repo test suite | 3 | Verification is an external golden run. Add characterization tests if this pipeline is maintained. |

## 5. Deliberately preserved legacy behavior

- `runs$eras` remains a global switch. The 73-cohort table excludes the ERAS
  cohort; run ERAS separately with its own table and `--eras TRUE`.
- `input$prevalence` still appears in output file names even when the filter is
  skipped (`Notprefilter` behavior).
- `output$suffix = NULL` derives the old directory names, including the
  `res_raclr_econ_filter` name for the confidence filter in read-count mode.
- ACA and ERAS factor orders are unchanged.
- The `config$abundance$labels` default strings match the old legend text.

## 6. Downstream analysis and figures

The nine final-figure/analysis scripts are refactored separately under
`analysis/`, with their own entry point (`analysis/run_analysis.R`), config
(`analysis/config/analysis.yaml`), and review (`analysis/REVIEW.md`). They were
verified with the same old-vs-new method: 288 files/plots, 0 failures. The
original tree was not modified.

## 7. Parameterised entry point

All switches from `01_config.R` are command-line flags, so a run no longer
requires uncommenting/copying a config block:

- The default config is `config/pipeline.yaml`; `--config` also accepts the
  legacy R-list format (`config/pipeline.R` is kept as an example).
- `config/eras.yaml` is the persistent ERAS variant: it points at
  `config/eras.tsv` and sets `runs.eras: true`, so the ERAS run needs no
  `--eras TRUE` flag. `run.sh eras` and `run.sh full` use it.
- `--data-type`, `--tool`, `--prevalence`, `--present`, `--floor`,
  `--fill-zero-counts`, `--skip-prevalence-filter`, `--filter-weight`,
  `--filter-confidence`, `--network-only`, `--core-network-only`,
  `--prjna398089-ibd-not`, `--spiec-easi-ra`, `--suffix`, `--flag`,
  `--merged`, `--eras`, `--cohorts`, `--outdir`, `--only`,
  `--hub-top-percent`, `--hub-min-degree`, `--seed`.
- `--set path/to/field=value` covers any nested field, e.g.
  `--set network/spiec_easi/pulsar/rep.num=50`.
- `--print-config` shows the effective config plus the derived output
  directory/suffix; `--save-config FILE` writes it in a form that can be
  passed back via `--config` and then continues the run.
- `run.sh` contains the common combinations as commented scenarios. It passes
  every option explicitly, writes the command and log to `runs/logs/`, and
  saves the effective config to `runs/<scenario>/configs/`. `DRY_RUN=1`
  prints the commands without executing them.
