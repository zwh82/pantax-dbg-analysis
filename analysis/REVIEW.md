# Analysis refactor review

Scope: the nine downstream analysis/figure scripts listed in
`scripts/final_figures/README.md`. The original `high_rare_net4` tree was not
modified.

## 1. File mapping

| Old | New |
|---|---|
| `stab_plot.R` | `R/01_stability.R` |
| `node_compare.R` | `R/02_node_overlap.R` |
| `degree_change_analysis_V2.R` | `R/03_degree_change.R` |
| `net_compare2_V3.R` | `R/04_net_compare.R` |
| `rho_plot.R` | `R/05_rho.R::run_rho_global` |
| `rho_edge_plot.R` | `R/05_rho.R::run_rho_by_edge_type` |
| `meta_consensushub_abundance_count_V3.R` | `R/06_consensus_hub.R::run_meta_consensus_hub` |
| `plot_consensus_hub_rare_abundant.R` | `R/06_consensus_hub.R::run_hub_summary_plot` |
| `final_main_figure_from_rds_V1.R` | `R/07_final_figure.R` |
| `02_utils.R` (selected helpers) | `R/00_utils.R` |
| — | `run_analysis.R`, `config/analysis.yaml` (new entry/config) |

3,594 legacy lines became 2,895 lines including the new entry and config.
The reduction comes from deleting commented-out alternatives, dead
computations, debug prints, and merging the shared hub plotting code.

## 2. Behavior-preservation proof

The legacy scripts and `run_analysis.R` were run on the same synthetic results
tree, then every output file was compared. Plot objects were captured at each
`ggsave()` call, re-rendered with a fixed seed, and compared pixel-for-pixel.

| Step | Files/plots compared | Result |
|---|---|---|
| stability | — | 0 failures |
| node overlap | — | 0 failures |
| degree change | — | 0 failures |
| edge retention (7 diseases + ERAS) | — | 0 failures |
| rho / rho by edge type | — | 0 failures |
| consensus hub (meta + summary) | — | 0 failures |
| final figure + split panels | — | 0 failures |

The full comparison reported 0 failures across 288 files/plots. Fixture and
harness details are in `verification/README.md`.

## 3. Behavior changes

1. **One entry point.** Nine standalone scripts became steps of
   `run_analysis.R`. Output paths and file names are unchanged.
2. **Separate network and analysis roots (optional).** `network_dir` (CLI
   `--network-dir`) points at the cohort networks while `results_dir` receives
   the analysis outputs. When it is `NULL` both are read/written under
   `results_dir`, reproducing the legacy layout exactly.
3. **`hub_summary` CLI.** The two positional arguments became `--input` and
   `--outdir`; the legacy defaults are kept in `config/analysis.yaml`.
4. **`pheatmap` is optional.** The legacy scripts failed at `library(pheatmap)`
   when the package was absent. The refactor keeps the heatmap branch when
   `pheatmap` is installed and skips it with a warning otherwise. When
   installed, the calls and outputs are identical.
5. **Invariant failure stops instead of quitting R.** The
   `PreservedRate + LossRate != 1` branch called `q()`. It now calls `stop()`
   with the same message. The invariant is mathematically guaranteed for
   valid input.
6. **Dead code removed** (no output change): the unused preserved/lost
   `p_combined`, `kw_result`, pairwise-Wilcoxon prints, dead `All` branch in
   `compare_network_edges`, dead node heatmap preparation, unused
   `make_loss_panel` in the final-figure script, and unused `log_p` column.
7. **Global theme restored.** The final-figure script used to call
   `theme_set()` and leave it changed. The refactor restores the previous
   theme on exit; within the step the plots are identical.
8. **Consensus-hub scan uses one suffix.** The legacy script looped over its
   suffix list with a fixed `TOOLS = "res_rc"`; scanning once with the first
   configured suffix produces the same `master_df`.
9. **Degree change loops over configured suffixes.** The legacy script used
   one hard-coded suffix (extension only; the default is `res_rc`).
10. **Disease mapping read.** `read_tsv(col_names = ...)` became
   `read.delim(header = FALSE, col.names = ...)`, preserving the "no header"
   layout of `disease_name_mapping.tsv`.
11. **No stray `Rplots.pdf`.** `grid.arrange()` also draws to the active
    device, so the legacy run left an `Rplots.pdf` in the working directory.
    The refactor uses `arrangeGrob()`, which returns the same gtable without
    drawing; the saved heatmap is unchanged.

## 4. code-humanizer scan (Mode A)

### Fixed

| Pattern | Location | Action |
|---|---|---|
| 1. Reimplemented helper | disease-row selection in both rho scripts | Extracted `pick_disease_row()` |
| 1. Reimplemented helper | first-title-per-row in both rho scripts | Extracted `keep_first_title_each_row()` |
| 1. Reimplemented helper | Rich -> Abundant relabeling | Extracted `rename_rich_to_abundant()` |
| 1. Reimplemented helper | hub positions / stacked hub plot duplicated in two scripts | Shared `make_hub_positions()`, `build_hub_grouped_plot()`, `hub_figure_width()` |
| 5. Dead code | preserved/lost plot overwritten before use | Removed |
| 5. Dead code | Kruskal-Wallis result, pairwise Wilcoxon block, `log_p`, control-network summaries | Removed |
| 5. Dead code | unused final-figure `make_loss_panel` | Removed |
| 5. Dead code | node heatmap data preparation in `node_compare.R` | Removed |
| 12/14. Noise | `print()` dumps of paths, tables, and intermediate results | Removed |
| 14. Dead imports | `rstatix`, `purrr`, `reshape2`, top-level `pheatmap` | Removed / namespaced |
| 7. Config sprawl | hub thresholds, colors, spacing, paths | Moved to `config/analysis.yaml` |

### Reported, not fixed

| Pattern | Location | Why left alone |
|---|---|---|
| 14. Destructive cleanup | `unlink(out_base/<disease>)` for single-cohort diseases | Pre-existing behavior; removing it changes the output tree. Recommend deleting this branch after confirming no workflow relies on it. |
| 14. Abrupt failure | `stop()` when a disease has cohorts but no network files | Pre-existing; keeps the run from silently producing partial figures. |
| Leaked loop variable | stability's `np_overall_<last method>` / `fig_net_1` | Deliberately reproduced: the legacy output always uses the last configured method. Documented in `R/01_stability.R`. |
| 14. Hidden constant | edge weight threshold `0.2` in `compare_network_edges` | Pre-existing analysis definition; moving it to config without a re-run would be a silent analysis change. Candidate for `config`. |
| 14. Hidden constants | figure sizes, axis limits, palettes | Pre-existing publication layout; many are tuned to the manuscript. |
| 1. Duplication | hub data preparation in meta vs summary scripts | They consume different inputs (networks vs a CSV); only the plotting half is shared. |
| 8. Defensive code | `tryCatch()` around hub reads | One failed cohort should not abort the whole scan; the legacy behavior was a message and skip. |
| Naming contract | disease label = last `_` token of the cohort directory | Pre-existing; works for the current cohort naming. |
| 15/16. Test gap | No in-repo tests | Verification is an external golden run; see `verification/README.md`. |

## 5. Deliberately preserved legacy behavior

- Output directory names, file names (including
  `<suffix>multiple_cohort_combined.png` without an underscore), and plot sizes.
- `fig_net_1`'s use of the last method for `np_overall` and the density panel.
- The per-script disease skip lists (`Health` vs `ERAS` vs both).
- `net_compare` removing the output directory of diseases with fewer than two
  cohorts.
- ERAS time points `D4, D8, D42, D180` and the `Control -> D0` relabeling in
  the final figure.
- Rich/Rare labels are renamed to Abundant/Rare only in the places the legacy
  code renamed them.
