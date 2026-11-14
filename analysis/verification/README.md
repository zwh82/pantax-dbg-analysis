# Verification

`compare_analysis.R` compares two result trees plus the plots captured by
`prelude.R`:

```bash
Rscript analysis/verification/compare_analysis.R \
  <old_results> <new_results> <old_plots> <new_plots>
```

- tables: parsed and compared by value (`all.equal`, 1e-12)
- RDS: data-frame columns, rows, factor levels, list structure
- text: exact lines
- recorded `ggsave()` plots: rendered with `set.seed(1)` at 6x6 in / 100 dpi
  and compared pixel-for-pixel (`png::readPNG`)
- PNG/PDF files: skipped (device-dependent)

A run with 0 failures is the acceptance criterion.

`prelude.R` captures each `ggsave()` plot object while still writing the
requested file:

```bash
PLOT_RECORD_DIR=/tmp/plots_old PLOT_RECORD_BASE=/tmp/results_old \
  Rscript -e 'source("analysis/verification/prelude.R"); source("<legacy script>")'
```

## Fixture

The synthetic results tree was built outside the repository from the
SPIEC-EASI and GML fixtures used for the pipeline verification:

1. Seven diseases with two cohorts each (`PRJTEST_A_TestDisease` ...
   `PRJTEST_N_Alzheimer`). `Alzheimer` is present in
   `disease_name_mapping.tsv`, so the label mapping (`Alzheimer -> AD`) is
   exercised.
2. Each cohort has `res_rc/g_Control_full_network.rds`,
   `g_<Disease>_full_network.rds`, `data/stability_df.csv`,
   `data/rho_results_multi.csv`, `data/rho_results_edge_type.csv`, and the two
   prevalence TSVs used by `stab_plot.R`.
3. Edge subsets differ per cohort (balanced deletion across Rich-Rich,
   Rare-Rare and Rich-Rare) so Friedman/Jonckheere and degree-change tests run
   on non-degenerate data.
4. `PRJEB20800_57_ERAS` has `g_Control` plus `g_D4/D8/D42/D180`, the real
   `data/combined_data.rds` and `data/edge_count_summary.rds`, and its own
   stability/rho/prevalence inputs. This exercises the ERAS branch of
   `net_compare2_V3.R` and the time-point panels of the final figure.
5. `ConsensusHub_master_summary.csv` is a small per-cohort table for the
   standalone hub-summary plot.

## Environment

`pheatmap` and `rstatix` are not installed in the R 4.3 environment used for
the run. The legacy scripts import both. Verification therefore installed two
minimal stub packages into a temporary library:

- `pheatmap`: returns a gtable-compatible object (the only use is
  `p_heatmap2[[4]]` for the control-network heatmap);
- `rstatix`: empty namespace (imported but never called).

Both the legacy run and the refactored run see the same stubs, so the
comparison is still valid. With a real `pheatmap` installed, the refactored
code calls `pheatmap::pheatmap()` exactly as the legacy code does.

## Final figure

`final_main_figure_from_rds_V1.R` is compared against the published input
objects (`net_compare2_V3/combined/res_rc_main_data.rds`,
`consensus_hub/Meta_ConsensusHub_Abundance_plot_data.rds`, the ERAS
`combined_data.rds` / `edge_count_summary.rds`, and
`net_compare2_V3/ERAS/res_rc/Preserved_zhexian.rds`) copied into both trees
before the final step. This guarantees the six-panel layout, the real hub
table, and the five ERAS time points are covered. The fixture-produced
`net_compare`/`hub` RDS (including the ERAS preserved table) are restored
afterwards so those objects are also part of the comparison.

## Cases

| Case | Result |
|---|---|
| stability + sample-size figure | 0 failures |
| node overlap + bar plot | 0 failures |
| degree change (multi/single-cohort) | 0 failures |
| edge retention, 7 diseases + ERAS | 0 failures |
| rho (global and by edge type) | 0 failures |
| consensus hub (meta + summary table) | 0 failures |
| final composite figure + split panels | 0 failures |

The plots compared include `fig_net_1` ... `fig_net_6`, `fig_net_4m`,
`fig_net_main_2`, every `main_split` panel, the per-disease and combined
stability/rho/edge figures, and the hub MAIN/SUPP figures.
