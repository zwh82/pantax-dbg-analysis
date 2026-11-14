<p align="center">
  <strong>Language:</strong> <a href="README.md">English</a> | <a href="README.zh.md">Chinese</a>
</p>

# pantax-dbg-analysis

Analysis workflows for profiling data produced by `pantax-dbg`. The repository
builds microbial co-occurrence networks for 74 gut microbiome cohorts
(73 disease/control cohorts plus the ERAS antibiotic time course), then runs
downstream analyses of network stability, edge retention, hub species, rare vs
abundant taxa, and ERAS recovery.

## What it does

- Builds per-group networks with SPIEC-EASI (default), FlashWeave or SparCC.
- Supports read counts and relative abundances, prevalence filtering, edge
  weight/confidence filters, and merged-cohort input.
- Produces per-cohort network objects, node tables, centrality metrics, edge
  statistics and abundance-class plots.
- Runs cross-cohort analyses: stability vs sample size, node overlap, degree
  change, edge preservation/loss, rho comparisons, consensus hubs and the
  final multi-panel figures.
- Handles the ERAS time course (D0/D4/D8/D42/D180) with D0 relabelled as
  Control.

## Requirements

- R 4.4 (tested with R 4.4.1).
- Pipeline packages: `SpiecEasi`, `igraph`, `RColorBrewer`, `data.table`,
  `ggplot2`, `reshape2`, `patchwork`, `tidyr`, `cowplot`, `dplyr`, `purrr`,
  `stringr`, `Matrix`, `ggsci`.
- Analysis packages: `igraph`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`,
  `data.table`, `stringr`, `ggsci`, `scales`, `gridExtra`, `clinfun`,
  `RColorBrewer`, and `pheatmap` (optional; heatmaps are skipped when absent).
- FlashWeave and SparCC modes require their respective tooling; the default
  SPIEC-EASI mode only needs the R packages above.

Library locations are configurable in `config/pipeline.yaml`
(`runtime.library_paths` and `runtime.spiec_easi_library`). On a cluster, set
`LD_LIBRARY_PATH`/`R_LIBS` as needed before running.

```bash
export LD_LIBRARY_PATH=/path/to/R/lib
export R_LIBS=/path/to/R/library
Rscript run_pipeline.R --help
```

## Repository layout

```
run_pipeline.R        pipeline entry point (CLI and config handling)
run.sh                ready-to-run scenarios and output layout
R/                    network construction and per-cohort workflow
config/               default configs, cohort tables and ERAS metadata
analysis/             downstream analysis and figure entry point
verification/         developer checks for the refactored code
reproduction/         recorded recipes and code snapshots for published runs
archive/              refactor-era documentation
```

## Input data

Each cohort is described by one row in a tab-separated table such as
`config/cohorts.tsv`:

```text
disease  cohort  metadata_file  species_relative_abund  species_abund  maaslin3_res_file
```

- `cohort`: directory name used for the output, e.g. `PRJNA398089_598_CD`.
- `metadata_file`: sample metadata with `Sample_ID`, `Cohort` and
  `Standard_Disease` (or `Group`) columns.
- `species_relative_abund`: relative-abundance matrix, samples × species.
- `species_abund`: read-count matrix with the same samples.
- `maaslin3_res_file`: optional differential-marker table; `NA` means none.

A five-column layout without `species_abund` is also accepted. The ERAS cohort
uses its own `config/eras.tsv` and metadata file.

## Quick start

Build networks for all non-ERAS cohorts:

```bash
Rscript run_pipeline.R \
  --config config/pipeline.yaml \
  --cohorts config/cohorts.tsv \
  --outdir results/network
```

Build the ERAS time course into the same tree:

```bash
Rscript run_pipeline.R \
  --config config/eras.yaml \
  --outdir results/network
```

Run the downstream analysis:

```bash
Rscript analysis/run_analysis.R \
  --network-dir results/network \
  --results results/analysis \
  --suffix res_rc \
  --step all
```

Run one cohort while debugging:

```bash
Rscript run_pipeline.R \
  --config config/pipeline.yaml \
  --cohorts config/cohorts.tsv \
  --only PRJNA398089_598_CD \
  --outdir results/network
```

The equivalent one-command recipes are in `run.sh`:

```bash
./run.sh spiec_rc          # read counts + SPIEC-EASI + 10% prevalence
./run.sh eras              # ERAS time course
./run.sh full              # 73 cohorts + ERAS, then all analysis steps
./run.sh analysis          # analysis on an existing result tree
DRY_RUN=1 ./run.sh full    # print commands only
```

`run.sh` writes `<scenario>/network/` for cohort networks and
`<scenario>/analysis/` for downstream results; logs and effective configs are
saved alongside them. Use `--only`, `--outdir` and the other flags directly
with `run_pipeline.R` when you need a custom run.

## Configuration

`config/pipeline.yaml` holds the defaults for:

- input type (`read_count` / `relative_abundance`), prevalence and filters;
- network tool and SPIEC-EASI parameters (including StARS/pulsar settings);
- abundance-class thresholds, colours and labels;
- hub definition;
- output directory/suffix and runtime library paths.

Command-line flags override the config; `--set path/to/field=value` covers
nested fields:

```bash
--set network/spiec_easi/pulsar/rep.num=50
--set abundance/thresholds/high=2e-3
```

Use `--print-config` to inspect the effective config without running, and
`--save-config FILE` to write it (YAML or R format, chosen by extension) before
continuing. The ERAS run uses `config/eras.yaml` and `config/eras.tsv`.

## Output layout

```
results/
  network/
    PRJNA398089_598_CD/res_rc/
      g_Control_full_network.rds
      g_CD_full_network.rds
      Control_node_info.csv
      Control_centrality_metrics.csv
      ...
  analysis/
    stab/
    net_compare/node_overlap/
    mean_degree_change/
    net_compare2_V3/
    rho_compare/
    consensus_hub/
    scripts/final_figures/
```

The suffix (`res_rc` by default) is derived from the input/filter/tool
combination, or set explicitly with `--suffix`.

## Downstream analysis

`analysis/run_analysis.R` runs one or more steps:

```text
stability   stability vs sample size (fig_net_1)
node        node overlap (fig_net_2)
degree      mean degree change (fig_net_3)
net_compare edge preservation/loss (fig_net_4, fig_net_4m)
rho         rho comparisons (fig_net_5)
rho_edge    rho by edge type (fig_net_6)
hub         consensus hubs
hub_summary standalone hub summary plot
final       combined main figure and split panels
```

`--step all` runs every step except `hub_summary`, which needs an external
summary table (`--input`/`--outdir`). `--network-dir` and `--results` may point
at the same tree for the legacy layout, or at separate network/analysis roots
as shown above.

## Reproducing a published run

Recorded recipes and snapshots live under `reproduction/`. The 2026-11-14 full
run is documented in
[`reproduction/2026-11-14-full-run/`](reproduction/2026-11-14-full-run/);
its `reproduce.sh` also saves the executed code/config files before launching
the workers.

## License

MIT License. See [LICENSE](LICENSE).
