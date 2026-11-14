<p align="center">
  <strong>Language:</strong> <a href="README.md">English</a> | <a href="README.zh.md">Chinese</a>
</p>

# pantax-dbg-analysis

Analysis workflows for profiling data produced by `pantax-dbg`. The repository
builds microbial co-occurrence networks for 74 gut microbiome cohorts
(73 disease/control cohorts plus the ERAS antibiotic time course), then runs
downstream analyses of network stability, edge retention, hub species, rare vs
abundant taxa, and ERAS recovery. It also builds Control-only core-species
networks for 20 countries and analyses the role of globally defined core
species across all cohorts.

## What it does

- Builds per-group networks with SPIEC-EASI (default), FlashWeave or SparCC.
- Supports read counts and relative abundances, prevalence filtering, edge
  weight/confidence filters, and merged-cohort input.
- Produces per-cohort network objects, node tables, centrality metrics, edge
  statistics and abundance-class plots.
- Builds one Control-only network per country from the core-species matrices
  (`config/cohorts_core.tsv`).
- Runs cross-cohort analyses: stability vs sample size, node overlap, degree
  change, edge preservation/loss, rho comparisons, consensus hubs and the
  final multi-panel figures.
- Analyses core-species roles across all networks: healthy edge enrichment,
  disease edge loss, targeted-removal robustness, core rare composition of
  abundant-rare edges, and ERAS antibiotic recovery.
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
- Core-species analysis: the same `igraph`/`dplyr`/`tidyr`/`purrr` stack plus
  `survival`, and `ragg`/`svglite` for figure export.
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
core_species/         core-species Control networks and role analysis
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

Build the core-species Control networks and analyse core-species roles:

```bash
# 20 Control-only core-species networks -> results/all_country_core/
Rscript run_pipeline.R \
  --config config/core_network.yaml \
  --outdir results

# role analysis over the 74-cohort networks in results/network
Rscript core_species/run_core_species.R \
  --base-dir results/network \
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
./run.sh core_network      # 20 Control-only core-species networks
./run.sh core_species      # core-species tables, figures and summary
./run.sh hadza             # Hadza five- and three-timepoint season networks
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
continuing. The ERAS run uses `config/eras.yaml` and `config/eras.tsv`. The
core-species Control networks use `config/core_network.yaml` with
`config/cohorts_core.tsv` and write to `<outdir>/all_country_core/`; the
core-species analysis has its own config at
`core_species/config/core_species.yaml` (see `core_species/README.md`). The
Hadza season runs use `config/hadza.yaml`, `config/hadza.tsv` and the two
metadata files in `config/`.

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

## Core-species analysis

`core_species/run_core_species.R` reads the same `PR*/res_rc` network objects
and writes:

```text
results/core_species/
  qc/ tables/ plot_data/          analysis tables
  figures/                        main and supplementary figures
  manuscript_stats/               summary tables and markdown
```

It uses fixed global core lists (143 rare, 35 abundant in
`core_species/data/`), keeps each network's local abundance class, and
requires a preserved edge to retain both endpoints and their abundance
classes. `--step` accepts `analysis`, `figures`, `manuscript` or `all`.
Permutation counts, seed and consensus-hub thresholds are configurable in
`core_species/config/core_species.yaml`. See
[`core_species/README.md`](core_species/README.md) for definitions and the
full output list.

## Hadza season runs

`config/hadza.tsv` holds two runs that share the same updated abundance and
count matrices (`/snode_data/jialu/hazha/...`) but use different season
metadata:

| Cohort directory | Disease label | Groups |
|---|---|---|
| `Hadza_five_timepoints` | `2014-LW` | Control (`2014-EW`), `2014-LW`, `2013-LD`, `2014-ED`, `2014-LD` |
| `Hadza_three_timepoints` | `Dry2013` | Control (`Wet2014`), `Dry2013`, `Dry2014` |

The pipeline requires a `Control` group; in these non-disease data it is a
code alias for the first season. Metadata files are
`config/hadza_metadata_five.tsv` and `config/hadza_metadata_three.tsv`.

```bash
./run.sh hadza
# or
Rscript run_pipeline.R --config config/hadza.yaml \
  --cohorts config/hadza.tsv --outdir results/hadza
```

Outputs go to `<outdir>/Hadza_five_timepoints/res_rc/` and
`<outdir>/Hadza_three_timepoints/res_rc/` with the default `res_rc` suffix.
The canonical run writes them to
`pantax-dbg-analysis_results_20261114/hadza/`; `reproduce.sh` runs the same
launcher as its final step.

## Reproducing a published run

Recorded recipes and snapshots live under `reproduction/`. The 2026-11-14 full
run is documented in
[`reproduction/2026-11-14-full-run/`](reproduction/2026-11-14-full-run/);
its `reproduce.sh` also saves the executed code/config files before launching
the workers.

## License

MIT License. See [LICENSE](LICENSE).
