# pantax-dbg-analysis

<p align="center">
  <strong>Language:</strong> <a href="README.md">English</a> | <a href="README.zh.md">Chinese</a>
</p>

Analysis of profiling data produced by `pantax-dbg`: co-occurrence network
construction for 74 gut microbiome cohorts, downstream statistics and figures,
and network-recovery analysis across the ERAS time points. The code was
refactored from `high_rare_net4`; the original tree
`/home/work/wenhai/meta/high_rare_net4` was not modified.

## Scope

The network-construction part covers the flow reachable from `00_main.R`:

```
00_main.R  ->  run_pipeline.R
01_config.R ->  config/pipeline.yaml
02_utils.R  ->  R/01_abundance.R
03_functions_network.R -> R/02_network.R
04_functions_analysis.R -> R/03_hub_analysis.R
05_functions_plot.R -> R/04_edge_plots.R
                       R/05_centrality_plots.R
                       R/06_network_plot.R
06_process_one.R -> R/07_process_cohort.R
```

Downstream analysis and plotting (`net_compare2_V3.R`, `rho_plot.R`,
`stab_plot.R`, and related scripts) were refactored separately under
`analysis/`; see `analysis/README.md`.

## Layout

```
run_pipeline.R          entry point: argument parsing, package loading, cohort loop
run.sh                  reproducible scenario recipes (explicit flags, logs, effective config)
config/
  pipeline.yaml         defaults for all run parameters (legacy 01_config.R switches)
  eras.yaml             ERAS-specific config (points at eras.tsv, runs.eras: true)
  pipeline.R            R-format example of the same config (only read via --config)
  read_config.R         config loader (YAML by default, R format also accepted)
  cohorts.tsv           input table for the 73 non-ERAS cohorts
  eras.tsv              single-row ERAS input table (with time-point data paths)
  PRJEB20800_57_ERAS_metadata.csv  ERAS sample metadata (shipped with the repo)
R/
  00_config.R           cohort-table parsing, output paths, config validation
  01_abundance.R        abundance classes, prevalence filter, matrix merging
  02_network.R          network construction (SPIEC-EASI/FlashWeave/SparCC), topology, hubs
  03_hub_analysis.R     hub and differential-species overlap
  04_edge_plots.R       edge count/strength/sign/interaction-type plots
  05_centrality_plots.R Rich vs Rare centrality comparisons
  06_network_plot.R     per-group network plots and node tables
  07_process_cohort.R   per-cohort main workflow
analysis/               downstream analysis and plotting (entry point and checks in analysis/README.md)
REVIEW.md               refactor differences and code-humanizer scan
```

## Running

All parameters are passed on the command line; `config/pipeline.yaml` only
provides defaults. There is no need to copy and edit config blocks, and no need
to keep one config copy per mode.

```bash
export LD_LIBRARY_PATH=/home/wenhai/miniconda3/envs/r_4.3/lib
export R_LIBS=/home/wenhai/miniconda3/envs/r_4.3/lib/R/library:/home/wenhai/application/Qmatey/tools/R
Rscript run_pipeline.R --help
```

Common flags (all can be combined with `--config`; command-line values win):

| Flag                                                    | Legacy variable                                             |
| ------------------------------------------------------- | ----------------------------------------------------------- |
| `--data-type read_count\|relative_abundance`           | `is_rc`                                                   |
| `--tool spiec-easi\|flashweave\|sparcc`                 | `net_tool`                                                |
| `--prevalence 0.1` / `--present 0`                  | `PREVALENCE_THRESHOLD` / `PRESENT_THRESHOLD`            |
| `--skip-prevalence-filter TRUE`                       | `Notprefilter`                                            |
| `--filter-weight TRUE` / `--filter-confidence TRUE` | `is_filter` / `is_edge_con_filter`                      |
| `--fill-zero-counts TRUE` / `--floor 1e-3`          | `is_rc_fill` / `FilterLow`                              |
| `--merged TRUE`                                       | `is_merge`                                                |
| `--eras TRUE`                                         | `ERAS`                                                    |
| `--network-only TRUE` / `--core-network-only TRUE`  | `build_net_only` / `build_core_net_only`                |
| `--spiec-easi-ra TRUE`                                | `spieceasi_ra`                                            |
| `--suffix NAME\|auto`                                  | `outdir_last`                                             |
| `--cohorts FILE` / `--outdir DIR` / `--only A,B`  | `config_file` / `all_out_dir` / single-cohort debugging |
| `--hub-top-percent` / `--hub-min-degree`            | hub definition                                              |

Any nested field can be set with `--set`; separate levels with `/`:

```bash
--set network/spiec_easi/pulsar/rep.num=50
--set abundance/thresholds/high=2e-3
```

`--print-config` prints the merged, effective config as YAML plus the derived
output directory/suffix. `--save-config FILE` writes the same config to a file
that can be passed back via `--config`, then continues the run. The format is
chosen by extension: `.yaml`/`.yml` writes YAML, anything else writes an R list;
`--config` accepts both.

```bash
# Show the effective config without running
Rscript run_pipeline.R --tool sparcc --outdir /tmp/out --print-config

# Save the effective config and run
Rscript run_pipeline.R --tool sparcc --outdir /tmp/out --save-config /tmp/out/config.yaml
```

### run.sh (reproducible recipes)

`run.sh` turns common combinations into commented scenarios. Every command
lists its flags explicitly, and the command, log and effective config are saved
next to the results:

```bash
./run.sh spiec_rc                  # read count + SPIEC-EASI + 0.1
./run.sh eras                      # ERAS table + time-point ordering
./run.sh relative                  # relative abundance + SPIEC-EASI
./run.sh flashweave                # read count + FlashWeave
./run.sh sparcc                    # read count + SparCC
./run.sh filter_weight             # edge-weight filter
./run.sh filter_conf               # edge-confidence filter
./run.sh network_only              # Control networks only
./run.sh full                      # 73 cohorts + ERAS, then all analysis steps
./run.sh analysis                  # analysis on an existing result tree
DRY_RUN=1 ./run.sh full            # print commands only
```

Outputs go to `./runs/<scenario>/` and never overwrite the original project:

```
runs/<scenario>/
  network/            one directory per cohort (PRJ.../res_rc/...)
    configs/<scenario>.yaml   effective config for this run
  analysis/           downstream tables and figures (full/analysis scenarios)
logs/<scenario>.log   full command and output
```

Use `OUT_ROOT`, `COHORTS`, `ERAS_COHORTS`, `MERGED_TABLE`,
`ANALYSIS_RESULTS`, `ANALYSIS_NETWORK`, `ANALYSIS_SUFFIX` and
`ANALYSIS_STEP` to override paths and steps. Extra `run_pipeline.R` flags can
be appended after the scenario name:

```bash
./run.sh spiec_rc --only PRJNA398089_598_CD
```

If `--outdir` is omitted, the pipeline falls back to the legacy tree
`/home/work/wenhai/meta/high_rare_net4/` and overwrites existing results;
`run.sh` always sets an explicit output directory.

### ERAS-specific config

`config/eras.yaml` is the ERAS version of `pipeline.yaml` with two changes:

- `cohorts_file: config/eras.tsv`
- `runs.eras: true` (D0 is renamed to Control; time points keep the
  D4/D8/D42/D180 order)

`config/eras.tsv` is a single-row ERAS cohort table that points at matrices that
actually exist under `themis_abundance_matrix_0508`. The metadata file
`config/PRJEB20800_57_ERAS_metadata.csv` is shipped with the repo, so the ERAS
run no longer depends on `high_rare_net3`. Run it directly:

```bash
Rscript run_pipeline.R --config config/eras.yaml --outdir /path/to/results
```

The `output.dir` in `config/eras.yaml` is still the legacy result root; without
`--outdir` it overwrites ERAS results there. `./run.sh eras` uses this cohort
table and writes to `$OUT_ROOT/eras/network/`.

## Configuration notes

`config/pipeline.yaml` keeps every mode from the legacy entry point and removes
the large, mutually conflicting comment blocks (`config/pipeline.R` is an
R-format example of the same config, read only when passed explicitly):

- `input$data_type`: `read_count` or `relative_abundance`
- `input$merged`: merged-cohort mode (legacy `is_merge`)
- `input$skip_prevalence_filter`: legacy `Notprefilter`
- `network$tool`: `spiec-easi`, `flashweave`, `sparcc`
- `network$filter_by_weight` / `filter_by_confidence`: legacy edge filters
- `runs$network_only` / `core_network_only`: legacy `build_net_only`
- `runs$eras`: global ERAS time-ordering switch
- `output$suffix`: when `NULL`, derived from the switches using the legacy
  rules (`res_rc`, `res_raclr_no_filter`, `res_flashweave_rc`, etc.)
- `abundance$thresholds` and `abundance$labels`: abundance-class thresholds and
  legend text; update the labels when changing the thresholds

`input$prevalence` is used both for the prevalence filter and in output file
names, preserving the legacy behaviour: even with
`skip_prevalence_filter = TRUE`, file names still contain `0.1`.

## Verification

The legacy and refactored code were compared file by file on synthetic inputs:

| Scenario                                    | Files compared | Result                                       |
| ------------------------------------------- | -------------- | -------------------------------------------- |
| read count + SPIEC-EASI + 10% prevalence    | 36             | match                                        |
| fixed GML network input + multi-group plots | 40             | match                                        |
| prevalence filter skipped                   | 36             | match                                        |
| relative-abundance input                    | 40             | match                                        |
| Control networks only                       | 6              | match                                        |
| merged cohorts                              | 36             | match (legacy branch was broken; see REVIEW) |
| edge-weight filter                          | 4              | match                                        |
| edge-confidence filter                      | 4              | match                                        |
| SparCC                                      | 43             | match                                        |
| FlashWeave                                  | 41             | match                                        |

Compared content: vertex/edge attributes of network RDS files, CSV/TSV tables,
statistical text, and ggplot build data and labels. PNG bytes are not compared
because they depend on the graphics device. See `REVIEW.md` for details.

## Reproduction

The full 2026-11-14 run (73 non-ERAS cohorts + ERAS, then analysis) is recorded
in [`reproduction/2026-11-14-full-run/`](reproduction/2026-11-14-full-run/).
Its `reproduce.sh` snapshots the executed code/config into
`control/executed_files_<UTC>/` before launching the 13 workers and the
downstream analysis.

## License

MIT License. See [LICENSE](LICENSE).
