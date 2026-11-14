# 2026-11-14 full run

This directory records how the 73 non-ERAS cohorts, the ERAS time course, the
20 core-species Control networks and the core-species role analysis were
rebuilt. The canonical result tree is:

```text
/home/work/wenhai/meta/pantax-dbg-analysis_results_20261114/
  network/       cohort networks
  analysis/      downstream tables and figures
  core_network/  Control-only core-species networks (20 countries)
  core_species/  core-species tables, figures and manuscript summary
  hadza/         Hadza five- and three-timepoint season networks
  control/       shard tables, launchers, logs, effective configs,
                 executed_files_<UTC>/ code and config snapshot
```

## Reproduce

```bash
bash reproduction/2026-11-14-full-run/reproduce.sh
```

The script:

1. copies the exact code, configs, cohort shards and launchers to
   `<run>/control/executed_files_<UTC>/`, together with `GIT_COMMIT`,
   `GIT_STATUS` and `COMMANDS.txt`;
2. launches 12 non-ERAS shards and the ERAS shard on node001/node003;
3. waits for all 13 workers and checks their exit codes;
4. runs all analysis steps on node001;
5. builds the 20 core-species Control networks with
   `config/core_network.yaml` and runs `core_species/run_core_species.R`
   (`--step all`, 1000 label permutations, seed 123) on node001;
6. runs the Hadza five- and three-timepoint season networks with
   `config/hadza.yaml` and `config/hadza.tsv`.

Useful overrides: `NODE_A`, `NODE_B`, `RUN`, `REPO`, `SKIP_ANALYSIS=1`,
`SKIP_CORE=1`, `SKIP_HADZA=1`, `DRY_RUN=1`, `SNAPSHOT_ONLY=1`.

## Launchers and `run_recorded`

The script does not rebuild the `Rscript` command for every step. Each step is
saved as a launcher under `<run>/control/launchers/`:

```text
control/launchers/
  w01_node001.sh ... w12_node003.sh   cohort shards
  e13_eras.sh                         ERAS shard
  analysis.sh                         downstream analysis
  core_network.sh                     20 core-species Control networks
  core_species.sh                     core-species tables/figures/summary
  hadza.sh                            Hadza season networks
```

A launcher contains the exact `cd`, environment variables, input/output paths
and `Rscript` arguments used for that step. It writes its status to
`control/status/<name>.state` and its exit code to
`control/status/<name>.exit`; stdout/stderr go to
`control/logs/<name>.log`.

`run_recorded <name> <timeout-seconds>` launches
`control/launchers/<name>.sh` on `NODE_A`, records the PID, polls for
`<name>.exit` every 20 s, and fails the whole run if the step exits non-zero or
exceeds the timeout.

| Call | Launcher | Timeout |
|---|---|---:|
| `run_recorded core_network 7200` | `core_network.sh` | 2 h |
| `run_recorded core_species 21600` | `core_species.sh` | 6 h |
| `run_recorded hadza 21600` | `hadza.sh` | 6 h |

The 13 cohort shards use the same launcher pattern through `launch_worker`;
`analysis.sh` is launched by the older `run_analysis` helper, which waits the
same way.

For this run the snapshot is:

```text
/home/work/wenhai/meta/pantax-dbg-analysis_results_20261114/control/executed_files_20261114T094014Z
```

It records git commit `2f8de29234713883a01040b2940f14751166dfa8` and contains
83 files: the pipeline code, analysis code, configs, cohort shard tables,
launchers, ERAS metadata and documentation.

The data matrices under `/snode_data/` and the node names are specific to this
cluster; edit `config/pipeline.yaml`, `config/eras.tsv` and `control/configs/`
when running elsewhere.
