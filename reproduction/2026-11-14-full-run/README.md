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
   (`--step all`, 1000 label permutations, seed 123) on node001.

Useful overrides: `NODE_A`, `NODE_B`, `RUN`, `REPO`, `SKIP_ANALYSIS=1`,
`SKIP_CORE=1`, `DRY_RUN=1`, `SNAPSHOT_ONLY=1`.

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
