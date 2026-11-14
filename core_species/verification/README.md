# Verification notes

Developer checks used while refactoring the core-species analysis. They are not
part of the analysis workflow.

## Equivalence with the original scripts

A two-cohort fixture was built by symlinking the `g_*_full_network.rds` files
of one disease/control cohort (`PRJNA398089_598_CD`) and the ERAS time course
(`PRJEB20800_57_ERAS`) into a temporary `PR*/res_rc` tree. Both the original
`core_species_network_analysis.R` and `run_core_species.R --step all` were run
on it with `--n-perm 10 --n-attack-perm 5 --seed 123`, then compared with
`compare_core_results.R`.

Result: all 44 common tables and source-data files match. `run_manifest.tsv`
and `run_summary.tsv` differ only in recorded paths and timestamps, which is
expected because the refactor points at the repository copies of the core
lists. The graph objects are read, never rewritten.

## Core Control networks

All 20 countries were rebuilt with

```bash
Rscript run_pipeline.R --config config/core_network.yaml --outdir /tmp/core_net_all
```

and compared with the legacy `all_country_core/` tree:

- identical file lists (140 files, 7 per country);
- all `*.csv` / `*.tsv` tables match within `1e-9`;
- for every country, identical vertex and edge counts, vertex names, adjacency
  matrix and absolute edge weights in `g_Control_full_network.rds`.

`Control_centrality_metrics.csv` can differ in the last digits of the
eigenvector column for the same ARPACK reason noted below. The Italy rebuild
done during development matched the legacy file exactly.

## Known nondeterminism

`igraph::eigen_centrality` uses an ARPACK start vector, so the last digits of
the eigenvector column in the unseeded consensus-hub source data and in
`Control_centrality_metrics.csv` can differ between runs. Running the original
script twice shows the same variation. Core-species analysis tables are
produced under a fixed seed and are stable.

## Re-running a comparison

```bash
Rscript core_species/verification/compare_core_results.R REFERENCE CANDIDATE
```

The comparator matches tables by shape and column names, compares numeric
columns with tolerance `1e-9`, and skips `run_manifest.tsv` /
`run_summary.tsv` unless `--strict` is passed.
