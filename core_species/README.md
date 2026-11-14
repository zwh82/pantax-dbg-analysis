# Core-species network analysis

Characterises the role of globally defined core species in the per-cohort
networks: edge enrichment in healthy controls, loss in disease, targeted
removal robustness, and recovery during the ERAS antibiotic time course.

The module reads existing `PR*/res_rc/g_*_full_network.rds` objects. It does
not rebuild networks and does not modify upstream results.

## Definitions

- **Core rare / core abundant**: fixed global lists in
  `data/core_rare_species.tsv` (143 species) and
  `data/core_abundant_species.tsv` (35 species).
- **Local abundance**: each network's own `abundance_category` is kept.
  `high` becomes *abundant*; `low` and `very_low` become *rare*; everything
  else is *other*.
- A species is counted as core only when it is in the global list **and** is
  locally abundant or rare. A globally core species that is locally `other` is
  treated as non-core.
- **Preserved edge**: the edge still exists in the disease network and both
  endpoint abundance classes are unchanged. `low ↔ very_low` is not a class
  change. Lost and gained edges are the complements.
- **ERAS recovery**: Control is the baseline, followed by D4, D8, D42 and
  D180. Strict recovery uses the same edge and abundance-class rule.
- **Consensus hub** (supplementary topology figure): a node in the top 10% of
  at least 3 of 5 centrality metrics (degree, betweenness, closeness,
  eigenvector, PageRank) with degree ≥ 3. The main topology panel uses the
  simpler top-10%-degree hub definition from the analysis tables.

ACA and PreHypertension enter single-network descriptions only; they are not
used as primary disease pairs. ERAS is excluded from all health/disease
summaries and appears only in the recovery panel.

## Layout

```
core_species/
  run_core_species.R        entry point (CLI, config, steps)
  config/core_species.yaml  input paths, steps, permutations, outputs
  data/                     global core species lists
  R/                        discovery, edges, permutation, roles,
                            disease comparison, recovery, figures, manuscript
  verification/             developer comparison utilities
```

## Requirements

`igraph`, `dplyr`, `tidyr`, `purrr`, `survival` for the analysis;
`ggplot2`, `patchwork`, `scales`, `ragg` and `svglite` for figures. R 4.4 is
tested. Library paths are configurable under `runtime.library_paths`.

## Quick start

Analyse an existing network tree:

```bash
Rscript core_species/run_core_species.R \
  --base-dir results/network \
  --step analysis
```

Run tables, figures and the manuscript summary in one go:

```bash
./run.sh core_species          # CORE_NETWORK=$OUT_ROOT/full/network by default
CORE_NETWORK=/path/to/network ./run.sh core_species
CORE_STEP=analysis ./run.sh core_species
```

The equivalent direct command:

```bash
Rscript core_species/run_core_species.R \
  --base-dir results/network \
  --outdir results/core_species \
  --figure-dir results/core_species/figures \
  --manuscript-dir results/core_species/manuscript_stats \
  --step all \
  --n-perm 1000 --n-attack-perm 100 --seed 123
```

Use `--n-perm 20 --n-attack-perm 10` for a quick check and `--print-config`
to inspect the resolved paths. `--save-config FILE` writes the effective
config before running. Any nested field can be set with
`--set analysis/n_perm=50`.

## Configuration

`config/core_species.yaml`:

- `network_dir`: directory holding `PR*/res_rc` networks;
- `core_rare`, `core_abundant`: core lists (relative to the repository root);
- `step`: `analysis`, `figures`, `manuscript` or `all` (comma-separated);
- `analysis.n_perm`, `analysis.n_attack_perm`, `analysis.seed`;
- `analysis.consensus_hub`: vote threshold and top fraction;
- `output.results`, `output.figures`, `output.manuscript`;
- `runtime.library_paths`.

Relative paths are resolved against the repository root, so commands can be
run from any directory.

## Outputs

```
results/core_species/
  run_manifest.tsv
  qc/
    network_manifest.tsv
    network_core_coverage.tsv
  tables/
    core_species_node_roles.tsv
    core_species_node_role_summary.tsv
    core_edge_enrichment_observed.tsv
    core_edge_enrichment_permutation.tsv
    core_species_network_summary.tsv
    core_species_targeted_removal.tsv
    core_edge_disease_loss.tsv
    core_edge_disease_loss_stats.tsv
    abundant_rare_core_rare_by_edge_status.tsv
    core_edge_weight_sign_rewiring.tsv
    core_edge_antibiotic_recovery_curve.tsv
    core_edge_first_recovery_time.tsv
    core_edge_antibiotic_recovery_stats.tsv
    abundant_rare_core_rare_fraction.tsv
  plot_data/
    plot_core_edge_enrichment.tsv
    plot_core_edge_loss_rate.tsv
    plot_core_edge_recovery_curve.tsv
    plot_abundant_rare_core_rare_fraction.tsv
  figures/
    Figure_core_species_network_role.{pdf,png}
    Figure_core_species_network_role2.{pdf,png}
    FigureS_core_species_topology.{pdf,svg,tiff,png}
    SourceData_*.tsv
  manuscript_stats/
    manuscript_statistics.md
    run_summary.tsv
    healthy_network_enrichment.tsv
    disease_edge_loss_by_category.tsv
    disease_edge_loss_contrasts.tsv
    abundant_rare_status_core_fraction.tsv
    abundant_rare_status_core_fraction_contrasts.tsv
    eras_recovery_curve.tsv
    eras_recovery_tests.tsv
    healthy_node_roles.tsv
```

`step: figures` writes the main figure, then redraws the topology supplement
with the consensus-hub definition, so the final
`FigureS_core_species_topology` uses consensus hubs. `SourceData_*.tsv` files
carry the values behind each panel.

Ratios with a zero denominator are written as `NA`. The analysis covers nodes
present in the constructed networks; isolates removed upstream are absent.

## Notes

- `eigen_centrality` uses an ARPACK start vector, so the last digits of
  eigenvector centrality in the consensus-hub source data can vary between
  runs of the unseeded figure step. Analysis tables are seeded and stable.
- The core lists are copied from the global core-species analysis so this
  module is self-contained.
