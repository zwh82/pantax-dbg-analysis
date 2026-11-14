# Verification

`compare_outputs.R` compares two result trees:

```bash
Rscript verification/compare_outputs.R <old_dir> <new_dir>
```

It compares:

- `g_*_full_network.rds`: sorted vertex/edge tables and graph attributes
- `.csv`/`.tsv`: parsed tables with `all.equal`
- `.txt`: exact lines
- `fig_rds/*.rds`: ggplot labels, mapping, and `ggplot_build()$data`
- `.png`: skipped

A run with 0 failures is the acceptance criterion.

## Fixtures used

The fixtures are synthetic and were kept outside the repository. Recreate them
with the recipes below.

1. **SPIEC-EASI fixture**: 30 species x 120 samples in three groups
   (`Control`, `TestDisease`, `Severe`); relative abundances in the
   high/low/very-low bands; Poisson read counts from the relative abundances;
   a small Maaslin3-like table with risk/protective markers.
2. **GML fixture**: 20 species (8 high, 6 low, 6 very-low) and fixed graphs
   for the three groups. The graphs include Rich-Rich, Rich-Rare, and
   Rare-Rare edges, positive and negative weights, and a degree-6 hub. This
   exercises the downstream plots without depending on SPIEC-EASI selection.
3. **Merged fixture**: two cohorts (`Control`, `MergedDisease`) under
   `meta_info/<cohort>/` paths, used for the `input$merged` branch.

## Cases

Each case ran the legacy `00_main.R` chain and `run_pipeline.R` on the same
fixture, then compared the result trees.

| Case | Legacy switch | New config |
|---|---|---|
| read count + SPIEC-EASI | `is_rc=TRUE`, defaults | `input$data_type="read_count"` |
| GML input | pre-placed `*.gml` | same |
| skip prevalence | `Notprefilter=TRUE` | `input$skip_prevalence_filter=TRUE` |
| relative abundance | `is_rc=FALSE` | `input$data_type="relative_abundance"` |
| network only | `build_net_only=TRUE` | `runs$network_only=TRUE` |
| merged cohorts | `is_merge=TRUE` | `input$merged=TRUE` |
| weight filter | `is_filter=TRUE` | `network$filter_by_weight=TRUE` |
| confidence filter | `is_edge_con_filter=TRUE` | `network$filter_by_confidence=TRUE` |
| SparCC | `net_tool="sparcc"` | `network$tool="sparcc"` |
| FlashWeave | `net_tool="flashweave"` | `network$tool="flashweave"` |

The legacy merged branch could not run as written (see `REVIEW.md` section 3).
For that one case the legacy copy was patched with the same four missing
values, so the comparison validates the intended merged behavior rather than
the broken original.
