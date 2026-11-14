# Core-species analysis orchestration.
#
# Rewrites the same result files as the original monolithic script.  The
# random seed is set once and the network loop keeps the original order:
# per-network permutation, then targeted removal of core-rare and
# core-abundant species.

run_core_species_analysis <- function(base_dir, core_rare, core_abundant, outdir,
                                      n_perm, n_attack_perm, seed,
                                      core_rare_file, core_abundant_file) {
  if (!dir.exists(base_dir)) stop("Base directory not found: ", base_dir)
  set.seed(seed)

  prepare_output_dirs(outdir)
  manifest <- discover_networks(base_dir)
  if (!nrow(manifest)) stop("No res_rc network RDS files were found.")
  write_tsv(manifest, file.path(outdir, "qc", "network_manifest.tsv"))

  network_objects <- list()
  node_detail <- list()
  node_summary <- list()
  observed_all <- list()
  permutation_all <- list()
  network_summary <- list()
  attack_all <- list()
  qc_all <- list()

  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]
    g <- readRDS(info$rds_file)
    if (!inherits(g, "igraph")) stop("Not an igraph object: ", info$rds_file)
    if (igraph::is_directed(g)) stop("Directed graph is not supported: ", info$rds_file)
    nodes <- make_node_meta(g, core_rare, core_abundant)
    edges <- make_edge_df(g, nodes)
    key <- paste(info$cohort, info$condition, sep = "::")
    network_objects[[key]] <- list(g = g, nodes = nodes, edges = edges, info = info)

    roles <- node_role_table(g, nodes)
    node_detail[[key]] <- dplyr::bind_cols(info[rep(1, nrow(roles)),
      c("cohort", "disease_name", "condition", "role", "time_order")
    ], roles)
    node_summary[[key]] <- dplyr::bind_cols(info[rep(1, 1),
      c("cohort", "disease_name", "condition", "role", "time_order")
    ], summarise_node_roles(roles))

    observed <- edge_count_table(edges)
    counts <- setNames(observed$observed, observed$metric)
    observed_all[[key]] <- dplyr::bind_cols(info[rep(1, nrow(observed)),
      c("cohort", "disease_name", "condition", "role", "time_order")
    ], observed) %>%
      dplyr::mutate(
        total_edges = counts["all_edges"],
        observed_fraction = safe_fraction(observed, total_edges)
      )

    perm <- permutation_enrichment(g, nodes, observed, n_perm)
    permutation_all[[key]] <- dplyr::bind_cols(info[rep(1, nrow(perm)),
      c("cohort", "disease_name", "condition", "role", "time_order")
    ], perm)

    network_summary[[key]] <- dplyr::tibble(
      cohort = info$cohort,
      disease_name = info$disease_name,
      condition = info$condition,
      role = info$role,
      time_order = info$time_order,
      n_nodes = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      n_abundance_edges = counts["all_edges"],
      n_core_rare_species = sum(nodes$core_type == "core_rare"),
      n_core_abundant_species = sum(nodes$core_type == "core_abundant"),
      n_core_associated_edges = counts["core_associated"],
      fraction_core_associated_edges =
        safe_fraction(counts["core_associated"], counts["all_edges"]),
      n_core_rare_associated_edges = counts["core_rare_associated"],
      fraction_core_rare_associated_edges =
        safe_fraction(counts["core_rare_associated"], counts["all_edges"]),
      n_core_abundant_associated_edges = counts["core_abundant_associated"],
      fraction_core_abundant_associated_edges =
        safe_fraction(counts["core_abundant_associated"], counts["all_edges"]),
      n_abundant_rare_edges = counts["abundant_rare"],
      n_core_rare_abundant_edges = counts["core_rare_abundant"],
      fraction_core_rare_among_abundant_rare =
        safe_fraction(counts["core_rare_abundant"], counts["abundant_rare"])
    )
    attack_all[[paste0(key, "::rare")]] <- dplyr::bind_cols(
      info[rep(1, 1), c("cohort", "condition")],
      attack_summary(g, nodes, "core_rare", n_attack_perm)
    )
    attack_all[[paste0(key, "::abundant")]] <- dplyr::bind_cols(
      info[rep(1, 1), c("cohort", "condition")],
      attack_summary(g, nodes, "core_abundant", n_attack_perm)
    )
    qc_all[[key]] <- dplyr::tibble(
      cohort = info$cohort,
      condition = info$condition,
      n_nodes = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      core_rare_coverage =
        safe_fraction(sum(nodes$core_type == "core_rare"), length(core_rare)),
      core_abundant_coverage =
        safe_fraction(sum(nodes$core_type == "core_abundant"), length(core_abundant)),
      core_rare_locally_abundant =
        sum(nodes$core_type == "core_rare" & nodes$abundance_class == "abundant"),
      core_abundant_locally_rare =
        sum(nodes$core_type == "core_abundant" & nodes$abundance_class == "rare")
    )
  }

  node_detail_df <- bind_rows_safe(node_detail)
  node_summary_df <- bind_rows_safe(node_summary)
  observed_df <- bind_rows_safe(observed_all)
  permutation_df <- bind_rows_safe(permutation_all)
  summary_df <- bind_rows_safe(network_summary)
  attack_df <- bind_rows_safe(attack_all)
  qc_df <- bind_rows_safe(qc_all)

  write_tsv(node_detail_df, file.path(outdir, "tables", "core_species_node_roles.tsv"))
  write_tsv(node_summary_df,
    file.path(outdir, "tables", "core_species_node_role_summary.tsv"))
  write_tsv(observed_df,
    file.path(outdir, "tables", "core_edge_enrichment_observed.tsv"))
  write_tsv(permutation_df,
    file.path(outdir, "tables", "core_edge_enrichment_permutation.tsv"))
  write_tsv(summary_df,
    file.path(outdir, "tables", "core_species_network_summary.tsv"))
  write_tsv(attack_df,
    file.path(outdir, "tables", "core_species_targeted_removal.tsv"))
  write_tsv(qc_df, file.path(outdir, "qc", "network_core_coverage.tsv"))

  disease_loss <- list()
  disease_status_fraction <- list()
  rewiring <- list()
  primary <- manifest %>% dplyr::filter(role %in% c("control", "disease"))
  for (cohort in unique(primary$cohort)) {
    rows <- primary[primary$cohort == cohort, ]
    if (!all(c("control", "disease") %in% rows$role)) next
    ci <- rows[rows$role == "control", ][1, ]
    di <- rows[rows$role == "disease", ][1, ]
    control <- network_objects[[paste(ci$cohort, ci$condition, sep = "::")]]$edges
    disease <- network_objects[[paste(di$cohort, di$condition, sep = "::")]]$edges
    compared <- compare_pair(control, disease)
    loss <- loss_summary(compared$control)
    disease_loss[[cohort]] <- dplyr::mutate(loss,
      disease_name = di$disease_name, cohort = cohort, .before = 1
    )
    statuses <- dplyr::bind_rows(compared$control, compared$gained)
    disease_status_fraction[[cohort]] <- status_core_fraction(statuses) %>%
      dplyr::mutate(disease_name = di$disease_name, cohort = cohort, .before = 1)
    common <- dplyr::inner_join(
      dplyr::select(control, edge_id, control_weight = weight, control_sign = sign),
      dplyr::select(disease, edge_id, disease_weight = weight, disease_sign = sign),
      by = "edge_id"
    )
    rewiring[[cohort]] <- common %>%
      dplyr::mutate(
        disease_name = di$disease_name,
        cohort = cohort,
        weight_change = disease_weight - control_weight,
        sign_flip = control_sign != disease_sign
      )
  }

  disease_loss_df <- bind_rows_safe(disease_loss)
  status_fraction_df <- bind_rows_safe(disease_status_fraction)
  rewiring_df <- bind_rows_safe(rewiring)
  loss_stats_df <- if (nrow(disease_loss_df)) {
    paired_loss_tests(disease_loss_df)
  } else {
    data.frame()
  }

  write_tsv(disease_loss_df, file.path(outdir, "tables", "core_edge_disease_loss.tsv"))
  write_tsv(loss_stats_df,
    file.path(outdir, "tables", "core_edge_disease_loss_stats.tsv"))
  write_tsv(status_fraction_df,
    file.path(outdir, "tables", "abundant_rare_core_rare_by_edge_status.tsv"))
  write_tsv(rewiring_df,
    file.path(outdir, "tables", "core_edge_weight_sign_rewiring.tsv"))

  eras_rows <- manifest %>%
    dplyr::filter(disease_name == "ERAS",
      role %in% c("antibiotic_baseline", "antibiotic_followup")) %>%
    dplyr::arrange(time_order)
  recovery_curve_df <- data.frame()
  first_recovery_df <- data.frame()
  recovery_stats_df <- data.frame()
  if (nrow(eras_rows)) {
    baseline_row <- eras_rows[eras_rows$role == "antibiotic_baseline", ][1, ]
    follow_rows <- eras_rows[eras_rows$role == "antibiotic_followup", ]
    baseline <- network_objects[[paste(baseline_row$cohort,
      baseline_row$condition, sep = "::"
    )]]$edges
    followups <- setNames(purrr::map(seq_len(nrow(follow_rows)), function(i) {
      network_objects[[paste(follow_rows$cohort[i],
        follow_rows$condition[i], sep = "::"
      )]]$edges
    }), follow_rows$condition)
    rec <- recovery_tables(baseline, followups, follow_rows$time_order)
    recovery_curve_df <- dplyr::mutate(rec$curve, cohort = baseline_row$cohort, .before = 1)
    first_recovery_df <- dplyr::mutate(rec$first, cohort = baseline_row$cohort, .before = 1)

    masks <- edge_masks(baseline)
    comparisons <- list(
      c("core_rare_associated", "noncore_rare_associated"),
      c("core_associated", "noncore_noncore"),
      c("core_rare_abundant", "other_abundant_rare"),
      c("core_rare_rare", "noncore_rare_rare")
    )
    recovery_stats_df <- bind_rows_safe(purrr::map(comparisons, function(pair) {
      group <- rep(NA_character_, nrow(baseline))
      group[masks[[pair[1]]]] <- pair[1]
      group[masks[[pair[2]]]] <- pair[2]
      keep <- !is.na(group)
      dat <- first_recovery_df[keep, ]
      dat$group <- group[keep]
      p_w <- if (sum(dat$is_recovered) > 1 &&
        length(unique(dat$group[dat$is_recovered])) == 2) {
        wilcox.test(first_recovery_order ~ group,
          data = dat[dat$is_recovered, ], exact = FALSE
        )$p.value
      } else {
        NA_real_
      }
      p_lr <- if (nrow(dat) > 1 && length(unique(dat$group)) == 2) {
        1 - pchisq(survival::survdiff(
          survival::Surv(first_recovery_order, is_recovered) ~ group, data = dat
        )$chisq, df = 1)
      } else {
        NA_real_
      }
      dplyr::tibble(
        cohort = baseline_row$cohort,
        group1 = pair[1],
        group2 = pair[2],
        n_edges = nrow(dat),
        wilcoxon_recovered_only_p = p_w,
        logrank_p = p_lr
      )
    })) %>%
      dplyr::mutate(
        wilcoxon_bh = p.adjust(wilcoxon_recovered_only_p, "BH"),
        logrank_bh = p.adjust(logrank_p, "BH")
      )
  }

  write_tsv(recovery_curve_df,
    file.path(outdir, "tables", "core_edge_antibiotic_recovery_curve.tsv"))
  write_tsv(first_recovery_df,
    file.path(outdir, "tables", "core_edge_first_recovery_time.tsv"))
  write_tsv(recovery_stats_df,
    file.path(outdir, "tables", "core_edge_antibiotic_recovery_stats.tsv"))

  abundant_fraction <- summary_df %>%
    dplyr::transmute(
      cohort, condition,
      timepoint = ifelse(role == "antibiotic_followup", condition, NA),
      n_abundant_rare_edges,
      n_core_rare_abundant_edges,
      n_noncore_rare_abundant_edges =
        n_abundant_rare_edges - n_core_rare_abundant_edges,
      fraction_core_rare_among_abundant_rare
    )
  write_tsv(abundant_fraction,
    file.path(outdir, "tables", "abundant_rare_core_rare_fraction.tsv"))

  write_tsv(dplyr::select(permutation_df,
    cohort, condition, time_order, metric, observed,
    random_mean, random_sd, z_score, empirical_p
  ), file.path(outdir, "plot_data", "plot_core_edge_enrichment.tsv"))
  write_tsv(dplyr::select(disease_loss_df,
    disease_name, cohort, edge_category, loss_rate, preservation_rate
  ), file.path(outdir, "plot_data", "plot_core_edge_loss_rate.tsv"))
  write_tsv(dplyr::select(recovery_curve_df,
    cohort, timepoint, time_order, edge_category, recovery_rate
  ), file.path(outdir, "plot_data", "plot_core_edge_recovery_curve.tsv"))
  write_tsv(status_fraction_df,
    file.path(outdir, "plot_data", "plot_abundant_rare_core_rare_fraction.tsv"))

  run_manifest <- dplyr::tibble(
    run_time = as.character(Sys.time()),
    base_dir = normalizePath(base_dir),
    core_rare_file = core_rare_file,
    core_abundant_file = core_abundant_file,
    n_core_rare = length(core_rare),
    n_core_abundant = length(core_abundant),
    n_networks = nrow(manifest),
    n_perm = n_perm,
    n_attack_perm = n_attack_perm,
    seed = seed,
    r_version = R.version.string,
    igraph_version = as.character(packageVersion("igraph"))
  )
  write_tsv(run_manifest, file.path(outdir, "run_manifest.tsv"))
  invisible(manifest)
}
