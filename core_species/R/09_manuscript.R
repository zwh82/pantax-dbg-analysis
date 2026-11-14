# Manuscript-ready summary tables and markdown.
#
# Reads the result TSVs only; it does not recompute networks.

core_safe_div <- function(x, y) ifelse(is.na(y) | y == 0, NA_real_, x / y)

core_median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

core_mean_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

core_min_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

core_wilcox_one_sample <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  wilcox.test(x, mu = 0, exact = FALSE)$p.value
}

core_fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

core_fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA",
    paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

core_fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 3)
}

core_get_row <- function(x, col, value) {
  out <- x[x[[col]] == value, , drop = FALSE]
  if (!nrow(out)) return(x[NA, , drop = FALSE][1, , drop = FALSE])
  out[1, , drop = FALSE]
}

run_core_species_manuscript <- function(result_dir, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  antibiotic_disease <- "ERAS"
  antibiotic_cohort <- "PRJEB20800_57_ERAS"

  n_distinct <- function(x) length(unique(x[!is.na(x)]))

  metric_labels <- c(
    core_associated = "All core-associated",
    core_rare_associated = "Core rare-associated",
    core_abundant_associated = "Core abundant-associated",
    core_rare_abundant = "Core rare-abundant",
    core_rare_rare = "Core rare-rare",
    core_abundant_abundant = "Core abundant-abundant"
  )

  edge_labels <- c(
    all_edges = "All abundance-classified edges",
    core_associated = "All core-associated",
    noncore_noncore = "Non-core edges",
    core_rare_associated = "Core rare-associated",
    noncore_rare_associated = "Non-core rare-associated",
    core_rare_abundant = "Core rare-abundant",
    other_abundant_rare = "Other abundant-rare",
    core_rare_rare = "Core rare-rare",
    noncore_rare_rare = "Non-core rare-rare",
    core_abundant_associated = "Core abundant-associated",
    noncore_abundant_associated = "Non-core abundant-associated"
  )

  node_role_labels <- c(
    core_rare = "Core rare",
    core_abundant = "Core abundant",
    noncore_rare = "Non-core rare",
    noncore_abundant = "Non-core abundant"
  )

  run_info <- read_core_tsv(file.path(result_dir, "run_manifest.tsv"))
  manifest <- read_core_tsv(file.path(result_dir, "qc", "network_manifest.tsv"))
  enrichment <- read_core_tsv(file.path(
    result_dir, "plot_data", "plot_core_edge_enrichment.tsv"
  ))
  loss <- read_core_tsv(file.path(
    result_dir, "plot_data", "plot_core_edge_loss_rate.tsv"
  ))
  status_fraction <- read_core_tsv(file.path(
    result_dir, "plot_data", "plot_abundant_rare_core_rare_fraction.tsv"
  ))
  recovery <- read_core_tsv(file.path(
    result_dir, "plot_data", "plot_core_edge_recovery_curve.tsv"
  ))
  recovery_tests <- read_core_tsv(file.path(
    result_dir, "tables", "core_edge_antibiotic_recovery_stats.tsv"
  ))
  node_roles <- read_core_tsv(file.path(
    result_dir, "tables", "core_species_node_role_summary.tsv"
  ))

  eras_manifest <- manifest[manifest$disease_name == antibiotic_disease, , drop = FALSE]
  eras_manifest <- eras_manifest[order(eras_manifest$time_order), , drop = FALSE]

  run_summary <- data.frame(
    statistic = c(
      "run_time",
      "base_dir",
      "n_core_rare_species",
      "n_core_abundant_species",
      "n_total_networks",
      "n_permutations_per_network",
      "n_targeted_removal_permutations",
      "random_seed",
      "n_non_ERAS_control_networks",
      "n_non_ERAS_disease_pairs",
      "n_ERAS_baseline_networks",
      "n_ERAS_followup_networks",
      "ERAS_time_order",
      "R_version",
      "igraph_version"
    ),
    value = c(
      as.character(run_info$run_time[1]),
      as.character(run_info$base_dir[1]),
      as.character(run_info$n_core_rare[1]),
      as.character(run_info$n_core_abundant[1]),
      as.character(run_info$n_networks[1]),
      as.character(run_info$n_perm[1]),
      as.character(run_info$n_attack_perm[1]),
      as.character(run_info$seed[1]),
      as.character(n_distinct(manifest$cohort[
        manifest$condition == "Control" & manifest$cohort != antibiotic_cohort
      ])),
      as.character(n_distinct(loss$cohort[loss$disease_name != antibiotic_disease])),
      as.character(sum(manifest$disease_name == antibiotic_disease &
        manifest$role == "antibiotic_baseline")),
      as.character(sum(manifest$disease_name == antibiotic_disease &
        manifest$role == "antibiotic_followup")),
      paste(eras_manifest$condition, collapse = " -> "),
      as.character(run_info$r_version[1]),
      as.character(run_info$igraph_version[1])
    ),
    stringsAsFactors = FALSE
  )

  enrichment$enrichment_ratio <- core_safe_div(
    enrichment$observed, enrichment$random_mean
  )
  enrichment_keep <- enrichment[
    enrichment$condition == "Control" &
      enrichment$cohort != antibiotic_cohort &
      enrichment$metric %in% names(metric_labels) &
      is.finite(enrichment$z_score),
    ,
    drop = FALSE
  ]

  healthy_enrichment <- do.call(rbind, lapply(names(metric_labels), function(m) {
    x <- enrichment_keep[enrichment_keep$metric == m, , drop = FALSE]
    data.frame(
      metric = m,
      metric_label = unname(metric_labels[m]),
      n_networks = n_distinct(x$cohort),
      median_observed_edges = core_median_safe(x$observed),
      median_random_mean = core_median_safe(x$random_mean),
      median_random_sd = core_median_safe(x$random_sd),
      median_z_score = core_median_safe(x$z_score),
      median_enrichment_ratio = core_median_safe(x$enrichment_ratio),
      prop_z_score_gt_0 = core_mean_safe(x$z_score > 0),
      prop_empirical_p_lt_0.05 = core_mean_safe(x$empirical_p < 0.05),
      min_empirical_p = core_min_safe(x$empirical_p),
      stringsAsFactors = FALSE
    )
  }))

  loss_keep <- loss[
    loss$disease_name != antibiotic_disease &
      loss$edge_category %in% names(edge_labels),
    ,
    drop = FALSE
  ]

  loss_by_category <- do.call(rbind, lapply(names(edge_labels), function(cat) {
    x <- loss_keep[loss_keep$edge_category == cat, , drop = FALSE]
    data.frame(
      edge_category = cat,
      edge_label = unname(edge_labels[cat]),
      n_cohorts = n_distinct(x$cohort),
      median_loss_rate = core_median_safe(x$loss_rate),
      median_preservation_rate = core_median_safe(x$preservation_rate),
      stringsAsFactors = FALSE
    )
  }))

  loss_wide <- reshape(
    loss[loss$disease_name != antibiotic_disease,
      c("disease_name", "cohort", "edge_category", "loss_rate")],
    idvar = c("disease_name", "cohort"),
    timevar = "edge_category",
    direction = "wide"
  )

  paired_loss_contrast <- function(data, group1, group2, contrast_label) {
    c1 <- paste0("loss_rate.", group1)
    c2 <- paste0("loss_rate.", group2)
    d <- data[[c1]] - data[[c2]]
    data.frame(
      contrast = contrast_label,
      group1 = group1,
      group2 = group2,
      n_pairs = sum(is.finite(d)),
      median_difference = core_median_safe(d),
      prop_difference_lt_0 = core_mean_safe(d < 0),
      prop_difference_gt_0 = core_mean_safe(d > 0),
      wilcoxon_p = core_wilcox_one_sample(d),
      stringsAsFactors = FALSE
    )
  }

  loss_contrasts <- rbind(
    paired_loss_contrast(
      loss_wide,
      "core_rare_associated", "noncore_rare_associated",
      "Core rare-associated minus non-core rare-associated loss rate"
    ),
    paired_loss_contrast(
      loss_wide,
      "core_associated", "noncore_noncore",
      "All core-associated minus non-core edge loss rate"
    ),
    paired_loss_contrast(
      loss_wide,
      "core_rare_abundant", "other_abundant_rare",
      "Core rare-abundant minus other abundant-rare loss rate"
    ),
    paired_loss_contrast(
      loss_wide,
      "core_abundant_associated", "noncore_abundant_associated",
      "Core abundant-associated minus non-core abundant-associated loss rate"
    )
  )
  loss_contrasts$wilcoxon_bh <- p.adjust(loss_contrasts$wilcoxon_p, method = "BH")

  status_keep <- status_fraction[
    status_fraction$disease_name != antibiotic_disease &
      status_fraction$status %in% c("gained", "lost", "preserved"),
    ,
    drop = FALSE
  ]

  status_summary <- do.call(rbind, lapply(c("gained", "lost", "preserved"), function(st) {
    x <- status_keep[status_keep$status == st, , drop = FALSE]
    data.frame(
      status = st,
      n_cohorts = n_distinct(x$cohort),
      median_n_abundant_rare_edges = core_median_safe(x$n_abundant_rare_edges),
      median_n_core_rare_abundant_edges = core_median_safe(x$n_core_rare_abundant_edges),
      median_fraction_core_rare = core_median_safe(x$fraction_core_rare),
      stringsAsFactors = FALSE
    )
  }))

  status_wide <- reshape(
    status_keep[, c("disease_name", "cohort", "status", "fraction_core_rare")],
    idvar = c("disease_name", "cohort"),
    timevar = "status",
    direction = "wide"
  )

  paired_status_contrast <- function(data, status1, status2, contrast_label) {
    c1 <- paste0("fraction_core_rare.", status1)
    c2 <- paste0("fraction_core_rare.", status2)
    d <- data[[c1]] - data[[c2]]
    data.frame(
      contrast = contrast_label,
      status1 = status1,
      status2 = status2,
      n_pairs = sum(is.finite(d)),
      median_difference = core_median_safe(d),
      prop_difference_gt_0 = core_mean_safe(d > 0),
      prop_difference_lt_0 = core_mean_safe(d < 0),
      wilcoxon_p = core_wilcox_one_sample(d),
      stringsAsFactors = FALSE
    )
  }

  status_contrasts <- rbind(
    paired_status_contrast(status_wide, "preserved", "lost",
      "Preserved minus lost core rare fraction"),
    paired_status_contrast(status_wide, "preserved", "gained",
      "Preserved minus gained core rare fraction"),
    paired_status_contrast(status_wide, "lost", "gained",
      "Lost minus gained core rare fraction")
  )
  status_contrasts$wilcoxon_bh <- p.adjust(status_contrasts$wilcoxon_p, method = "BH")

  recovery_keep <- c(
    "core_rare_associated",
    "noncore_rare_associated",
    "core_rare_abundant",
    "other_abundant_rare"
  )

  eras_recovery <- recovery[
    recovery$cohort == antibiotic_cohort & recovery$edge_category %in% recovery_keep,
    ,
    drop = FALSE
  ]
  eras_recovery$edge_label <- unname(edge_labels[eras_recovery$edge_category])
  eras_recovery <- eras_recovery[
    order(eras_recovery$edge_category, eras_recovery$time_order), ,
    drop = FALSE
  ]

  eras_recovery_tests <- recovery_tests
  eras_recovery_tests$group1_label <- unname(edge_labels[eras_recovery_tests$group1])
  eras_recovery_tests$group2_label <- unname(edge_labels[eras_recovery_tests$group2])

  node_keep <- node_roles[
    node_roles$condition == "Control" &
      node_roles$disease_name != antibiotic_disease &
      node_roles$node_group %in% names(node_role_labels),
    ,
    drop = FALSE
  ]

  healthy_node_roles <- do.call(rbind, lapply(names(node_role_labels), function(g) {
    x <- node_keep[node_keep$node_group == g, , drop = FALSE]
    data.frame(
      node_group = g,
      node_group_label = unname(node_role_labels[g]),
      n_networks = n_distinct(x$cohort),
      median_degree = core_median_safe(x$median_degree),
      median_strength = core_median_safe(x$median_strength),
      median_betweenness = core_median_safe(x$median_betweenness),
      median_eigenvector = core_median_safe(x$median_eigenvector),
      median_pagerank = core_median_safe(x$median_pagerank),
      median_coreness = core_median_safe(x$median_coreness),
      median_hub_fraction = core_median_safe(x$hub_fraction),
      median_bridge_fraction = core_median_safe(x$bridge_fraction),
      stringsAsFactors = FALSE
    )
  }))

  write_tsv(run_summary, file.path(out_dir, "run_summary.tsv"))
  write_tsv(healthy_enrichment, file.path(out_dir, "healthy_network_enrichment.tsv"))
  write_tsv(loss_by_category, file.path(out_dir, "disease_edge_loss_by_category.tsv"))
  write_tsv(loss_contrasts, file.path(out_dir, "disease_edge_loss_contrasts.tsv"))
  write_tsv(status_summary, file.path(out_dir, "abundant_rare_status_core_fraction.tsv"))
  write_tsv(status_contrasts,
    file.path(out_dir, "abundant_rare_status_core_fraction_contrasts.tsv"))
  write_tsv(eras_recovery, file.path(out_dir, "eras_recovery_curve.tsv"))
  write_tsv(eras_recovery_tests, file.path(out_dir, "eras_recovery_tests.tsv"))
  write_tsv(healthy_node_roles, file.path(out_dir, "healthy_node_roles.tsv"))

  hr <- core_get_row(healthy_enrichment, "metric", "core_rare_abundant")
  ha <- core_get_row(healthy_enrichment, "metric", "core_associated")
  l_cr <- core_get_row(loss_by_category, "edge_category", "core_rare_associated")
  l_nr <- core_get_row(loss_by_category, "edge_category", "noncore_rare_associated")
  l_cra <- core_get_row(loss_by_category, "edge_category", "core_rare_abundant")
  l_oar <- core_get_row(loss_by_category, "edge_category", "other_abundant_rare")
  lc_cr <- core_get_row(loss_contrasts, "group1", "core_rare_associated")
  lc_cra <- core_get_row(loss_contrasts, "group1", "core_rare_abundant")
  s_g <- core_get_row(status_summary, "status", "gained")
  s_l <- core_get_row(status_summary, "status", "lost")
  s_p <- core_get_row(status_summary, "status", "preserved")
  sc_pl <- core_get_row(status_contrasts, "contrast",
    "Preserved minus lost core rare fraction")
  sc_pg <- core_get_row(status_contrasts, "contrast",
    "Preserved minus gained core rare fraction")

  rec_sentence <- function(category) {
    x <- eras_recovery[eras_recovery$edge_category == category, , drop = FALSE]
    x <- x[order(x$time_order), , drop = FALSE]
    paste0(
      unname(edge_labels[category]), ": ",
      paste0(x$timepoint, "=", core_fmt_pct(x$recovery_rate, 2), collapse = ", ")
    )
  }

  md <- c(
    "# Manuscript-ready statistics",
    "",
    paste0("Generated from `", normalizePath(result_dir, mustWork = FALSE), "` on ",
      Sys.time(), "."),
    "",
    "## Analysis scope",
    "",
    paste0(
      "- Core rare species: ", run_info$n_core_rare[1],
      "; core abundant species: ", run_info$n_core_abundant[1], "."
    ),
    paste0(
      "- Networks analysed: ", run_info$n_networks[1],
      "; non-ERAS healthy control networks: ",
      run_summary$value[run_summary$statistic == "n_non_ERAS_control_networks"],
      "; non-ERAS matched disease pairs: ",
      run_summary$value[run_summary$statistic == "n_non_ERAS_disease_pairs"],
      "."
    ),
    paste0(
      "- Permutations per network: ", run_info$n_perm[1],
      "; targeted-removal permutations: ", run_info$n_attack_perm[1],
      "; seed: ", run_info$seed[1], "."
    ),
    paste0("- ERAS time order: ",
      run_summary$value[run_summary$statistic == "ERAS_time_order"], "."),
    "",
    "## Healthy-network enrichment",
    "",
    paste0(
      "- Core rare-abundant edges: median z-score = ", core_fmt_num(hr$median_z_score, 2),
      "; median enrichment ratio = ", core_fmt_num(hr$median_enrichment_ratio, 3),
      "; z-score > 0 in ", core_fmt_pct(hr$prop_z_score_gt_0, 1),
      " of healthy networks; empirical P < 0.05 in ",
      core_fmt_pct(hr$prop_empirical_p_lt_0.05, 1), " of healthy networks."
    ),
    paste0(
      "- All core-associated edges: median z-score = ", core_fmt_num(ha$median_z_score, 2),
      "; median enrichment ratio = ", core_fmt_num(ha$median_enrichment_ratio, 3), "."
    ),
    "",
    "## Disease-associated edge loss",
    "",
    paste0(
      "- Core rare-associated edges: median loss rate = ",
      core_fmt_pct(l_cr$median_loss_rate, 1),
      "; non-core rare-associated edges: median loss rate = ",
      core_fmt_pct(l_nr$median_loss_rate, 1),
      "; paired median difference = ", core_fmt_num(lc_cr$median_difference, 4),
      "; lower loss for core rare-associated edges in ",
      core_fmt_pct(lc_cr$prop_difference_lt_0, 1), " of cohorts; Wilcoxon P = ",
      core_fmt_p(lc_cr$wilcoxon_p), "."
    ),
    paste0(
      "- Core rare-abundant edges: median loss rate = ",
      core_fmt_pct(l_cra$median_loss_rate, 1),
      "; other abundant-rare edges: median loss rate = ",
      core_fmt_pct(l_oar$median_loss_rate, 1),
      "; paired median difference = ", core_fmt_num(lc_cra$median_difference, 4),
      "; lower loss for core rare-abundant edges in ",
      core_fmt_pct(lc_cra$prop_difference_lt_0, 1), " of cohorts; Wilcoxon P = ",
      core_fmt_p(lc_cra$wilcoxon_p), "."
    ),
    "",
    "## Core rare fraction among abundant-rare edge fates",
    "",
    paste0(
      "- Median core rare fraction: gained = ", core_fmt_pct(s_g$median_fraction_core_rare, 1),
      "; lost = ", core_fmt_pct(s_l$median_fraction_core_rare, 1),
      "; preserved = ", core_fmt_pct(s_p$median_fraction_core_rare, 1), "."
    ),
    paste0(
      "- Preserved minus lost: median difference = ",
      core_fmt_num(sc_pl$median_difference, 4),
      "; positive in ", core_fmt_pct(sc_pl$prop_difference_gt_0, 1),
      " of cohorts; Wilcoxon P = ", core_fmt_p(sc_pl$wilcoxon_p), "."
    ),
    paste0(
      "- Preserved minus gained: median difference = ",
      core_fmt_num(sc_pg$median_difference, 4),
      "; positive in ", core_fmt_pct(sc_pg$prop_difference_gt_0, 1),
      " of cohorts; Wilcoxon P = ", core_fmt_p(sc_pg$wilcoxon_p), "."
    ),
    "",
    "## ERAS antibiotic recovery",
    "",
    paste0("- ", rec_sentence("core_rare_associated"), "."),
    paste0("- ", rec_sentence("noncore_rare_associated"), "."),
    paste0("- ", rec_sentence("core_rare_abundant"), "."),
    paste0("- ", rec_sentence("other_abundant_rare"), "."),
    "",
    "## Output files",
    "",
    "- `run_summary.tsv`",
    "- `healthy_network_enrichment.tsv`",
    "- `disease_edge_loss_by_category.tsv`",
    "- `disease_edge_loss_contrasts.tsv`",
    "- `abundant_rare_status_core_fraction.tsv`",
    "- `abundant_rare_status_core_fraction_contrasts.tsv`",
    "- `eras_recovery_curve.tsv`",
    "- `eras_recovery_tests.tsv`",
    "- `healthy_node_roles.tsv`"
  )

  writeLines(md, file.path(out_dir, "manuscript_statistics.md"))

  message("Manuscript statistics written to: ", normalizePath(out_dir, mustWork = FALSE))
}
