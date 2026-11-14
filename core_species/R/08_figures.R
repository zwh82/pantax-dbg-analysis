# Publication figures for the core-species results.
#
# `run_core_species_figures` reproduces the main four-panel figure and the
# degree-hub topology supplement.  `run_core_species_consensus_figures`
# redraws the topology supplement with the consensus-hub definition and is run
# afterwards, matching the original two-script workflow.

core_palette <- function() {
  c(
    navy = "#355C7D",
    blue = "#4C78A8",
    teal = "#2A9D8F",
    gold = "#D89C2B",
    coral = "#D6655A",
    purple = "#7A6FAC",
    grey = "#9AA0A6",
    light_grey = "#D9DDE1",
    dark = "#252525"
  )
}

core_theme_paper <- function(base_size = 7, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.4, colour = "#555555"),
      plot.tag = ggplot2::element_text(size = 9, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

save_core_bundle <- function(plot, stem, width_mm = 183, height_mm = 145, dpi = 600,
                             formats = c("pdf", "png")) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  if ("pdf" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".pdf"), plot,
      width = width_in, height = height_in,
      device = cairo_pdf, bg = "white"
    )
  }
  if ("svg" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".svg"), plot,
      width = width_in, height = height_in,
      device = svglite::svglite, bg = "white"
    )
  }
  if ("tiff" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".tiff"), plot,
      width = width_in, height = height_in,
      device = ragg::agg_tiff, dpi = dpi, compression = "lzw", bg = "white"
    )
  }
  if ("png" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".png"), plot,
      width = width_in, height = height_in,
      device = ragg::agg_png, dpi = 300, bg = "white"
    )
  }
}

run_core_species_figures <- function(result_dir, figure_dir) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- core_palette()
  ggplot2::theme_set(core_theme_paper())

  label_metric <- c(
    core_associated = "All core-associated",
    core_rare_associated = "Core rare-associated",
    core_abundant_associated = "Core abundant-associated",
    core_rare_abundant = "Core rare–Abundant",
    core_rare_rare = "Core rare–Rare",
    core_abundant_abundant = "Core abundant–Abundant"
  )

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
  node_roles <- read_core_tsv(file.path(
    result_dir, "tables", "core_species_node_role_summary.tsv"
  ))
  run_info <- read_core_tsv(file.path(result_dir, "run_manifest.tsv"))

  # ERAS is the dedicated antibiotic cohort. Exclude it from all general
  # health/disease summaries and retain it only in the recovery panel.
  antibiotic_disease <- "ERAS"
  antibiotic_cohort <- "PRJEB20800_57_ERAS"

  # Panel a: healthy-network enrichment.
  enrich_control <- enrichment %>%
    dplyr::filter(
      condition == "Control",
      cohort != antibiotic_cohort,
      metric %in% names(label_metric),
      is.finite(z_score)
    ) %>%
    dplyr::mutate(
      metric_label = factor(label_metric[metric], levels = rev(unname(label_metric))),
      emphasis = ifelse(metric == "core_rare_abundant", "Core rare–abundant", "Other")
    )

  p_a <- ggplot2::ggplot(enrich_control, ggplot2::aes(z_score, metric_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2,
      colour = pal["grey"]) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = emphasis), width = 0.58, outlier.shape = NA,
      linewidth = 0.35, alpha = 0.88
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = emphasis),
      position = ggplot2::position_jitter(height = 0.12, width = 0),
      size = 0.8, alpha = 0.42
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Core rare–abundant" = unname(pal["coral"]),
      "Other" = unname(pal["light_grey"])
    )) +
    ggplot2::scale_colour_manual(values = c(
      "Core rare–abundant" = unname(pal["coral"]),
      "Other" = unname(pal["grey"])
    )) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Healthy-network edge enrichment",
      subtitle = paste0("Exploratory; ", run_info$n_perm[1], " permutations per network"),
      x = "Permutation z-score", y = NULL
    ) +
    ggplot2::theme(legend.position = "none")

  p_a2 <- ggplot2::ggplot(enrich_control, ggplot2::aes(z_score, metric_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2,
      colour = pal["grey"]) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = emphasis), width = 0.58, outlier.shape = NA,
      linewidth = 0.35, alpha = 0.88
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = emphasis),
      position = ggplot2::position_jitter(height = 0.12, width = 0),
      size = 0.8, alpha = 0.42
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Core rare–abundant" = unname(pal["coral"]),
      "Other" = unname(pal["light_grey"])
    )) +
    ggplot2::scale_colour_manual(values = c(
      "Core rare–abundant" = unname(pal["coral"]),
      "Other" = unname(pal["grey"])
    )) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Permutation z-score", y = NULL) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  # Panel b: paired within-cohort loss-rate differences.
  loss_wide <- loss %>%
    dplyr::filter(
      disease_name != antibiotic_disease,
      edge_category %in% c(
        "core_rare_associated", "noncore_rare_associated",
        "core_associated", "noncore_noncore",
        "core_rare_abundant", "other_abundant_rare"
      )
    ) %>%
    dplyr::select(disease_name, cohort, edge_category, loss_rate) %>%
    tidyr::pivot_wider(names_from = edge_category, values_from = loss_rate)

  loss_delta <- dplyr::bind_rows(
    dplyr::transmute(loss_wide,
      disease_name, cohort,
      contrast = "Core rare assoc.\n− non-core rare assoc.",
      delta = core_rare_associated - noncore_rare_associated
    ),
    dplyr::transmute(loss_wide,
      disease_name, cohort,
      contrast = "All core assoc.\n− non-core edges",
      delta = core_associated - noncore_noncore
    ),
    dplyr::transmute(loss_wide,
      disease_name, cohort,
      contrast = "Core rare–abundant\n− other A–R",
      delta = core_rare_abundant - other_abundant_rare
    )
  ) %>%
    dplyr::filter(is.finite(delta)) %>%
    dplyr::mutate(contrast = factor(contrast, levels = unique(contrast)))

  p_b <- ggplot2::ggplot(loss_delta, ggplot2::aes(contrast, delta)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2,
      colour = pal["grey"]) +
    ggplot2::geom_violin(fill = pal["light_grey"], colour = NA, width = 0.82,
      trim = TRUE) +
    ggplot2::geom_boxplot(
      width = 0.28, outlier.shape = NA, fill = "white",
      colour = pal["dark"], linewidth = 0.35
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      size = 0.75, alpha = 0.42, colour = pal["navy"]
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      title = "Disease-associated edge loss",
      subtitle = "Paired difference within each cohort",
      x = NULL, y = "Difference in loss rate"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 22, hjust = 1))

  # Panel c: rare-core fraction among abundant-rare edges by edge status.
  status_fraction <- status_fraction %>%
    dplyr::filter(
      disease_name != antibiotic_disease,
      status %in% c("gained", "lost", "preserved")
    ) %>%
    dplyr::mutate(status = factor(status,
      levels = c("gained", "lost", "preserved"),
      labels = c("Gained", "Lost", "Preserved")
    ))

  p_c <- ggplot2::ggplot(status_fraction,
    ggplot2::aes(status, fraction_core_rare, fill = status)) +
    ggplot2::geom_violin(width = 0.82, trim = TRUE, alpha = 0.72, colour = NA) +
    ggplot2::geom_boxplot(width = 0.27, outlier.shape = NA, fill = "white",
      linewidth = 0.35) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      size = 0.75, alpha = 0.38, colour = pal["dark"]
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Gained" = unname(pal["gold"]),
      "Lost" = unname(pal["coral"]),
      "Preserved" = unname(pal["teal"])
    )) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA)) +
    ggplot2::labs(
      title = "Core rare contribution to abundant–rare edges",
      subtitle = "Each point represents one cohort",
      x = NULL, y = "Core rare fraction"
    ) +
    ggplot2::theme(legend.position = "none")

  p_c2 <- ggplot2::ggplot(status_fraction,
    ggplot2::aes(status, fraction_core_rare, fill = status)) +
    ggplot2::geom_violin(width = 0.82, trim = TRUE, alpha = 0.72, colour = NA) +
    ggplot2::geom_boxplot(width = 0.27, outlier.shape = NA, fill = "white",
      linewidth = 0.35) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      size = 0.75, alpha = 0.38, colour = pal["dark"]
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Gained" = unname(pal["gold"]),
      "Lost" = unname(pal["coral"]),
      "Preserved" = unname(pal["teal"])
    )) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      labels = scales::label_number(scale = 100)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Core rare fraction (%)"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  # Panel d: ERAS antibiotic recovery. Impose the biological order explicitly.
  recovery_keep <- c(
    "core_rare_associated", "noncore_rare_associated",
    "core_rare_abundant", "other_abundant_rare"
  )
  recovery_labels <- c(
    core_rare_associated = "Core rare-associated",
    noncore_rare_associated = "Non-core rare-associated",
    core_rare_abundant = "Core rare–abundant",
    other_abundant_rare = "Other abundant–rare"
  )
  recovery_plot <- recovery %>%
    dplyr::filter(cohort == antibiotic_cohort, edge_category %in% recovery_keep) %>%
    dplyr::mutate(
      timepoint = factor(timepoint, levels = c("D4", "D8", "D42", "D180")),
      category = factor(recovery_labels[edge_category],
        levels = unname(recovery_labels)
      )
    )

  p_d <- ggplot2::ggplot(recovery_plot, ggplot2::aes(timepoint, recovery_rate,
    colour = category, group = category
  )) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.7, stroke = 0.3) +
    ggplot2::scale_colour_manual(values = c(
      "Core rare-associated" = unname(pal["coral"]),
      "Non-core rare-associated" = unname(pal["grey"]),
      "Core rare–abundant" = unname(pal["purple"]),
      "Other abundant–rare" = unname(pal["gold"])
    )) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.1),
      limits = c(0, NA)) +
    ggplot2::labs(
      title = "ERAS edge recovery",
      subtitle = "Control network used as baseline",
      x = "Time after antibiotic exposure", y = "Strict recovery rate",
      colour = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = grid::unit(8, "pt")
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  main_figure <- (p_a | p_b) / (p_c | p_d) +
    patchwork::plot_layout(guides = "keep", heights = c(1, 1.05)) +
    patchwork::plot_annotation(tag_levels = "a")

  save_core_bundle(
    main_figure,
    file.path(figure_dir, "Figure_core_species_network_role"),
    width_mm = 183, height_mm = 150
  )

  main_figure2 <- (p_a2 | p_c2) +
    patchwork::plot_layout() +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1)
    )

  save_core_bundle(
    main_figure2,
    file.path(figure_dir, "Figure_core_species_network_role2"),
    width_mm = 183, height_mm = 70
  )

  # Supplementary figure: topology of core and non-core species in Control networks.
  role_labels <- c(
    core_rare = "Core rare",
    core_abundant = "Core abundant",
    noncore_rare = "Non-core rare",
    noncore_abundant = "Non-core abundant"
  )
  metric_labels <- c(
    median_degree = "Median degree",
    median_betweenness = "Median betweenness",
    median_coreness = "Median coreness",
    hub_fraction = "Hub fraction"
  )

  role_long <- node_roles %>%
    dplyr::filter(
      condition == "Control",
      disease_name != antibiotic_disease,
      node_group %in% names(role_labels)
    ) %>%
    dplyr::select(cohort, node_group, dplyr::all_of(names(metric_labels))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(metric_labels)),
      names_to = "metric", values_to = "value"
    ) %>%
    dplyr::mutate(
      node_group = factor(role_labels[node_group], levels = unname(role_labels)),
      metric = factor(metric_labels[metric], levels = unname(metric_labels))
    )

  p_role <- ggplot2::ggplot(role_long,
    ggplot2::aes(node_group, value, fill = node_group)) +
    ggplot2::geom_violin(width = 0.85, trim = TRUE, alpha = 0.68, colour = NA) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white",
      linewidth = 0.35) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      size = 0.65, alpha = 0.30, colour = pal["dark"]
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = c(
      "Core rare" = unname(pal["coral"]),
      "Core abundant" = unname(pal["purple"]),
      "Non-core rare" = unname(pal["light_grey"]),
      "Non-core abundant" = unname(pal["gold"])
    )) +
    ggplot2::labs(
      title = "Topological roles in healthy networks",
      subtitle = "Network-level summaries; each point represents one cohort",
      x = NULL, y = NULL
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 28, hjust = 1)
    )

  save_core_bundle(
    p_role,
    file.path(figure_dir, "FigureS_core_species_topology"),
    width_mm = 183, height_mm = 70
  )

  write.table(loss_delta,
    file.path(figure_dir, "SourceData_Figure_panel_b.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(enrich_control,
    file.path(figure_dir, "SourceData_Figure_panel_a.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )

  message("Figures written to: ", normalizePath(figure_dir))
}

# ---- consensus-hub topology supplement -------------------------------------
# Consensus hub definition follows scripts/meta_consensushub_abundance_count_V3.R:
#   - compute Degree, Betweenness, Closeness, Eigenvector and PageRank
#   - each metric contributes one vote if the node is in the top 10%
#   - consensus hub = at least 3 votes and degree >= 3

core_norm01 <- function(v) {
  v <- as.numeric(v)
  if (all(is.na(v))) return(rep(0, length(v)))
  r <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  if (!is.finite(r) || r == 0) rep(0, length(v)) else (v - min(v, na.rm = TRUE)) / r
}

core_is_top_metric <- function(v, top_percent = 0.10) {
  nv <- core_norm01(v)
  cutoff <- quantile(nv, 1 - top_percent, na.rm = TRUE)
  nv >= cutoff
}

core_consensus_hub_table <- function(rds_file, top_percent = 0.10,
                                     min_degree = 3, vote_thresh = 3) {
  g <- readRDS(rds_file)
  if (!inherits(g, "igraph")) stop("Not an igraph object: ", rds_file)
  g <- igraph::as_undirected(g, mode = "collapse")

  aw <- if ("weight" %in% igraph::edge_attr_names(g)) abs(igraph::E(g)$weight) else
    rep(1, igraph::ecount(g))
  if (!length(aw)) aw <- NULL
  distance_w <- if (is.null(aw)) NULL else 1 / (aw + 1e-10)

  cent <- dplyr::tibble(
    species = as.character(igraph::V(g)$name),
    consensus_degree = igraph::degree(g),
    consensus_betweenness = igraph::betweenness(g, weights = aw),
    consensus_closeness = igraph::closeness(g, weights = distance_w),
    consensus_eigenvector = igraph::eigen_centrality(g, weights = aw)$vector,
    consensus_pagerank = igraph::page_rank(g, weights = aw)$vector
  )

  votes <- Reduce(`+`, list(
    core_is_top_metric(cent$consensus_degree, top_percent),
    core_is_top_metric(cent$consensus_betweenness, top_percent),
    core_is_top_metric(cent$consensus_closeness, top_percent),
    core_is_top_metric(cent$consensus_eigenvector, top_percent),
    core_is_top_metric(cent$consensus_pagerank, top_percent)
  ))

  cent %>%
    dplyr::mutate(
      consensus_votes = as.integer(votes),
      is_consensus_hub = consensus_votes >= vote_thresh &
        consensus_degree >= min_degree
    )
}

run_core_species_consensus_figures <- function(result_dir, figure_dir,
                                               top_percent = 0.10,
                                               min_degree = 3,
                                               vote_thresh = 3) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- core_palette()
  ggplot2::theme_set(core_theme_paper())

  antibiotic_disease <- "ERAS"

  manifest <- read_core_tsv(file.path(result_dir, "qc", "network_manifest.tsv"))
  node_roles <- read_core_tsv(file.path(
    result_dir, "tables", "core_species_node_roles.tsv"
  ))
  node_summary <- read_core_tsv(file.path(
    result_dir, "tables", "core_species_node_role_summary.tsv"
  ))

  control_manifest <- manifest %>%
    dplyr::filter(condition == "Control", disease_name != antibiotic_disease)

  if (!nrow(control_manifest)) stop("No non-ERAS Control networks found in manifest.")

  consensus_detail <- purrr::map_dfr(seq_len(nrow(control_manifest)), function(i) {
    info <- control_manifest[i, ]
    core_consensus_hub_table(info$rds_file, top_percent, min_degree, vote_thresh) %>%
      dplyr::mutate(
        cohort = info$cohort,
        disease_name = info$disease_name,
        condition = info$condition,
        .before = 1
      )
  })

  role_detail <- node_roles %>%
    dplyr::filter(condition == "Control", disease_name != antibiotic_disease) %>%
    dplyr::mutate(species = as.character(species)) %>%
    dplyr::select(cohort, disease_name, condition, species, core_type, abundance_class)

  consensus_role_detail <- role_detail %>%
    dplyr::left_join(consensus_detail,
      by = c("cohort", "disease_name", "condition", "species")) %>%
    dplyr::mutate(
      node_group = dplyr::case_when(
        core_type == "core_rare" ~ "core_rare",
        core_type == "core_abundant" ~ "core_abundant",
        core_type == "noncore_rare" ~ "noncore_rare",
        core_type == "noncore_abundant" ~ "noncore_abundant",
        TRUE ~ "other"
      )
    )

  consensus_summary <- consensus_role_detail %>%
    dplyr::group_by(cohort, disease_name, condition, node_group) %>%
    dplyr::summarise(
      n_nodes = dplyr::n(),
      n_consensus_hubs = sum(is_consensus_hub, na.rm = TRUE),
      consensus_hub_fraction = ifelse(n_nodes == 0, NA_real_,
        n_consensus_hubs / n_nodes),
      median_consensus_votes = median(consensus_votes, na.rm = TRUE),
      .groups = "drop"
    )

  role_labels <- c(
    core_rare = "Core rare",
    core_abundant = "Core abundant",
    noncore_rare = "Non-core rare",
    noncore_abundant = "Non-core abundant"
  )

  metric_labels <- c(
    median_degree = "Median degree",
    median_betweenness = "Median betweenness",
    median_coreness = "Median coreness",
    hub_fraction = "Degree hub fraction",
    consensus_hub_fraction = "Consensus hub fraction"
  )

  role_long_base <- node_summary %>%
    dplyr::filter(
      condition == "Control",
      disease_name != antibiotic_disease,
      node_group %in% names(role_labels)
    ) %>%
    dplyr::select(
      cohort, disease_name, condition, node_group,
      median_degree, median_betweenness, median_coreness, hub_fraction
    )

  role_long_consensus <- consensus_summary %>%
    dplyr::filter(node_group %in% names(role_labels)) %>%
    dplyr::select(cohort, disease_name, condition, node_group, consensus_hub_fraction)

  role_long <- role_long_base %>%
    dplyr::left_join(role_long_consensus,
      by = c("cohort", "disease_name", "condition", "node_group")) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(metric_labels)),
      names_to = "metric", values_to = "value"
    ) %>%
    dplyr::mutate(
      node_group = factor(role_labels[node_group], levels = unname(role_labels)),
      metric = factor(metric_labels[metric], levels = unname(metric_labels))
    )

  p_role <- ggplot2::ggplot(role_long,
    ggplot2::aes(node_group, value, fill = node_group)) +
    ggplot2::geom_violin(width = 0.85, trim = TRUE, alpha = 0.68, colour = NA) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white",
      linewidth = 0.35) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      size = 0.65, alpha = 0.30, colour = pal["dark"]
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = c(
      "Core rare" = unname(pal["coral"]),
      "Core abundant" = unname(pal["purple"]),
      "Non-core rare" = unname(pal["light_grey"]),
      "Non-core abundant" = unname(pal["gold"])
    )) +
    ggplot2::labs(
      title = "Topological roles in healthy networks",
      subtitle = paste0(
        "Degree hub: top 10% by degree; consensus hub: top ", top_percent * 100,
        "% in at least ", vote_thresh,
        " of 5 centrality metrics and degree ≥ ", min_degree
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 28, hjust = 1)
    )

  write.table(consensus_detail,
    file.path(figure_dir, "SourceData_FigureS_consensus_hub_node_detail.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(consensus_summary,
    file.path(figure_dir, "SourceData_FigureS_consensus_hub_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(role_long,
    file.path(figure_dir, "SourceData_FigureS_core_species_topology_consensus_hub.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )

  save_core_bundle(
    p_role,
    file.path(figure_dir, "FigureS_core_species_topology"),
    width_mm = 220, height_mm = 75,
    formats = c("pdf", "svg", "tiff", "png")
  )

  message("Consensus-hub topology figure written to: ", normalizePath(figure_dir))
}
