# Final composite main figure (fig_net_main_2) and its split panels.
#
# Inputs are produced by earlier steps:
#   <net_compare>/combined/<suffix>_main_data.rds
#   <hub>/Meta_ConsensusHub_Abundance_plot_data.rds
#   <eras_cohort>/<suffix>/data/combined_data.rds
#   <eras_cohort>/<suffix>/data/edge_count_summary.rds
#   <net_compare>/ERAS/<suffix>/Preserved_zhexian.rds

theme_pub <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.7, colour = "black"),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.6),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      plot.title = element_text(size = base_size + 1.2, face = "bold", hjust = 0.5),
      plot.margin = margin(2, 3, 2, 3),
      panel.grid = element_blank()
    )
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  as.character(signif(p, 2))
}

make_loss_overview_panel <- function(data_list) {
  edge_cols <- c(
    "R-R" = "#3F8C5A",
    "A-R" = "#D09A38",
    "A-A" = "#7E6AA2"
  )

  plot_data <- bind_rows(lapply(data_list, function(bundle) {
    bundle$plot_data %>%
      mutate(
        Disease = bundle$disease,
        friedman_p = bundle$friedman_p,
        jonck_p = bundle$jonck_p
      )
  })) %>%
    mutate(
      Disease = factor(Disease, levels = names(data_list)),
      EdgeType = factor(
        EdgeType,
        levels = c("Rare-Rare", "Abundant-Rare", "Abundant-Abundant"),
        labels = c("R-R", "A-R", "A-A")
      )
    )

  stat_df <- plot_data %>%
    distinct(Disease, friedman_p, jonck_p) %>%
    mutate(
      label = paste0(
        "Friedman p=", vapply(friedman_p, format_p, character(1)),
        "\nJT p=", vapply(jonck_p, format_p, character(1))
      ),
      x = 0.65,
      y = 0.705
    )

  ggplot(plot_data, aes(x = EdgeType, y = Rate)) +
    geom_line(
      aes(group = Cohort),
      colour = "#A9A9A9",
      alpha = 0.58,
      linewidth = 0.28
    ) +
    geom_point(
      aes(colour = EdgeType),
      fill = "white",
      shape = 21,
      size = 0.9,
      stroke = 0.25,
      alpha = 0.72,
      position = position_jitter(width = 0.08, height = 0, seed = 1)
    ) +
    geom_boxplot(
      aes(colour = EdgeType),
      fill = NA,
      width = 0.50,
      linewidth = 0.55,
      outlier.shape = NA
    ) +
    geom_text(
      data = stat_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0,
      size = 1.75,
      lineheight = 0.92,
      fontface = "italic",
      colour = "#333333"
    ) +
    facet_wrap(~ Disease, nrow = 1) +
    coord_cartesian(ylim = c(0.68, 1.01), clip = "off") +
    scale_colour_manual(values = edge_cols, guide = "none") +
    labs(x = "Interaction type", y = "Loss rate") +
    theme_pub(base_size = 6.2) +
    theme(
      panel.spacing.x = unit(4, "pt"),
      strip.text = element_text(size = 6.8, face = "bold"),
      strip.background = element_rect(fill = "#ECECEC", colour = NA),
      axis.text.x = element_text(size = 5.7, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 2)),
      axis.title.y = element_text(margin = margin(r = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
}

make_final_hub_positions <- function(master_df, disease_seq, step_within, step_between) {
  pos_list <- list()
  centers <- numeric(length(disease_seq))
  names(centers) <- disease_seq
  x_cursor <- 1

  for (disease in disease_seq) {
    cohort_order <- master_df %>%
      filter(Disease == disease) %>%
      group_by(Cohort) %>%
      summarise(Total_Hub = max(Total_Hub, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Total_Hub), Cohort) %>%
      pull(Cohort)

    if (length(cohort_order) == 0) next

    xpos <- x_cursor + (seq_along(cohort_order) - 1) * step_within
    pos_list[[disease]] <- data.frame(
      Cohort = cohort_order,
      xpos = xpos,
      Disease = disease,
      stringsAsFactors = FALSE
    )
    centers[disease] <- mean(xpos)
    x_cursor <- max(xpos) + step_between
  }

  list(pos = bind_rows(pos_list), centers = centers[centers > 0])
}

make_hub_panel <- function(hub_data) {
  hub_cols <- c(
    "rare" = "#379fb4",
    "abundant" = "#fcae59"
  )

  master_df <- hub_data$master_df
  diseases_use <- hub_data$multi_dis
  disease_seq <- hub_data$main_seq
  pinfo <- make_final_hub_positions(
    master_df = master_df,
    disease_seq = disease_seq,
    step_within = hub_data$step_within,
    step_between = hub_data$step_between
  )

  posmap <- pinfo$pos
  centers <- pinfo$centers
  df <- master_df %>%
    filter(Disease %in% diseases_use) %>%
    left_join(posmap %>% select(Cohort, Disease, xpos), by = c("Cohort", "Disease")) %>%
    mutate(
      Abundance = factor(as.character(Abundance), levels = c("rare", "abundant")),
      Group = factor(as.character(Group), levels = c("Control", "Disease")),
      prop = ifelse(Total_Hub > 0, n / Total_Hub, NA_real_)
    )

  ggplot(df, aes(x = xpos, y = prop, fill = Abundance)) +
    geom_col(width = hub_data$bar_w, linewidth = 0.12, colour = "white") +
    geom_text(
      aes(label = ifelse(n > 0, n, "")),
      position = position_stack(vjust = 0.5),
      size = 1.75,
      colour = "white",
      fontface = "bold"
    ) +
    facet_wrap(~ Group, ncol = 1) +
    scale_x_continuous(
      breaks = unname(centers),
      labels = names(centers),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) paste0(x * 100, "%"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_fill_manual(
      values = hub_cols,
      drop = FALSE,
      guide = guide_legend(
        title.position = "top",
        keywidth = unit(7, "pt"),
        keyheight = unit(5, "pt"),
        override.aes = list(linewidth = 0.15)
      )
    ) +
    labs(x = "Disease", y = "Consensus hub proportion", fill = "Abundance stratification") +
    theme_pub(base_size = 7.4) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.title = element_text(size = 6.5, face = "bold"),
      legend.text = element_text(size = 6.1),
      legend.key.width = unit(8, "pt"),
      legend.key.height = unit(5.5, "pt"),
      axis.text.x = element_text(size = 6.8, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 6.8, face = "bold"),
      axis.title.x = element_text(size = 8.0, face = "bold", margin = margin(t = 3)),
      axis.title.y = element_text(size = 8.0, face = "bold", margin = margin(r = 3)),
      strip.text = element_text(size = 7.4, face = "bold"),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA),
      panel.spacing.y = unit(2.5, "pt")
    )
}

make_eras_centrality_panel <- function(p1_data) {
  plot_data <- p1_data %>%
    mutate(
      AbundanceType = rename_rich_to_abundant(AbundanceType),
      Group = str_replace_all(as.character(Group), "Control", "D0"),
      Group = factor(Group, levels = c("D0", "D4", "D8", "D42", "D180"), ordered = TRUE),
      AbundanceType = factor(AbundanceType, levels = c("Rare", "Abundant"))
    )

  ggplot(plot_data, aes(x = AbundanceType, y = Degree, fill = Group)) +
    geom_boxplot(
      alpha = 0.85, outlier.size = 0.35, linewidth = 0.22, width = 0.68,
      position = position_dodge(0.78)
    ) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
    labs(x = "Abundance type", y = "Degree", fill = "Time point") +
    theme_pub(base_size = 6.4) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold")
    )
}

make_eras_edge_panel <- function(p2_data) {
  plot_data <- p2_data %>%
    mutate(
      Edge_Type = rename_rich_to_abundant(Edge_Type),
      Edge_Type = factor(
        Edge_Type,
        levels = c("All", "Rare-Rare", "Abundant-Rare", "Abundant-Abundant"),
        labels = c("All", "R-R", "A-R", "A-A"),
        ordered = TRUE
      ),
      Group = str_replace_all(as.character(Group), "Control", "D0"),
      Group = factor(Group, levels = c("D0", "D4", "D8", "D42", "D180"), ordered = TRUE)
    )

  ggplot(plot_data, aes(x = Edge_Type, y = Count, fill = Group)) +
    geom_col(
      position = position_dodge(0.78), width = 0.68,
      linewidth = 0.12, colour = "white"
    ) +
    geom_text(
      aes(label = Count),
      position = position_dodge(0.78),
      vjust = -0.35,
      size = 1.3
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
    labs(x = "Interaction type", y = "Edge count", fill = "Time point") +
    theme_pub(base_size = 6.4) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold")
    )
}

make_eras_preserved_panel <- function(p3_data) {
  plot_data <- p3_data %>%
    mutate(
      EdgeType = factor(
        as.character(EdgeType),
        levels = c("All", "Rare-Rare", "Abundant-Rare", "Abundant-Abundant"),
        labels = c("All", "R-R", "A-R", "A-A"),
        ordered = TRUE
      )
    )

  ggplot(plot_data, aes(x = timepoint, y = Preserved, colour = EdgeType, group = EdgeType)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.5) +
    geom_text(
      aes(label = Preserved),
      vjust = -0.6,
      size = 1.9,
      show.legend = FALSE,
      position = position_jitter(width = 0.08, height = 0)
    ) +
    coord_cartesian(
      ylim = c(
        min(plot_data$Preserved, na.rm = TRUE),
        max(plot_data$Preserved, na.rm = TRUE) * 1.13
      )
    ) +
    labs(x = "Time point", y = "Preserved edge", colour = "Interaction type") +
    theme_pub(base_size = 6.4) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold")
    )
}

save_split_panel <- function(plot, filename, width, height, split_out_dir) {
  ggsave(
    filename = file.path(split_out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "mm",
    device = cairo_pdf,
    limitsize = FALSE,
    bg = "white"
  )
  message("Saved split panel: ", file.path(split_out_dir, filename))
}

run_final_figure <- function(cfg) {
  network_suffix <- analysis_suffixes(cfg)[1]
  out_dir <- ensure_dir(file.path(cfg$results_dir, cfg$paths$final))
  split_out_dir <- ensure_dir(file.path(cfg$results_dir, cfg$paths$final_split))

  net_main_rds <- file.path(
    cfg$results_dir, cfg$paths$net_compare, "combined",
    paste0(network_suffix, "_main_data.rds")
  )
  hub_rds <- file.path(
    cfg$results_dir, cfg$paths$hub, "Meta_ConsensusHub_Abundance_plot_data.rds"
  )
  eras_centrality_rds <- file.path(
    network_dir(cfg), cfg$eras_cohort, network_suffix, "data", "combined_data.rds"
  )
  eras_edge_rds <- file.path(
    network_dir(cfg), cfg$eras_cohort, network_suffix, "data", "edge_count_summary.rds"
  )
  eras_preserved_rds <- file.path(
    cfg$results_dir, cfg$paths$net_compare, "ERAS", network_suffix,
    "Preserved_zhexian.rds"
  )

  old_theme <- theme_set(theme_pub())
  on.exit(theme_set(old_theme), add = TRUE)

  net_main <- readRDS(net_main_rds)
  hub_data <- readRDS(hub_rds)
  eras_centrality <- readRDS(eras_centrality_rds)
  eras_edge <- readRDS(eras_edge_rds)
  eras_preserved <- readRDS(eras_preserved_rds)

  if (length(net_main) != 6) {
    warning("Expected 6 network comparison panels, found ", length(net_main),
            ". The layout will still be generated.")
  }

  net_panel <- make_loss_overview_panel(net_main)
  hub_panel <- make_hub_panel(hub_data)

  eras_panels <- list(
    make_eras_centrality_panel(eras_centrality),
    make_eras_edge_panel(eras_edge),
    make_eras_preserved_panel(eras_preserved)
  )
  names(eras_panels) <- c("C", "D", "E")

  split_eras_legend_theme <- theme(
    text = element_text(family = "Arial"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(family = "Arial", size = 5.7, face = "bold"),
    legend.text = element_text(family = "Arial", size = 5.3),
    legend.key.width = unit(4.8, "pt"),
    legend.key.height = unit(4.8, "pt"),
    legend.spacing.x = unit(0.8, "pt"),
    legend.box.spacing = unit(0, "pt"),
    axis.title = element_text(family = "Arial", size = 7.2, face = "bold"),
    axis.text = element_text(family = "Arial", size = 6.4, colour = "black")
  )

  eras_split_d <- eras_panels[["D"]]
  eras_split_d$layers[[1]]$geom_params$width <- 0.78
  eras_split_d$layers[[2]]$aes_params$family <- "Arial"
  eras_split_d$layers[[2]]$aes_params$size <- 1.65

  eras_split_panels <- list(
    C = eras_panels[["C"]] +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
      split_eras_legend_theme,
    D = eras_split_d +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
      split_eras_legend_theme,
    E = eras_panels[["E"]] +
      guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
      split_eras_legend_theme
  )
  eras_split_panels[["E"]]$layers[[3]]$aes_params$family <- "Arial"

  net_split_panel <- net_panel +
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(family = "Arial", size = 7.2, face = "bold"),
      axis.text = element_text(family = "Arial", size = 6.3, colour = "black")
    )

  save_split_panel(net_split_panel, "panel_A_loss_overview.pdf", 180, 55, split_out_dir)
  save_split_panel(hub_panel, "panel_B_consensus_hub.pdf", 180, 85, split_out_dir)
  save_split_panel(eras_split_panels[["C"]], "panel_C_eras_degree_centrality.pdf", 80, 80, split_out_dir)
  save_split_panel(eras_split_panels[["D"]], "panel_D_eras_edge_count.pdf", 80, 80, split_out_dir)
  save_split_panel(eras_split_panels[["E"]], "panel_E_eras_preserved_edges.pdf", 80, 80, split_out_dir)

  all_panels <- c(list(net_panel, hub_panel), eras_panels)
  names(all_panels) <- LETTERS[seq_along(all_panels)]
  layout_design <- "
AAAAAA
BBBBBB
CCDDEE
"

  final_plot <- wrap_plots(all_panels) +
    plot_layout(
      design = layout_design,
      heights = c(1.05, 1.55, 1.12),
      guides = "keep"
    ) +
    plot_annotation(tag_levels = "a") &
    theme(
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.3),
      legend.key.width = unit(7, "pt"),
      legend.key.height = unit(5, "pt"),
      legend.spacing.x = unit(2.5, "pt"),
      legend.box.spacing = unit(1.5, "pt"),
      plot.tag = element_text(size = 8, face = "bold", colour = "black"),
      plot.tag.position = c(0, 1)
    )

  ggsave(
    filename = file.path(out_dir, "fig_net_main_2.pdf"),
    plot = final_plot,
    width = 180,
    height = 220,
    units = "mm",
    device = cairo_pdf,
    limitsize = FALSE,
    bg = "white"
  )
  message("Saved: ", file.path(out_dir, "fig_net_main_2.pdf"))

  invisible(final_plot)
}
