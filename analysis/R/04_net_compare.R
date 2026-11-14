# Figures 4 and 4m: Control-vs-Disease edge retention.
#
# Also writes the per-disease comparison tables and the data bundles consumed
# by the final main figure.

get_edge_df <- function(g) {
  node_names <- V(g)$name
  node_category <- V(g)$abundance_category
  node_types <- ifelse(
    node_category == "high", "Rich",
    ifelse(node_category %in% c("low", "very_low"), "Rare", "Other")
  )

  edge_mat <- as_edgelist(g, names = FALSE)
  node1 <- node_names[edge_mat[, 1]]
  node2 <- node_names[edge_mat[, 2]]
  t1 <- node_types[edge_mat[, 1]]
  t2 <- node_types[edge_mat[, 2]]

  edge_id <- paste(pmin(node1, node2), pmax(node1, node2), sep = "|")

  edge_type <- ifelse(
    is.na(t1) | is.na(t2), "Others",
    ifelse(
      t1 == "Rich" & t2 == "Rich", "Rich-Rich",
      ifelse(t1 == "Rare" & t2 == "Rare", "Rare-Rare", "Rich-Rare")
    )
  )

  edge_weights <- E(g)$weight
  if (is.null(edge_weights)) {
    edge_weights <- rep(NA_real_, length(edge_type))
  }

  data.frame(
    edge_id = edge_id,
    edge_type = edge_type,
    weight = edge_weights,
    stringsAsFactors = FALSE
  ) %>%
    distinct(edge_id, .keep_all = TRUE)
}

compare_network_edges <- function(g_control, g_disease) {
  control_df <- get_edge_df(g_control)
  disease_df <- get_edge_df(g_disease)

  results <- list()
  for (etype in c("Rich-Rich", "Rare-Rare", "Rich-Rare")) {
    control_edges <- control_df %>% filter(edge_type == etype) %>% pull(edge_id)
    disease_edges <- disease_df %>% filter(edge_type == etype) %>% pull(edge_id)
    control_edges_with_weight <- control_df %>%
      filter(edge_type == etype) %>%
      dplyr::select(edge_id, weight)

    preserved <- intersect(control_edges, disease_edges)
    lost <- setdiff(control_edges, disease_edges)
    gained <- setdiff(disease_edges, control_edges)

    preserved_rate <- ifelse(
      length(control_edges) == 0, NA, length(preserved) / length(control_edges)
    )
    loss_rate <- ifelse(
      length(control_edges) == 0, NA, length(lost) / length(control_edges)
    )

    if (length(control_edges) > 0) {
      sum_rate <- preserved_rate + loss_rate
      if (abs(sum_rate - 1) > 1e-6) {
        stop(
          "PreservedRate + LossRate != 1 for edge type ", etype,
          "; PreservedRate = ", preserved_rate,
          "; LossRate = ", loss_rate,
          "; Sum = ", sum_rate
        )
      }
    }

    preserved_high_weight <- control_edges_with_weight %>%
      filter(edge_id %in% preserved, abs(weight) > 0.2) %>%
      nrow()
    control_high_weight <- control_edges_with_weight %>%
      filter(abs(weight) > 0.2) %>%
      nrow()

    results[[etype]] <- data.frame(
      EdgeType = etype,
      Control = length(control_edges),
      Disease = length(disease_edges),
      Preserved = length(preserved),
      Lost = length(lost),
      Gained = length(gained),
      PreservedRate = preserved_rate,
      LossRate = loss_rate,
      Control_HighWeight = control_high_weight,
      Preserved_HighWeight = preserved_high_weight,
      HighWeight_PreservedRate = ifelse(
        control_high_weight == 0, NA, preserved_high_weight / control_high_weight
      ),
      HighWeight_Ratio_InPreserved = ifelse(
        length(preserved) == 0, NA, preserved_high_weight / length(preserved)
      )
    )
  }

  type_df <- bind_rows(results)
  all_row <- type_df %>%
    summarise(
      EdgeType = "All",
      Control = sum(Control),
      Disease = sum(Disease),
      Preserved = sum(Preserved),
      Lost = sum(Lost),
      Gained = sum(Gained),
      PreservedRate = ifelse(Control == 0, NA, Preserved / Control),
      LossRate = ifelse(Control == 0, NA, Lost / Control),
      Control_HighWeight = sum(Control_HighWeight),
      Preserved_HighWeight = sum(Preserved_HighWeight),
      HighWeight_PreservedRate = ifelse(
        Control_HighWeight == 0, NA, Preserved_HighWeight / Control_HighWeight
      ),
      HighWeight_Ratio_InPreserved = ifelse(
        Preserved == 0, NA, Preserved_HighWeight / Preserved
      )
    )

  bind_rows(type_df, all_row)
}

compare_control_networks <- function(g_control1, g_control2, cohort1, cohort2) {
  control1_df <- get_edge_df(g_control1)
  control2_df <- get_edge_df(g_control2)

  results <- list()
  for (etype in c("Rich-Rich", "Rare-Rare", "Rich-Rare")) {
    edges1 <- control1_df %>% filter(edge_type == etype) %>% pull(edge_id)
    edges2 <- control2_df %>% filter(edge_type == etype) %>% pull(edge_id)
    common <- intersect(edges1, edges2)

    overlap_rate <- ifelse(
      length(edges1) == 0, NA, length(common) / length(edges1)
    )
    jaccard <- ifelse(
      length(union(edges1, edges2)) == 0,
      NA,
      length(common) / length(union(edges1, edges2))
    )

    results[[etype]] <- data.frame(
      EdgeType = etype,
      Control1_Edges = length(edges1),
      Control2_Edges = length(edges2),
      Common_Edges = length(common),
      Overlap_Rate = overlap_rate,
      Jaccard_Index = jaccard,
      Cohort1 = cohort1,
      Cohort2 = cohort2,
      Comparison = paste(cohort1, "vs", cohort2, sep = " ")
    )
  }
  bind_rows(results)
}

make_loss_panel <- function(bundle, show_y = TRUE) {
  stat_label1 <- paste0("Friedman p = ", signif(bundle$friedman_p, 3))
  stat_label2 <- paste0("Jonckheere p = ", signif(bundle$jonck_p, 3))

  ggplot(bundle$plot_data, aes(x = EdgeType, y = Rate)) +
    geom_line(
      aes(group = Cohort, color = Cohort),
      alpha = 0.6, size = 1, position = position_dodge(width = 0.3)
    ) +
    geom_boxplot(fill = NA, width = 0.5) +
    geom_point(
      aes(color = Cohort), size = 3, alpha = 0.9,
      position = position_dodge(width = 0.3)
    ) +
    coord_cartesian(ylim = c(0.7, 1)) +
    theme_bw() +
    labs(
      x = NULL,
      y = if (show_y) "Loss Rate" else NULL,
      title = bundle$disease
    ) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = if (show_y) element_text(size = 16, face = "bold") else element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    scale_color_brewer(palette = "Dark2", name = "Cohort") +
    annotate("text", x = -Inf, y = -Inf, label = stat_label1,
             hjust = -0.2, vjust = -1.7, size = 3.2, fontface = "italic") +
    annotate("text", x = -Inf, y = -Inf, label = stat_label2,
             hjust = -0.17, vjust = -0.5, size = 3.2, fontface = "italic")
}

compare_eras_timepoints <- function(cfg, outdir_last, out_base) {
  cohort <- cfg$eras_cohort
  output_dir <- ensure_dir(file.path(out_base, "ERAS", outdir_last))
  g_control_path <- file.path(
    network_dir(cfg), cohort, outdir_last, "g_Control_full_network.rds"
  )
  times_list <- c(
    D4 = file.path(network_dir(cfg), cohort, outdir_last, "g_D4_full_network.rds"),
    D8 = file.path(network_dir(cfg), cohort, outdir_last, "g_D8_full_network.rds"),
    D42 = file.path(network_dir(cfg), cohort, outdir_last, "g_D42_full_network.rds"),
    D180 = file.path(network_dir(cfg), cohort, outdir_last, "g_D180_full_network.rds")
  )

  all_comparison_results <- list()
  for (i in seq_along(times_list)) {
    timepoint <- names(times_list)[i]
    g_control <- readRDS(g_control_path)
    g_time <- readRDS(times_list[i])
    res <- compare_network_edges(g_control, g_time)
    res$timepoint <- timepoint
    all_comparison_results[[timepoint]] <- res
  }

  combined_results <- bind_rows(all_comparison_results) %>%
    mutate(EdgeType = rename_rich_to_abundant(EdgeType)) %>%
    mutate(
      EdgeType = factor(
        EdgeType,
        levels = c("All", "Rare-Rare", "Abundant-Rare", "Abundant-Abundant"),
        ordered = TRUE
      )
    )
  combined_results$timepoint <- factor(
    combined_results$timepoint, levels = c("D4", "D8", "D42", "D180")
  )

  p <- ggplot(
    combined_results,
    aes(x = timepoint, y = Preserved, color = EdgeType, group = EdgeType)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_text(
      aes(label = Preserved),
      vjust = -0.8, size = 4,
      position = position_jitter(width = 0.2, height = 0.1),
      show.legend = FALSE
    ) +
    theme_bw() +
    labs(x = "Timepoint", y = "Preserved Edge") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8)
    )
  ggsave(file.path(output_dir, "Preserved_zhexian.png"), p, width = 6, height = 6)
  saveRDS(combined_results, file.path(output_dir, "Preserved_zhexian.rds"))

  p2 <- ggplot(
    combined_results,
    aes(x = timepoint, y = PreservedRate, color = EdgeType, group = EdgeType)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_text(aes(label = round(PreservedRate, 3)), vjust = -0.6, size = 6) +
    theme_bw() +
    labs(
      x = "Timepoint",
      y = "PreservedRate",
      title = "Network Edge Preservation Over Time"
    )
  ggsave(file.path(output_dir, "Preserved_rate_zhexian.png"), p2,
         width = 10, height = 10)
}

analysis_disease <- function(cfg, disease, cohorts, outdir_last, out_base,
                             disease_map) {
  control_networks <- list()
  all_comparison_results <- list()
  for (cohort in cohorts) {
    g_control_path <- file.path(
      network_dir(cfg), cohort, outdir_last, "g_Control_full_network.rds"
    )
    disease_net <- cohort_disease_token(cohort)
    g_disease_path <- file.path(
      network_dir(cfg), cohort, outdir_last,
      paste0("g_", disease_net, "_full_network.rds")
    )

    if (file.exists(g_control_path) && file.exists(g_disease_path)) {
      g_control <- readRDS(g_control_path)
      g_disease <- readRDS(g_disease_path)
      control_networks[[cohort]] <- g_control
      res <- compare_network_edges(g_control, g_disease)
      res$Cohort <- cohort
      all_comparison_results[[cohort]] <- res
    } else {
      warning("Files missing for cohort: ", cohort)
    }
  }

  if (length(all_comparison_results) == 0) {
    stop("No valid cohorts found for disease: ", disease)
  }

  output_dir <- ensure_dir(file.path(out_base, disease, outdir_last))
  disease <- replace_single_disease_name(disease, disease_map)

  combined_results <- bind_rows(all_comparison_results)

  loss_data <- combined_results %>%
    dplyr::select(Cohort, EdgeType, LossRate) %>%
    dplyr::rename(Rate = LossRate)

  friedman_res <- friedman.test(Rate ~ EdgeType | Cohort, data = loss_data)
  loss_long <- loss_data %>%
    filter(EdgeType != "All") %>%
    mutate(EdgeType = factor(
      EdgeType, levels = c("Rich-Rare", "Rare-Rare", "Rich-Rich")
    ))

  trend_data <- loss_data %>%
    filter(EdgeType != "All") %>%
    mutate(EdgeType = factor(
      EdgeType,
      levels = c("Rich-Rare", "Rare-Rare", "Rich-Rich"),
      ordered = TRUE
    ))
  jonck_res <- jonckheere.test(
    trend_data$Rate, trend_data$EdgeType, alternative = "decreasing"
  )

  friedman_p <- friedman_res$p.value
  jonck_p <- jonck_res$p.value
  stat_label1 <- paste0("Friedman p = ", signif(friedman_p, 3))
  stat_label2 <- paste0("Jonckheere p = ", signif(jonck_p, 3))

  edge_colors <- c(
    "Abundant-Abundant" = "#2E8B57",
    "Rare-Rare" = "#4682B4",
    "Abundant-Rare" = "#DC143C"
  )

  plot_data_loss <- combined_results %>%
    filter(EdgeType != "All") %>%
    mutate(EdgeType = rename_rich_to_abundant(EdgeType)) %>%
    dplyr::select(Cohort, EdgeType, LossRate) %>%
    dplyr::rename(Rate = LossRate) %>%
    mutate(Metric = "Lost") %>%
    mutate(
      EdgeType = factor(
        EdgeType,
        levels = c("Rare-Rare", "Abundant-Rare", "Abundant-Abundant"),
        ordered = TRUE
      )
    )

  p_combined <- ggplot(plot_data_loss, aes(x = EdgeType, y = Rate)) +
    geom_line(
      aes(group = Cohort, color = Cohort),
      alpha = 0.6, size = 1, position = position_dodge(width = 0.3)
    ) +
    geom_boxplot(aes(fill = EdgeType), alpha = 0.7, width = 0.5) +
    geom_point(
      aes(color = Cohort), size = 3, alpha = 0.9,
      position = position_dodge(width = 0.3)
    ) +
    coord_cartesian(ylim = c(0.5, 1)) +
    theme_bw() +
    labs(x = "Edge Type", y = "Loss Rate", title = paste(disease, "Network")) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    scale_fill_manual(values = edge_colors, name = "Edge Type") +
    scale_color_brewer(palette = "Dark2", name = "Cohort") +
    annotate("text", x = -Inf, y = -Inf, label = stat_label1,
             hjust = -0.2, vjust = -1.7, size = 3.2, fontface = "italic",
             color = "black") +
    annotate("text", x = -Inf, y = -Inf, label = stat_label2,
             hjust = -0.17, vjust = -0.5, size = 3.2, fontface = "italic",
             color = "black")

  ggsave(file.path(output_dir, "Combined_Edge_Preservation_Boxplot.png"),
         p_combined, width = 8, height = 6, dpi = 300)

  p_connected <- loss_long %>%
    ggplot(aes(x = EdgeType, y = Rate, group = Cohort, color = Cohort)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3) +
    stat_summary(aes(group = 1), fun = mean, geom = "line",
                 color = "black", size = 1.2, linetype = "dashed") +
    theme_bw() +
    labs(y = "Loss Rate", x = "Edge Type", color = "Cohort") +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set1") +
    annotate("text", x = Inf, y = Inf,
             label = paste("Friedman p =", signif(friedman_p, 3)),
             hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "italic")
  ggsave(file.path(output_dir, "Combined_Edge_Loss_Boxplot_connect.png"),
         p_connected, width = 8, height = 6, dpi = 300)

  bundle <- list(
    disease = disease,
    plot_data = plot_data_loss,
    friedman_p = friedman_p,
    jonck_p = jonck_p
  )

  write.csv(combined_results,
            file.path(output_dir, "Control_vs_Disease_Comparison_Summary.csv"),
            row.names = FALSE)

  if (length(control_networks) >= 2) {
    cohort_names <- names(control_networks)
    control_comp_list <- list()
    for (i in 1:(length(cohort_names) - 1)) {
      for (j in (i + 1):length(cohort_names)) {
        cohort1 <- cohort_names[i]
        cohort2 <- cohort_names[j]
        control_comp_list[[paste(cohort1, cohort2, sep = "_vs_")]] <-
          compare_control_networks(
            control_networks[[cohort1]],
            control_networks[[cohort2]],
            cohort1, cohort2
          )
      }
    }
    control_comparison_results <- bind_rows(control_comp_list)

    p_control_overlap <- ggplot(
      control_comparison_results,
      aes(x = EdgeType, y = Overlap_Rate, fill = EdgeType)
    ) +
      geom_boxplot(alpha = 0.7) +
      geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
      theme_bw() +
      labs(x = "Edge Type", y = "Overlap Rate",
           title = "Control Network Consistency Across Cohorts") +
      ylim(0, 1) +
      theme(legend.position = "bottom") +
      scale_fill_brewer(palette = "Set2")
    ggsave(file.path(output_dir, "Control_Network_Overlap_Boxplot.png"),
           p_control_overlap, width = 6, height = 6, dpi = 300)

    p_jaccard <- ggplot(
      control_comparison_results,
      aes(x = EdgeType, y = Jaccard_Index, fill = EdgeType)
    ) +
      geom_boxplot(alpha = 0.7) +
      geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
      theme_bw() +
      labs(x = "Edge Type", y = "Jaccard Index",
           title = "Control Network Similarity (Jaccard Index) Across Cohorts") +
      ylim(0, 1) +
      theme(legend.position = "bottom") +
      scale_fill_brewer(palette = "Set2")
    ggsave(file.path(output_dir, "Control_Network_Jaccard_Boxplot.png"),
           p_jaccard, width = 6, height = 6, dpi = 300)

    write.csv(control_comparison_results,
              file.path(output_dir, "Control_Network_Comparisons.csv"),
              row.names = FALSE)

    heatmaps <- list()
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      edge_type <- "Rare-Rare"
      richrich_overlap <- control_comparison_results %>%
        filter(EdgeType == edge_type) %>%
        select(Cohort1, Cohort2, Overlap_Rate)

      if (nrow(richrich_overlap) > 0) {
        cohort_list <- unique(c(richrich_overlap$Cohort1, richrich_overlap$Cohort2))
        overlap_matrix <- matrix(
          NA, nrow = length(cohort_list), ncol = length(cohort_list)
        )
        rownames(overlap_matrix) <- cohort_list
        colnames(overlap_matrix) <- cohort_list
        for (i in seq_len(nrow(richrich_overlap))) {
          overlap_matrix[richrich_overlap$Cohort1[i], richrich_overlap$Cohort2[i]] <-
            richrich_overlap$Overlap_Rate[i]
          overlap_matrix[richrich_overlap$Cohort2[i], richrich_overlap$Cohort1[i]] <-
            richrich_overlap$Overlap_Rate[i]
        }
        diag(overlap_matrix) <- 1

        pheatmap::pheatmap(
          overlap_matrix,
          main = sprintf("Control Network Overlap Rate (%s Edges)", edge_type),
          cluster_rows = TRUE, cluster_cols = TRUE,
          display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
          filename = file.path(
            output_dir, paste0(edge_type, "Control_Network_Overlap_Heatmap.png")
          ),
          width = 10, height = 8
        )

        for (cohort in names(all_comparison_results)) {
          rr_preserved <- all_comparison_results[[cohort]] %>%
            filter(EdgeType == edge_type) %>%
            pull(PreservedRate)
          overlap_matrix[cohort, cohort] <- rr_preserved
        }

        p_heatmap2 <- pheatmap::pheatmap(
          overlap_matrix,
          main = sprintf("Control Network Overlap Rate (%s Edges)", edge_type),
          cluster_rows = TRUE, cluster_cols = TRUE,
          display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
          color = colorRampPalette(c("#EAF2FF", "#e81527ff"))(50),
          filename = file.path(
            output_dir, paste0(edge_type, "Control_Network_Overlap_Heatmap2.png")
          ),
          width = 10, height = 8
        )
        heatmaps[[1]] <- p_heatmap2[[4]]
      }
    } else {
      warning("pheatmap is not installed; control-network heatmaps are skipped")
    }
  } else {
    heatmaps <- list()
    message("Skipping control network comparison: need at least 2 cohorts (found ",
            length(control_networks), ")")
  }

  per_cohort_stats <- combined_results %>%
    select(Cohort, EdgeType, Control, Disease, Preserved, Lost, Gained,
           PreservedRate, LossRate) %>%
    arrange(Cohort, EdgeType)
  write.csv(per_cohort_stats,
            file.path(output_dir, "Per_Cohort_Details.csv"),
            row.names = FALSE)

  list(
    bundle = bundle,
    heatmaps = heatmaps
  )
}

draw_combined_loss <- function(data_list, out_file, n_col) {
  if (length(data_list) == 0) return(invisible())
  first_in_row <- ((seq_along(data_list) - 1) %% n_col) == 0
  panels <- Map(function(b, sy) make_loss_panel(b, show_y = sy),
                data_list, first_in_row)
  n_row <- ceiling(length(panels) / n_col)
  combined <- wrap_plots(panels, ncol = n_col, nrow = n_row) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  ggsave(out_file, combined, width = 3 * n_col, height = 5 * n_row)
}

run_net_compare <- function(cfg) {
  suffixes <- analysis_suffixes(cfg)
  out_base <- ensure_dir(file.path(cfg$results_dir, cfg$paths$net_compare))
  sig_alpha <- 0.05
  disease_map <- read_disease_map(cfg)

  for (outdir_last in suffixes) {
    main_data_list <- list()
    supp_data_list <- list()
    heatmaps <- list()
    disease_cohorts <- get_disease_cohorts(network_dir(cfg))

    for (disease in names(disease_cohorts)) {
      if (disease %in% cfg$skip_diseases$net_compare) next

      if (disease == "ERAS") {
        compare_eras_timepoints(cfg, outdir_last, out_base)
      }

      cohorts <- disease_cohorts[[disease]]
      if (length(cohorts) < 2) {
        if (disease == "ERAS") next
        unlink(file.path(out_base, disease), recursive = TRUE, force = TRUE)
        message("Removed empty directory: ", file.path(out_base, disease))
        next
      }

      result <- analysis_disease(
        cfg, disease, cohorts, outdir_last, out_base,
        disease_map = disease_map
      )
      if (isTRUE(result$bundle$friedman_p < sig_alpha) &&
          isTRUE(result$bundle$jonck_p < sig_alpha)) {
        main_data_list[[result$bundle$disease]] <- result$bundle
      } else {
        supp_data_list[[result$bundle$disease]] <- result$bundle
      }
      heatmaps <- c(heatmaps, result$heatmaps)
    }

    combined_output_dir <- ensure_dir(file.path(out_base, "combined"))
    saveRDS(main_data_list,
            file.path(combined_output_dir, paste0(outdir_last, "_main_data.rds")))
    saveRDS(supp_data_list,
            file.path(combined_output_dir, paste0(outdir_last, "_supp_data.rds")))

    message("Main-figure diseases (Friedman & Jonckheere < ", sig_alpha, "): ",
            paste(names(main_data_list), collapse = ", "))
    message("Supplementary diseases: ",
            paste(names(supp_data_list), collapse = ", "))

    draw_combined_loss(
      main_data_list,
      file.path(cfg$results_dir, cfg$paths$final, "fig_net_4m.pdf"),
      n_col = 3
    )
    draw_combined_loss(
      supp_data_list,
      file.path(cfg$results_dir, cfg$paths$final, "fig_net_4.pdf"),
      n_col = 4
    )

    if (length(heatmaps) > 0) {
      n_plots <- length(heatmaps)
      n_col <- ceiling(sqrt(n_plots))
      n_row <- ceiling(n_plots / n_col)
      combined <- do.call(
        gridExtra::arrangeGrob,
        c(heatmaps, ncol = n_col, nrow = n_row,
          top = paste("Network Heatmaps -", outdir_last))
      )
      ggsave(
        file.path(combined_output_dir, paste0(outdir_last, "_net_heatmap.png")),
        combined, width = 32, height = 24, dpi = 300, limitsize = FALSE, bg = "white"
      )
    }
  }

  invisible(TRUE)
}
