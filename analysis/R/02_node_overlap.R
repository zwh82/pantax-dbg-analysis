# Figure 2: node overlap between Control and Disease networks.

get_node_df <- function(g) {
  node_names <- V(g)$name
  node_category <- V(g)$abundance_category
  node_types <- ifelse(
    node_category == "high", "Rich",
    ifelse(node_category %in% c("low", "very_low"), "Rare", "Other")
  )
  data.frame(
    node_id = node_names,
    node_type = node_types,
    stringsAsFactors = FALSE
  )
}

compare_node_overlap <- function(g_control, g_disease) {
  control_nodes <- get_node_df(g_control)
  disease_nodes <- get_node_df(g_disease)

  control_all <- control_nodes$node_id
  disease_all <- disease_nodes$node_id
  common_all <- intersect(control_all, disease_all)
  preserved_rate_all <- length(common_all) / length(control_all)

  results <- list()
  for (ntype in c("Rich", "Rare", "Other")) {
    control_subset <- control_nodes %>% filter(node_type == ntype) %>% pull(node_id)
    disease_subset <- disease_nodes %>% filter(node_type == ntype) %>% pull(node_id)

    if (length(control_subset) == 0) {
      results[[ntype]] <- data.frame(
        NodeType = ntype,
        Control_Count = 0,
        Disease_Count = length(disease_subset),
        Common_Count = 0,
        PreservedRate = NA_real_
      )
    } else {
      common <- intersect(control_subset, disease_subset)
      results[[ntype]] <- data.frame(
        NodeType = ntype,
        Control_Count = length(control_subset),
        Disease_Count = length(disease_subset),
        Common_Count = length(common),
        PreservedRate = length(common) / length(control_subset)
      )
    }
  }

  results[["All"]] <- data.frame(
    NodeType = "All",
    Control_Count = length(control_all),
    Disease_Count = length(disease_all),
    Common_Count = length(common_all),
    PreservedRate = preserved_rate_all
  )

  bind_rows(results)
}

collect_node_overlap <- function(cfg, outdir_last) {
  disease_cohorts <- get_disease_cohorts(network_dir(cfg))
  all_disease_results <- list()

  for (disease in names(disease_cohorts)) {
    if (disease %in% cfg$skip_diseases$node) next
    cohorts <- disease_cohorts[[disease]]
    if (length(cohorts) < 1) next

    disease_results <- list()
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
        overlap_res <- compare_node_overlap(
          readRDS(g_control_path),
          readRDS(g_disease_path)
        )
        overlap_res$Cohort <- cohort
        overlap_res$Disease <- disease
        disease_results[[cohort]] <- overlap_res
      } else {
        warning("Files missing for cohort: ", cohort)
      }
    }

    if (length(disease_results) > 0) {
      all_disease_results[[disease]] <- bind_rows(disease_results)
    }
  }
  all_disease_results
}

make_node_barplot <- function(bar_data) {
  ggplot(bar_data, aes(x = Disease, y = MeanPreservedRate, fill = MeanPreservedRate)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_errorbar(
      aes(ymin = MeanPreservedRate - SD, ymax = MeanPreservedRate + SD),
      width = 0.2
    ) +
    geom_text(
      aes(label = sprintf("%.2f", MeanPreservedRate)),
      vjust = -0.5,
      size = 2.5,
      color = "black"
    ) +
    scale_fill_gradient(low = "#6ea4da", high = "#f09faa", name = "Rate") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    labs(x = "Disease", y = "Species Overlap Rate")
}

run_node_overlap <- function(cfg) {
  suffixes <- analysis_suffixes(cfg)
  disease_map <- read_disease_map(cfg)

  for (outdir_last in suffixes) {
    node_overlap_dir <- ensure_dir(
      file.path(cfg$results_dir, cfg$paths$node_overlap, outdir_last)
    )
    all_disease_results <- collect_node_overlap(cfg, outdir_last)

    if (length(all_disease_results) == 0) {
      message("No node-overlap results found for ", outdir_last)
      next
    }

    for (disease in names(all_disease_results)) {
      write.csv(
        all_disease_results[[disease]],
        file.path(node_overlap_dir, paste0(disease, "_node_overlap_results.csv")),
        row.names = FALSE
      )
    }

    bar_data <- bind_rows(all_disease_results) %>%
      filter(NodeType %in% c("All"), !is.na(PreservedRate)) %>%
      replace_disease_name(disease_map) %>%
      group_by(Disease) %>%
      summarise(
        n = n(),
        MeanPreservedRate = mean(PreservedRate, na.rm = TRUE),
        SD = sd(PreservedRate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n)) %>%
      mutate(Disease = factor(Disease, levels = Disease))

    if (nrow(bar_data) > 0) {
      p_bar <- make_node_barplot(bar_data)
      ggsave(file.path(node_overlap_dir, "Node_Preservation_Barplot.png"),
             p_bar, width = 12, height = 8, dpi = 300)
      ggsave(file.path(cfg$results_dir, cfg$paths$final, "fig_net_2.pdf"),
             p_bar, width = 9, height = 6)
    }

    summary_nodes <- bind_rows(all_disease_results) %>%
      group_by(Disease) %>%
      summarise(
        Mean_PreservedRate = mean(PreservedRate, na.rm = TRUE),
        SD_PreservedRate = sd(PreservedRate, na.rm = TRUE),
        Min_PreservedRate = min(PreservedRate, na.rm = TRUE),
        Max_PreservedRate = max(PreservedRate, na.rm = TRUE),
        NumberOfCohorts = n(),
        .groups = "drop"
      )
    write.csv(summary_nodes,
              file.path(node_overlap_dir, "Node_Overlap_Summary.csv"),
              row.names = FALSE)
  }

  invisible(TRUE)
}
