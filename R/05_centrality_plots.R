# Centrality comparisons between Rich and Rare species.

perform_multi_group_centrality_test <- function(data, metric, group_names) {
  results <- list()
  abundance_types <- c("Rich", "Rare")

  global_tests <- data.frame(
    AbundanceType = character(),
    Metric = character(),
    Statistic = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )

  pairwise_p <- list()

  for (atype in abundance_types) {
    subset_data <- data[data$AbundanceType == atype, ]

    if (length(unique(subset_data$Group)) >= 2) {
      if (length(unique(subset_data$Group)) > 2) {
        kw_res <- kruskal.test(subset_data[[metric]] ~ Group, data = subset_data)
      } else {
        kw_res <- wilcox.test(subset_data[[metric]] ~ Group, data = subset_data)
      }

      global_tests <- rbind(global_tests, data.frame(
        AbundanceType = atype,
        Metric = metric,
        Statistic = kw_res$statistic,
        P_value = kw_res$p.value
      ))

      groups_present <- unique(subset_data$Group)
      pairwise_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names))
      rownames(pairwise_matrix) <- group_names
      colnames(pairwise_matrix) <- group_names

      for (i in 1:(length(groups_present) - 1)) {
        for (j in (i + 1):length(groups_present)) {
          grp1 <- groups_present[i]
          grp2 <- groups_present[j]

          val1 <- if (metric == "Betweenness") {
            subset_data[[metric]][subset_data$Group == grp1] + 1
          } else {
            subset_data[[metric]][subset_data$Group == grp1]
          }
          val2 <- if (metric == "Betweenness") {
            subset_data[[metric]][subset_data$Group == grp2] + 1
          } else {
            subset_data[[metric]][subset_data$Group == grp2]
          }

          if (length(val1) > 0 && length(val2) > 0) {
            wilcox_res <- wilcox.test(val1, val2)
            pairwise_matrix[grp1, grp2] <- pairwise_matrix[grp2, grp1] <- wilcox_res$p.value
          }
        }
      }

      pairwise_p[[atype]] <- pairwise_matrix
    }
  }

  results$global <- global_tests
  results$pairwise_p <- pairwise_p
  results
}

add_significance_labels <- function(plot, pairwise_p_list, abundance_types, y_pos) {
  for (atype in abundance_types) {
    if (!is.null(pairwise_p_list[[atype]])) {
      p_matrix <- pairwise_p_list[[atype]]
      min_p <- min(p_matrix, na.rm = TRUE)

      if (min_p < 0.001) {
        label <- "***"
      } else if (min_p < 0.01) {
        label <- "**"
      } else if (min_p < 0.05) {
        label <- "*"
      } else {
        label <- sprintf("p = %.2e", min_p)
      }

      x_pos <- which(abundance_types == atype)
      plot <- plot + annotate("text", x = x_pos, y = y_pos,
                              label = label, size = 3.5, fontface = "italic")
    }
  }
  plot
}

generate_centrality_summary <- function(processed_data_list, group_names) {
  summary_stats <- data.frame()

  for (grp in group_names) {
    if (!is.null(processed_data_list[[grp]])) {
      df <- processed_data_list[[grp]]

      for (atype in c("Rich", "Rare")) {
        subset_df <- df[df$AbundanceType == atype, ]
        if (nrow(subset_df) > 0) {
          summary_stats <- rbind(summary_stats, data.frame(
            Group = grp,
            AbundanceType = atype,
            Degree_Mean = mean(subset_df$Degree, na.rm = TRUE),
            Degree_Median = median(subset_df$Degree, na.rm = TRUE),
            Degree_SD = sd(subset_df$Degree, na.rm = TRUE),
            Betweenness_Mean = mean(subset_df$Betweenness, na.rm = TRUE),
            Betweenness_Median = median(subset_df$Betweenness, na.rm = TRUE),
            Betweenness_SD = sd(subset_df$Betweenness, na.rm = TRUE),
            Closeness_Mean = mean(subset_df$Closeness, na.rm = TRUE),
            Closeness_Median = median(subset_df$Closeness, na.rm = TRUE),
            Closeness_SD = sd(subset_df$Closeness, na.rm = TRUE),
            N = nrow(subset_df)
          ))
        }
      }
    }
  }

  summary_stats
}

output_test_results <- function(test_results, outdir, group_names) {
  sink(paste0(outdir, "Multi_Group_Centrality_Statistics.txt"))
  cat("Multi-Group Centrality Analysis by Abundance Type\n")
  cat("==================================================\n\n")
  cat("Groups analyzed:", paste(group_names, collapse = ", "), "\n\n")

  if (!is.null(test_results$degree)) {
    cat("\n>>> DEGREE CENTRALITY <<<\n")
    print(test_results$degree$global)

    cat("\nPairwise comparisons (Wilcoxon test):\n")
    for (atype in names(test_results$degree$pairwise_p)) {
      cat(sprintf("\n  %s species:\n", atype))
      p_mat <- test_results$degree$pairwise_p[[atype]]
      for (i in 1:(nrow(p_mat) - 1)) {
        for (j in (i + 1):ncol(p_mat)) {
          if (!is.na(p_mat[i, j])) {
            cat(sprintf("    %s vs %s: p = %.2e\n",
                        rownames(p_mat)[i], colnames(p_mat)[j], p_mat[i, j]))
          }
        }
      }
    }
  }

  if (!is.null(test_results$betweenness)) {
    cat("\n\n>>> BETWEENNESS CENTRALITY <<<\n")
    print(test_results$betweenness$global)

    cat("\nPairwise comparisons (Wilcoxon test):\n")
    for (atype in names(test_results$betweenness$pairwise_p)) {
      cat(sprintf("\n  %s species:\n", atype))
      p_mat <- test_results$betweenness$pairwise_p[[atype]]
      for (i in 1:(nrow(p_mat) - 1)) {
        for (j in (i + 1):ncol(p_mat)) {
          if (!is.na(p_mat[i, j])) {
            cat(sprintf("    %s vs %s: p = %.2e\n",
                        rownames(p_mat)[i], colnames(p_mat)[j], p_mat[i, j]))
          }
        }
      }
    }
  }

  if (!is.null(test_results$closeness)) {
    cat("\n\n>>> CLOSENESS CENTRALITY <<<\n")
    print(test_results$closeness$global)

    cat("\nPairwise comparisons (Wilcoxon test):\n")
    for (atype in names(test_results$closeness$pairwise_p)) {
      cat(sprintf("\n  %s species:\n", atype))
      p_mat <- test_results$closeness$pairwise_p[[atype]]
      for (i in 1:(nrow(p_mat) - 1)) {
        for (j in (i + 1):ncol(p_mat)) {
          if (!is.na(p_mat[i, j])) {
            cat(sprintf("    %s vs %s: p = %.2e\n",
                        rownames(p_mat)[i], colnames(p_mat)[j], p_mat[i, j]))
          }
        }
      }
    }
  }

  sink()
}

plot_multi_abundance_centrality_boxplots <- function(centrality_data_list, group_names,
                                                     outdir, eras = FALSE) {
  cat("\n========== Plotting Multi-Group Centrality Boxplots by Abundance Category ==========\n")

  classify_abundance_type <- function(df) {
    df$AbundanceType <- abundance_type_from_category(df$AbundanceCategory)
    df
  }

  processed_data_list <- list()
  combined_data <- data.frame()

  n_groups <- length(group_names)
  if (n_groups <= 3) {
    colors <- c("#4DAF4A", "#E41A1C", "#377EB8")
  } else if (n_groups <= 8) {
    colors <- RColorBrewer::brewer.pal(n_groups, "Set1")
  } else {
    colors <- rainbow(n_groups)
  }
  names(colors) <- group_names

  test_results <- list()

  for (grp in group_names) {
    if (!is.null(centrality_data_list[[grp]])) {
      df <- centrality_data_list[[grp]]
      df <- classify_abundance_type(df)
      df <- df[df$AbundanceType != "Other", ]
      df$Group <- grp
      processed_data_list[[grp]] <- df
      combined_data <- rbind(combined_data, df)
    }
  }

  if (length(unique(combined_data$Group)) < 2) {
    cat("[WARNING] Need at least 2 groups for comparison\n")
    return(NULL)
  }

  combined_data$AbundanceType <- factor(combined_data$AbundanceType,
                                        levels = c("Rich", "Rare"))
  if ("ACA" %in% group_names) {
    combined_data$Group <- factor(combined_data$Group, levels = c("Control", "CRC", "ACA"))
  } else if (eras) {
    combined_data$Group <- factor(combined_data$Group,
                                  levels = c("Control", "D4", "D8", "D42", "D180"))
  }

  cat("\nSample sizes:\n")
  for (grp in group_names) {
    if (grp %in% names(processed_data_list)) {
      df <- processed_data_list[[grp]]
      cat(sprintf("  %s - Rich: %d, Rare: %d\n",
                  grp, sum(df$AbundanceType == "Rich"), sum(df$AbundanceType == "Rare")))
    }
  }

  p1 <- ggplot(combined_data, aes(x = AbundanceType, y = Degree, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, width = 0.7,
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = colors) +
    labs(title = "Degree Centrality: Rich vs Rare Species Across Groups",
         x = "Abundance Type",
         y = "Degree Centrality (log10 scale)",
         fill = "Group") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray50", linewidth = 0.5)
    ) +
    scale_y_log10() +
    annotation_logticks(sides = "l")

  degree_test_results <- perform_multi_group_centrality_test(combined_data, "Degree", group_names)
  test_results$degree <- degree_test_results

  max_y_degree <- max(log10(combined_data$Degree + 1), na.rm = TRUE)
  y_pos_degree <- 10^max_y_degree * 1.2

  p1 <- add_significance_labels(p1, degree_test_results$pairwise_p,
                                levels(combined_data$AbundanceType), y_pos_degree)

  saveRDS(p1, file = paste0(outdir, "/fig_rds/Multi_Group_Degree_Centrality_Boxplot.rds"))
  saveRDS(combined_data, file = paste0(outdir, "/data/combined_data.rds"))

  p2 <- ggplot(combined_data, aes(x = AbundanceType, y = Betweenness + 1, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, width = 0.7,
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = colors) +
    labs(title = "Betweenness Centrality: Rich vs Rare Species Across Groups",
         x = "Abundance Type",
         y = "Betweenness Centrality (log10 scale)",
         fill = "Group") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray50", linewidth = 0.5)
    ) +
    scale_y_log10() +
    annotation_logticks(sides = "l")

  betweenness_test_results <- perform_multi_group_centrality_test(
    combined_data, "Betweenness", group_names
  )
  test_results$betweenness <- betweenness_test_results

  max_y_between <- max(log10(combined_data$Betweenness + 1), na.rm = TRUE)
  y_pos_between <- 10^max_y_between * 1.2

  p2 <- add_significance_labels(p2, betweenness_test_results$pairwise_p,
                                levels(combined_data$AbundanceType), y_pos_between)

  saveRDS(p2, file = paste0(outdir, "/fig_rds/Multi_Group_Betweenness_Centrality_Boxplot.rds"))

  thresholds <- combined_data %>%
    filter(!is.na(Group)) %>%
    group_by(Group) %>%
    summarise(
      Q1 = quantile(Closeness, 0.25, na.rm = TRUE),
      Q3 = quantile(Closeness, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - 1.5 * IQR,
      upper_bound = Q3 + 1.5 * IQR,
      .groups = "drop"
    )

  ylim_lower <- min(thresholds$lower_bound, na.rm = TRUE)
  ylim_upper <- max(thresholds$upper_bound, na.rm = TRUE)

  p3 <- ggplot(combined_data %>% filter(!is.na(Group)),
               aes(x = AbundanceType, y = Closeness, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.shape = NA, width = 0.7,
                 position = position_dodge(0.8)) +
    scale_y_continuous(limits = c(ylim_lower, ylim_upper), oob = scales::squish) +
    scale_fill_manual(values = colors) +
    labs(title = "Closeness Centrality: Rich vs Rare Species Across Groups",
         x = "Abundance Type",
         y = "Closeness Centrality",
         fill = "Group") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray50", linewidth = 0.5)
    )

  closeness_test_results <- perform_multi_group_centrality_test(
    combined_data, "Closeness", group_names
  )
  test_results$closeness <- closeness_test_results

  y_pos_close <- ylim_upper * 1.1
  p3 <- add_significance_labels(p3, closeness_test_results$pairwise_p,
                                levels(combined_data$AbundanceType), y_pos_close)

  saveRDS(p3, file = paste0(outdir, "/fig_rds/Multi_Group_Closeness_Centrality_Boxplot.rds"))

  plot_width <- 22
  plot_height <- 5
  combined_plot <- plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, align = "h")

  ggsave(paste0(outdir, "Multi_Group_Centrality_Rich_vs_Rare_Boxplot.png"),
         combined_plot, width = plot_width, height = plot_height, dpi = 300)

  summary_stats <- generate_centrality_summary(processed_data_list, group_names)
  write.csv(summary_stats,
            file = paste0(outdir, "Multi_Group_Centrality_Rich_vs_Rare_Summary.csv"),
            row.names = FALSE)

  cat("\nSummary statistics saved to:",
      paste0(outdir, "Multi_Group_Centrality_Rich_vs_Rare_Summary.csv"), "\n")
  print(summary_stats)

  output_test_results(test_results, outdir, group_names)

  list(
    test_results = test_results,
    summary_stats = summary_stats,
    plots = list(degree = p1, betweenness = p2, closeness = p3, combined = combined_plot),
    combined_data = combined_data
  )
}
