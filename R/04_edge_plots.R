# Edge-level summaries and multi-group edge plots.

extract_edge_info <- function(g, group_name) {
  edge_weights <- E(g)$weight
  edge_abs_weights <- abs(edge_weights)
  cat(group_name, "max edge weight: ", max(edge_abs_weights), "\n")

  edge_df <- data.frame(
    Group = group_name,
    Edge_Weight = edge_weights,
    Absolute_Weight = edge_abs_weights,
    Edge_Type = ifelse(edge_weights > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )

  edge_df$Strength_Category <- cut(
    edge_df$Absolute_Weight,
    breaks = c(0, 0.1, 0.3, 0.5, 1),
    labels = c("Very Weak", "Weak", "Moderate", "Strong")
  )

  edge_df
}

classify_abundance_for_edges <- function(g) {
  V(g)$AbundanceType <- abundance_type_from_category(V(g)$abundance_category)
  g
}

classify_edges <- function(g, group_name) {
  node_types <- V(g)$AbundanceType
  edge_mat <- as_edgelist(g, names = FALSE)

  t1 <- node_types[edge_mat[, 1]]
  t2 <- node_types[edge_mat[, 2]]

  edge_type <- ifelse(
    is.na(t1) | is.na(t2), "Others",
    ifelse(t1 == "Rich" & t2 == "Rich", "Rich-Rich",
           ifelse(t1 == "Rare" & t2 == "Rare", "Rare-Rare", "Rich-Rare"))
  )

  data.frame(
    Group = group_name,
    Edge_Type = edge_type,
    Weight = E(g)$weight,
    Abs_Weight = abs(E(g)$weight),
    stringsAsFactors = FALSE
  )
}

calculate_rho_by_edgetype <- function(group_edges_list, group_names) {
  target_types <- c("Rich-Rich", "Rich-Rare", "Rare-Rare")
  rho_results <- data.frame(
    Group = character(),
    Rich_Rich = numeric(),
    Rich_Rare = numeric(),
    Rare_Rare = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp in group_names) {
    if (!is.null(group_edges_list[[grp]])) {
      edges <- group_edges_list[[grp]]

      result_row <- data.frame(Group = grp, Rich_Rich = NA, Rich_Rare = NA, Rare_Rare = NA)

      for (etype in target_types) {
        subset_edges <- edges[edges$Edge_Type == etype, ]

        if (nrow(subset_edges) > 0) {
          pos_sum <- sum(subset_edges$Weight[subset_edges$Weight > 0])
          neg_sum <- sum(subset_edges$Weight[subset_edges$Weight < 0])

          rho <- if (pos_sum == 0 && neg_sum == 0) {
            NA
          } else {
            (pos_sum + neg_sum) / (pos_sum - neg_sum)
          }

          result_row[[gsub("-", "_", etype)]] <- rho
        }
      }

      rho_results <- rbind(rho_results, result_row)
    }
  }

  rho_results
}

plot_multi_edges <- function(graphs, group_names, outdir, eras = FALSE) {
  colors <- setNames(ggsci::pal_npg("nrc")(length(group_names)), group_names)

  group_edges_list <- list()

  for (grp in group_names) {
    if (!is.null(graphs[[grp]]) && ecount(graphs[[grp]]) > 0) {
      group_edges_list[[grp]] <- extract_edge_info(graphs[[grp]], grp)
    }
  }

  if (length(group_edges_list) < 2) {
    cat("[WARNING] Need at least 2 groups for comparison\n")
    return(NULL)
  }

  all_edges <- do.call(rbind, group_edges_list)
  if ("ACA" %in% group_names) {
    all_edges$Group <- factor(all_edges$Group, levels = c("Control", "ACA", "CRC"))
  } else if (eras) {
    all_edges$Group <- factor(all_edges$Group, levels = c("Control", "D4", "D8", "D42", "D180"))
  }

  rho_results <- data.frame(Group = character(), Rho = numeric(), stringsAsFactors = FALSE)

  for (grp in group_names) {
    if (!is.null(group_edges_list[[grp]])) {
      edges <- group_edges_list[[grp]]
      pos_sum <- sum(edges$Edge_Weight[edges$Edge_Type == "Positive"])
      neg_sum <- sum(edges$Edge_Weight[edges$Edge_Type == "Negative"])
      rho <- (pos_sum + neg_sum) / (pos_sum - neg_sum)
      rho_results <- rbind(rho_results, data.frame(Group = grp, Rho = rho))
    }
  }

  write.csv(rho_results, paste0(outdir, "data/rho_results_multi.csv"), row.names = FALSE)

  edge_counts <- data.frame()
  for (grp in group_names) {
    if (!is.null(group_edges_list[[grp]])) {
      edges <- group_edges_list[[grp]]
      edge_counts <- rbind(edge_counts,
        data.frame(
          Group = grp,
          Edge_Type = c("Positive", "Negative"),
          Count = c(sum(edges$Edge_Type == "Positive"),
                    sum(edges$Edge_Type == "Negative"))
        )
      )
    }
  }
  if ("ACA" %in% group_names) {
    edge_counts$Group <- factor(edge_counts$Group, levels = c("Control", "ACA", "CRC"))
  } else if (eras) {
    edge_counts$Group <- factor(edge_counts$Group, levels = c("Control", "D4", "D8", "D42", "D180"))
  }
  write.csv(edge_counts, paste0(outdir, "data/edge_counts_multi.csv"), row.names = FALSE)

  contingency_table <- matrix(edge_counts$Count, nrow = length(group_names), byrow = TRUE)
  rownames(contingency_table) <- unique(edge_counts$Group)
  colnames(contingency_table) <- c("Positive", "Negative")
  chi_res <- chisq.test(contingency_table)
  p_val_edge <- chi_res$p.value
  cat("\nEdge Composition Chi-square test p-value:",
      format(p_val_edge, scientific = TRUE), "\n")

  p_count <- ggplot(edge_counts, aes(x = Edge_Type, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = colors) +
    labs(title = "Comparison of Interaction Types Across Groups",
         x = "Interaction Type (Edge Type)",
         y = "Number of Edges",
         fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "top") +
    geom_text(aes(label = Count), position = position_dodge(0.8), vjust = -0.5, size = 3.5) +
    annotate("text", x = 1.5, y = max(edge_counts$Count) * 1.1,
             label = paste0("Chi-square test\np = ", format(p_val_edge, scientific = TRUE)),
             size = 3, fontface = "italic") +
    expand_limits(y = max(edge_counts$Count) * 1.2)

  saveRDS(p_count, file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Count_Comparison_plot.rds"))

  all_edges$Absolute_Weight <- abs(all_edges$Edge_Weight)

  p_strength <- ggplot(all_edges, aes(x = Group, y = Absolute_Weight, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
    scale_fill_manual(values = colors) +
    scale_y_log10() +
    labs(title = "Interaction Strength Across Groups",
         x = "Group",
         y = "Absolute Edge Weight (log10)",
         fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 10),
          legend.position = "top") +
    annotation_logticks(sides = "l")

  if (length(group_names) > 2) {
    test_name <- "KW"
    kw_test <- kruskal.test(Absolute_Weight ~ Group, data = all_edges)
  } else {
    test_name <- "U"
    kw_test <- wilcox.test(Absolute_Weight ~ Group, data = all_edges)
  }

  p_strength <- p_strength +
    annotate("text", x = length(group_names) / 2 + 0.5,
             y = max(all_edges$Absolute_Weight),
             label = sprintf("%s test: p = %.2e", test_name, kw_test$p.value),
             size = 4)

  saveRDS(p_strength, file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Strength_Boxplot.rds"))

  positive_edges <- all_edges[all_edges$Edge_Type == "Positive", ]
  if (nrow(positive_edges) > 0) {
    p_pos <- ggplot(positive_edges, aes(x = Group, y = Edge_Weight, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
      scale_fill_manual(values = colors) +
      labs(title = "Positive Interaction Strength Across Groups",
           x = "Group",
           y = "Edge Weight",
           fill = "Group") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 10),
            legend.position = "top")

    if (length(unique(positive_edges$Group)) > 2) {
      pos_test <- kruskal.test(Edge_Weight ~ Group, data = positive_edges)
      test_name <- "KW"
    } else if (length(unique(positive_edges$Group)) == 2) {
      pos_test <- wilcox.test(Edge_Weight ~ Group, data = positive_edges)
      test_name <- "U"
    } else {
      pos_test <- list(p.value = NA)
      test_name <- "N"
    }
    p_pos <- p_pos + annotate("text", x = length(unique(positive_edges$Group)) / 2 + 0.5,
                              y = max(positive_edges$Edge_Weight),
                              label = sprintf("%s test: p = %.2e", test_name, pos_test$p.value),
                              size = 4)
    saveRDS(p_pos, file = paste0(outdir, "/fig_rds/Multi_Group_Positive_Edge_Strength_Boxplot.rds"))
  } else {
    p_pos <- NULL
  }

  negative_edges <- all_edges[all_edges$Edge_Type == "Negative", ]
  if (nrow(negative_edges) > 0) {
    negative_edges$Abs_Weight <- abs(negative_edges$Edge_Weight)

    p_neg <- ggplot(negative_edges, aes(x = Group, y = Abs_Weight, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
      scale_fill_manual(values = colors) +
      labs(title = "Negative Interaction Strength (Absolute) Across Groups",
           x = "Group",
           y = "Absolute Edge Weight",
           fill = "Group") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 10),
            legend.position = "top")

    if (length(unique(negative_edges$Group)) > 2) {
      neg_test <- kruskal.test(Abs_Weight ~ Group, data = negative_edges)
      test_name <- "kw"
    } else if (length(unique(negative_edges$Group)) == 2) {
      neg_test <- wilcox.test(Abs_Weight ~ Group, data = negative_edges)
      test_name <- "U"
    } else {
      neg_test <- list(p.value = NA)
      test_name <- "N"
    }
    p_neg <- p_neg + annotate("text", x = length(unique(negative_edges$Group)) / 2 + 0.5,
                              y = max(negative_edges$Abs_Weight),
                              label = sprintf("%s test: p = %.2e", test_name, neg_test$p.value),
                              size = 4)
    saveRDS(p_neg, file = paste0(outdir, "/fig_rds/Multi_Group_Negative_Edge_Strength_Boxplot.rds"))
  } else {
    p_neg <- NULL
  }

  p_density <- ggplot(all_edges, aes(x = Absolute_Weight, color = Group, fill = Group)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_log10() +
    labs(title = "Edge Strength Distribution Across Groups",
         x = "Absolute Edge Weight (log10)",
         y = "Density",
         color = "Group",
         fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "top")

  saveRDS(p_density, file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Strength_Density.rds"))

  strength_df_list <- list()
  for (grp in group_names) {
    if (!is.null(group_edges_list[[grp]])) {
      edges <- group_edges_list[[grp]]
      strength_table <- table(edges$Strength_Category)
      temp_df <- as.data.frame.matrix(t(strength_table))
      temp_df$Group <- grp
      strength_df_list[[grp]] <- temp_df
    }
  }

  if (length(strength_df_list) > 0) {
    strength_df <- do.call(rbind, strength_df_list)
    if ("ACA" %in% group_names) {
      strength_df$Group <- factor(strength_df$Group, levels = c("Control", "CRC", "ACA"))
    } else if (eras) {
      strength_df$Group <- factor(strength_df$Group, levels = c("Control", "D4", "D8", "D42", "D180"))
    }
    strength_long <- melt(strength_df, id.vars = "Group",
                          variable.name = "Strength", value.name = "Count")

    p_strength_cat <- ggplot(strength_long, aes(x = Group, y = Count, fill = Strength)) +
      geom_bar(stat = "identity", position = "fill", width = 0.7) +
      scale_fill_brewer(palette = "RdYlGn", direction = -1) +
      labs(title = "Edge Strength Distribution Across Groups",
           x = "Group",
           y = "Proportion",
           fill = "Strength Category") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 10),
            legend.position = "top")

    saveRDS(p_strength_cat, file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Strength_Categories.rds"))
  } else {
    p_strength_cat <- NULL
  }

  pairwise_results <- list()

  if (length(group_names) >= 2) {
    edge_count_matrix <- matrix(0, nrow = length(group_names), ncol = 2)
    rownames(edge_count_matrix) <- group_names
    colnames(edge_count_matrix) <- c("Positive", "Negative")

    for (i in seq_along(group_names)) {
      grp <- group_names[i]
      if (!is.null(group_edges_list[[grp]])) {
        edges <- group_edges_list[[grp]]
        edge_count_matrix[i, 1] <- sum(edges$Edge_Type == "Positive")
        edge_count_matrix[i, 2] <- sum(edges$Edge_Type == "Negative")
      }
    }

    pairwise_edge_p <- matrix(NA, nrow = length(group_names), ncol = length(group_names))
    rownames(pairwise_edge_p) <- group_names
    colnames(pairwise_edge_p) <- group_names

    for (i in 1:(length(group_names) - 1)) {
      for (j in (i + 1):length(group_names)) {
        cont_table <- rbind(edge_count_matrix[i, ], edge_count_matrix[j, ])
        fisher_res <- fisher.test(cont_table, simulate.p.value = TRUE)
        pairwise_edge_p[i, j] <- pairwise_edge_p[j, i] <- fisher_res$p.value
      }
    }

    pairwise_results$edge_counts <- pairwise_edge_p
    write.csv(pairwise_edge_p, paste0(outdir, "data/pairwise_edge_count_pvalues.csv"))
  }

  pairwise_strength_p <- matrix(NA, nrow = length(group_names), ncol = length(group_names))
  rownames(pairwise_strength_p) <- group_names
  colnames(pairwise_strength_p) <- group_names

  for (i in seq_along(group_names)) {
    for (j in seq_along(group_names)) {
      if (i < j) {
        grp1 <- group_names[i]
        grp2 <- group_names[j]
        if (!is.null(group_edges_list[[grp1]]) && !is.null(group_edges_list[[grp2]])) {
          weights1 <- abs(group_edges_list[[grp1]]$Edge_Weight)
          weights2 <- abs(group_edges_list[[grp2]]$Edge_Weight)
          wilcox_res <- wilcox.test(weights1, weights2)
          pairwise_strength_p[i, j] <- pairwise_strength_p[j, i] <- wilcox_res$p.value
        }
      }
    }
  }

  pairwise_results$edge_strength <- pairwise_strength_p
  write.csv(pairwise_strength_p, paste0(outdir, "data/pairwise_edge_strength_pvalues.csv"))

  valid_plots <- list()
  if (!is.null(p_count)) valid_plots[["Count"]] <- p_count
  if (!is.null(p_strength)) valid_plots[["Strength"]] <- p_strength
  if (!is.null(p_pos)) valid_plots[["Positive"]] <- p_pos
  if (!is.null(p_neg)) valid_plots[["Negative"]] <- p_neg
  if (!is.null(p_density)) valid_plots[["Density"]] <- p_density
  if (!is.null(p_strength_cat)) valid_plots[["Categories"]] <- p_strength_cat

  n_plots <- length(valid_plots)
  if (n_plots >= 2) {
    if (n_plots <= 3) {
      ncol <- n_plots
      nrow <- 1
    } else if (n_plots <= 6) {
      ncol <- 3
      nrow <- ceiling(n_plots / 3)
    } else {
      ncol <- 4
      nrow <- ceiling(n_plots / 4)
    }

    combined_plots <- wrap_plots(valid_plots, ncol = ncol, nrow = nrow,
                                 guides = "collect") &
      theme(legend.position = "bottom")

    combined_plots <- combined_plots +
      plot_annotation(
        title = paste0("Comprehensive Network Analysis: ",
                       paste(group_names, collapse = " vs ")),
        subtitle = "Comparison of Edge Counts, Strengths, and Density",
        caption = "Statistical significance: Chi-square test (Counts) / Kruskal-Wallis test (Others)",
        tag_levels = "A"
      )

    ggsave(paste0(outdir, "multi_group_network_comprehensive_analysis.png"),
           combined_plots,
           width = min(20, 5 * ncol),
           height = min(15, 5 * nrow),
           dpi = 300, limitsize = FALSE)
  }

  edge_summary <- data.frame()
  for (grp in group_names) {
    if (!is.null(group_edges_list[[grp]])) {
      edges <- group_edges_list[[grp]]
      edge_summary <- rbind(edge_summary, data.frame(
        Group = grp,
        Total_Edges = ecount(graphs[[grp]]),
        Positive_Edges = sum(edges$Edge_Type == "Positive"),
        Negative_Edges = sum(edges$Edge_Type == "Negative"),
        Mean_Edge_Strength = round(mean(abs(edges$Edge_Weight)), 4),
        Median_Edge_Strength = round(median(abs(edges$Edge_Weight)), 4),
        Mean_Positive_Strength = round(mean(edges$Edge_Weight[edges$Edge_Type == "Positive"]), 4),
        Mean_Negative_Strength = round(mean(abs(edges$Edge_Weight[edges$Edge_Type == "Negative"])), 4),
        Max_Edge_Strength = round(max(abs(edges$Edge_Weight)), 4),
        Min_Edge_Strength = round(min(abs(edges$Edge_Weight)), 4),
        Rho = rho_results$Rho[rho_results$Group == grp]
      ))
    }
  }

  if (nrow(edge_summary) >= 2) {
    base_values <- edge_summary[1, ]

    for (i in 2:nrow(edge_summary)) {
      for (col in c("Total_Edges", "Positive_Edges", "Negative_Edges",
                    "Mean_Edge_Strength", "Median_Edge_Strength")) {
        change_col <- paste0(col, "_Change_%")
        edge_summary[i, change_col] <- round(
          (edge_summary[i, col] - base_values[[col]]) / base_values[[col]] * 100, 1
        )
      }
    }
  }

  write.csv(edge_summary, file = paste0(outdir, "Multi_Group_Edge_Statistics_Summary.csv"),
            row.names = FALSE)

  cat("\nEdge Statistics Summary:\n")
  print(edge_summary)

  sink(paste0(outdir, "Multi_Group_Edge_Statistics_Tests.txt"))
  cat("Multi-Group Edge Interaction Statistics\n")
  cat("========================================\n\n")

  cat("Groups analyzed:", paste(group_names, collapse = ", "), "\n\n")

  cat("Edge Count Comparison (Chi-square test):\n")
  cat("  Chi-squared =", chi_res$statistic, "\n")
  cat("  df =", chi_res$parameter, "\n")
  cat("  p-value =", chi_res$p.value, "\n\n")

  cat("Edge Strength Comparison (Kruskal-Wallis test):\n")
  print(kw_test)
  cat("\n")

  if (exists("pos_test") && !is.null(pos_test)) {
    cat("Positive Edge Strength Comparison:\n")
    print(pos_test)
    cat("\n")
  }

  if (exists("neg_test") && !is.null(neg_test)) {
    cat("Negative Edge Strength Comparison:\n")
    print(neg_test)
    cat("\n")
  }

  cat("\nPairwise Comparisons (Wilcoxon test for edge strength):\n")
  print(pairwise_strength_p)

  sink()

  list(
    all_edges = all_edges,
    group_edges = group_edges_list,
    edge_summary = edge_summary,
    rho_results = rho_results,
    pairwise_pvalues = pairwise_results,
    statistical_tests = list(
      edge_count_chi = chi_res,
      edge_strength_kw = kw_test
    )
  )
}

plot_multi_edges_connect <- function(graphs, group_names, outdir, eras = FALSE) {
  n_groups <- length(group_names)
  if (n_groups <= 3) {
    colors <- c("#f783ac", "#3bc9db", "#ffa500")
  } else {
    colors <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set1")
  }
  colors <- setNames(colors[1:n_groups], group_names)

  all_edges_classified_list <- list()
  group_edges_list <- list()

  for (grp in group_names) {
    if (!is.null(graphs[[grp]]) && vcount(graphs[[grp]]) > 0) {
      g_classified <- classify_abundance_for_edges(graphs[[grp]])
      edges_classified <- classify_edges(g_classified, grp)
      group_edges_list[[grp]] <- edges_classified
      all_edges_classified_list[[grp]] <- edges_classified
    }
  }

  rho_results <- calculate_rho_by_edgetype(group_edges_list, names(group_edges_list))
  write.csv(rho_results, paste0(outdir, "data/rho_results_edge_type.csv"), row.names = FALSE)

  if (length(all_edges_classified_list) < 2) {
    cat("[WARNING] Need at least 2 groups for comparison\n")
    return(NULL)
  }

  all_edges_classified <- do.call(rbind, all_edges_classified_list)
  all_edges_classified$Edge_Type <- factor(
    all_edges_classified$Edge_Type,
    levels = c("Rich-Rich", "Rich-Rare", "Rare-Rare")
  )
  if ("ACA" %in% group_names) {
    all_edges_classified$Group <- factor(
      all_edges_classified$Group, levels = c("Control", "CRC", "ACA")
    )
  } else if (eras) {
    all_edges_classified$Group <- factor(
      all_edges_classified$Group, levels = c("Control", "D4", "D8", "D42", "D180")
    )
  }

  edge_count_summary <- aggregate(Weight ~ Group + Edge_Type,
                                  data = all_edges_classified, FUN = length)
  colnames(edge_count_summary)[3] <- "Count"
  edge_count_summary$Edge_Type <- factor(
    edge_count_summary$Edge_Type, levels = c("Rich-Rich", "Rich-Rare", "Rare-Rare")
  )
  group_summary <- data.frame()
  for (grp in group_names) {
    if (!is.null(graphs[[grp]])) {
      group_summary <- rbind(group_summary, data.frame(
        Group = grp,
        Total_Nodes = vcount(graphs[[grp]]),
        Total_Edges = ecount(graphs[[grp]])
      ))
    }
  }

  edge_types <- levels(edge_count_summary$Edge_Type)
  p_values_matrix <- matrix(NA, nrow = length(edge_types), ncol = 1)
  rownames(p_values_matrix) <- edge_types
  colnames(p_values_matrix) <- "Global_p_value"

  result_table <- data.frame()
  for (et in edge_types) {
    counts <- sapply(group_names, function(grp) {
      cnt <- edge_count_summary$Count[
        edge_count_summary$Edge_Type == et & edge_count_summary$Group == grp
      ]
      ifelse(length(cnt) == 0, 0, cnt)
    })

    for (i in seq_along(group_names)) {
      result_table <- rbind(result_table, data.frame(
        Edge_Type = et,
        Group = group_names[i],
        Value = counts[i]
      ))
    }

    if (length(group_names) >= 2) {
      cont_matrix <- cbind(counts, group_summary$Total_Edges - counts)
      cont_matrix <- cont_matrix[
        rowSums(cont_matrix) > 0,
        colSums(cont_matrix) > 0,
        drop = FALSE
      ]
      if (nrow(cont_matrix) < 2 || ncol(cont_matrix) < 2) {
        p_values_matrix[et, "Global_p_value"] <- NA
      } else {
        fisher_res <- fisher.test(cont_matrix, simulate.p.value = TRUE, B = 10000)
        p_values_matrix[et, "Global_p_value"] <- fisher_res$p.value
      }
    }
  }

  p_labels <- edge_count_summary %>%
    group_by(Edge_Type) %>%
    summarize(max_count = max(Count), .groups = "drop")

  p_labels$p_val_text <- sprintf("p = %.2e",
                                 p_values_matrix[p_labels$Edge_Type, "Global_p_value"])

  p_count <- ggplot(edge_count_summary, aes(x = Edge_Type, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = colors) +
    labs(title = "Number of Interactions by Species Type Across Groups",
         x = "Interaction Type", y = "Edge Count", fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    geom_text(aes(label = Count), position = position_dodge(0.8), vjust = -0.5, size = 3) +
    geom_text(data = p_labels,
              aes(x = Edge_Type, y = max_count, label = p_val_text),
              inherit.aes = FALSE,
              vjust = -2.5, size = 3.5, fontface = "italic") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

  write.csv(result_table,
            paste0(outdir, "/data/Multi_Group_Edge_Count_by_Interaction_Type.csv"),
            row.names = FALSE, quote = FALSE)
  saveRDS(p_count,
          file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Count_by_Interaction_Type.rds"))

  if (eras) {
    saveRDS(edge_count_summary, file = paste0(outdir, "/data/edge_count_summary.rds"))
  }

  p_labels_strength <- data.frame()
  strength_stats <- list()

  for (et in edge_types) {
    strength_list <- list()
    for (grp in group_names) {
      vals <- all_edges_classified$Abs_Weight[
        all_edges_classified$Group == grp & all_edges_classified$Edge_Type == et
      ]
      if (length(vals) > 0) {
        strength_list[[grp]] <- vals
      } else {
        strength_list[[grp]] <- 0
      }
    }

    strength_df <- data.frame()
    for (grp in names(strength_list)) {
      strength_df <- rbind(strength_df, data.frame(
        Group = grp,
        Abs_Weight = strength_list[[grp]]
      ))
    }
    if (length(strength_list) > 2) {
      w_res <- kruskal.test(Abs_Weight ~ Group, data = strength_df)
    } else {
      w_res <- wilcox.test(Abs_Weight ~ Group, data = strength_df)
    }
    global_p <- w_res$p.value

    pairwise_p <- matrix(NA, nrow = length(group_names), ncol = length(group_names))
    rownames(pairwise_p) <- group_names
    colnames(pairwise_p) <- group_names

    for (i in 1:(length(group_names) - 1)) {
      for (j in (i + 1):length(group_names)) {
        grp1 <- group_names[i]
        grp2 <- group_names[j]
        if (grp1 %in% names(strength_list) && grp2 %in% names(strength_list)) {
          wilcox_res <- wilcox.test(strength_list[[grp1]], strength_list[[grp2]])
          pairwise_p[i, j] <- pairwise_p[j, i] <- wilcox_res$p.value
        }
      }
    }

    strength_stats[[et]] <- list(global_p = global_p, pairwise = pairwise_p)

    max_val <- max(all_edges_classified$Abs_Weight[all_edges_classified$Edge_Type == et],
                   na.rm = TRUE)
    p_labels_strength <- rbind(p_labels_strength, data.frame(
      Edge_Type = et,
      max_y = max_val,
      p_text = sprintf("Global p = %.2e", global_p)
    ))
  }

  p_strength <- ggplot(all_edges_classified, aes(x = Edge_Type, y = Abs_Weight, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    scale_y_log10() +
    labs(title = "Interaction Strength by Species Type Across Groups",
         x = "Interaction Type",
         y = "Absolute Edge Weight (log10)",
         fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    annotation_logticks(sides = "l") +
    geom_text(data = p_labels_strength,
              aes(x = Edge_Type, y = max_y, label = p_text),
              inherit.aes = FALSE,
              vjust = -1.2, size = 3, fontface = "italic") +
    expand_limits(y = max(all_edges_classified$Abs_Weight) * 1.5)

  saveRDS(p_strength,
          file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Strength_by_Interaction_Type.rds"))

  all_edges_classified$Sign <- ifelse(all_edges_classified$Weight > 0, "Positive", "Negative")

  sign_summary <- aggregate(Weight ~ Group + Edge_Type + Sign,
                            data = all_edges_classified, FUN = length)
  colnames(sign_summary)[4] <- "Count"

  p_sign <- ggplot(sign_summary, aes(x = Edge_Type, y = Count, fill = Sign)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    facet_wrap(~Group, nrow = 1) +
    scale_fill_manual(values = c("Positive" = "#ac0945", "Negative" = "#36de15")) +
    labs(title = "Proportion of Positive/Negative Interactions by Species Type Across Groups",
         x = "Interaction Type",
         y = "Proportion",
         fill = "Interaction Sign") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          strip.background = element_rect(fill = "lightgray"),
          strip.text = element_text(size = 11, face = "bold"))

  saveRDS(p_sign, file = paste0(outdir, "/fig_rds/Multi_Group_Edge_Sign_by_Interaction_Type.rds"))

  plot_width <- 20
  plot_height <- 6
  combined_main_plots <- p_count | p_strength | p_sign

  combined_final <- combined_main_plots +
    plot_annotation(
      title = paste0("Key Network Characteristics: ",
                     paste(group_names, collapse = " vs ")),
      subtitle = "Comparing Interaction Counts, Strengths, and Positive/Negative Proportions",
      tag_levels = "A"
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(paste0(outdir, "Multi_Group_Network_edge_connect_Combined.png"),
         combined_final, width = plot_width, height = plot_height, dpi = 300)

  interaction_summary <- aggregate(Abs_Weight ~ Group + Edge_Type,
                                   data = all_edges_classified,
                                   FUN = function(x) c(Mean = mean(x), Median = median(x),
                                                       SD = sd(x), N = length(x)))
  interaction_summary <- do.call(data.frame, interaction_summary)
  colnames(interaction_summary) <- c("Group", "Interaction_Type", "Mean_Strength",
                                     "Median_Strength", "SD_Strength", "Count")

  write.csv(interaction_summary,
            file = paste0(outdir, "Multi_Group_Interaction_Type_Summary.csv"),
            row.names = FALSE)

  cat("\nInteraction Type Summary:\n")
  print(interaction_summary)

  sink(paste0(outdir, "Multi_Group_Interaction_Type_Statistics.txt"))
  cat("Multi-Group Rich-Rare Interaction Analysis\n")
  cat("===========================================\n\n")
  cat("Groups analyzed:", paste(group_names, collapse = ", "), "\n\n")

  for (et in edge_types) {
    cat(sprintf("\n>>> %s <<<\n", et))
    cat(sprintf("  Global comparison (Kruskal-Wallis): p = %.2e\n",
                strength_stats[[et]]$global_p))

    if (!is.null(strength_stats[[et]]$pairwise)) {
      cat("  Pairwise comparisons (Wilcoxon):\n")
      pairwise_mat <- strength_stats[[et]]$pairwise
      for (i in 1:(nrow(pairwise_mat) - 1)) {
        for (j in (i + 1):ncol(pairwise_mat)) {
          if (!is.na(pairwise_mat[i, j])) {
            cat(sprintf("    %s vs %s: p = %.2e\n",
                        rownames(pairwise_mat)[i], colnames(pairwise_mat)[j],
                        pairwise_mat[i, j]))
          }
        }
      }
    }
    cat("--------------------------------------------------\n")
  }

  cat("\n\nEdge Count Proportion Analysis (Fisher's exact test):\n")
  for (et in edge_types) {
    cat(sprintf("  %s: p = %.2e\n", et, p_values_matrix[et, "Global_p_value"]))
  }

  sink()

  list(
    edge_count_summary = edge_count_summary,
    interaction_summary = interaction_summary,
    strength_stats = strength_stats,
    edge_type_pvalues = p_values_matrix
  )
}
