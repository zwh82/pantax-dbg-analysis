# Figures 5 and 6: rho (net interaction) Control vs Disease.
#
# 05a  <results>/rho_compare/<disease>/<suffix>/rho_boxplot_connected_*.png
#      <results>/rho_compare/combined/<suffix>_all_plots_combined.png
#      <results>/rho_compare/combined/<suffix>_single_cohort_combined.png
#      <results>/rho_compare/combined/<suffix>multiple_cohort_combined.png
#      <results>/scripts/final_figures/fig_net_5.pdf
# 05b  <results>/rho_compare/<disease>/<suffix>/rho_boxplot_facet_<disease>.png
#      <results>/rho_compare/combined/<suffix>_multiple_cohort_faceted.png
#      <results>/scripts/final_figures/fig_net_6.pdf

run_rho_global <- function(cfg) {
  disease_map <- read_disease_map(cfg)
  suffixes <- analysis_suffixes(cfg)

  for (outdir_last in suffixes) {
    plot_list <- list()
    plot_list_single <- list()
    plot_list_multiple <- list()
    disease_cohorts <- get_disease_cohorts(network_dir(cfg))

    for (disease in names(disease_cohorts)) {
      if (disease %in% cfg$skip_diseases$rho) next

      all_rho_data <- data.frame()
      cohorts <- disease_cohorts[[disease]]
      output_dir <- ensure_dir(
        file.path(cfg$results_dir, cfg$paths$rho, disease, outdir_last)
      )

      for (cohort in cohorts) {
        rho_file <- file.path(
          network_dir(cfg), cohort, outdir_last, "data", "rho_results_multi.csv"
        )
        if (!file.exists(rho_file)) next

        rho_data <- read.csv(rho_file, stringsAsFactors = FALSE)
        control_val <- rho_data %>%
          filter(grepl("Control", Group, ignore.case = TRUE)) %>%
          pull(Rho)
        if (length(control_val) == 0) next

        disease_rows <- rho_data %>%
          filter(!grepl("Control", Group, ignore.case = TRUE))
        if (nrow(disease_rows) == 0) next

        if (disease == "IBD") {
          for (i in seq_len(nrow(disease_rows))) {
            disease_row <- disease_rows[i, ]
            all_rho_data <- bind_rows(
              all_rho_data,
              data.frame(
                Cohort = cohort,
                Disease = "IBD",
                Outdir = outdir_last,
                Group = c("Control", "IBD"),
                Rho = c(control_val, disease_row$Rho),
                Pair_ID = paste(cohort, "Control", "IBD", sep = "_"),
                stringsAsFactors = FALSE
              )
            )
          }
        } else {
          disease_row <- pick_disease_row(disease_rows, disease, cohort)
          all_rho_data <- bind_rows(
            all_rho_data,
            data.frame(
              Cohort = cohort,
              Disease = disease,
              Outdir = outdir_last,
              Group = c("Control", disease_row$Group),
              Rho = c(control_val, disease_row$Rho),
              Pair_ID = paste(cohort, "Control", disease_row$Group, sep = "_"),
              stringsAsFactors = FALSE
            )
          )
        }
      }

      all_rho_data <- replace_disease_name2(
        all_rho_data, disease_map, cols = c("Disease", "Group")
      )
      if (nrow(all_rho_data) == 0) {
        message("No rho data found for ", disease)
        next
      }

      rho_df <- all_rho_data
      rho_df$Group <- factor(
        rho_df$Group,
        levels = c("Control", unique(rho_df$Group[rho_df$Group != "Control"]))
      )

      p1 <- ggplot(rho_df, aes(x = Group, y = Rho, fill = Group)) +
        geom_boxplot(alpha = 0.7, width = 0.6,
                     outlier.shape = 21, outlier.size = 1.5) +
        geom_line(
          data = rho_df %>% group_by(Pair_ID) %>% filter(n() == 2),
          aes(group = Pair_ID),
          color = "gray50", alpha = 0.6, linewidth = 0.5
        ) +
        geom_point(
          data = rho_df %>% group_by(Pair_ID) %>% filter(n() == 2),
          aes(color = Group), size = 2, alpha = 0.8
        ) +
        labs(x = "", y = "Net interaction") +
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "none"
        )

      ggsave(
        file.path(output_dir,
                  paste0("rho_boxplot_connected_", paste(disease, collapse = "_"), ".png")),
        p1, width = 8, height = 6
      )

      plot_list[[length(plot_list) + 1]] <- p1
      if (length(cohorts) == 1) {
        plot_list_single[[length(plot_list_single) + 1]] <- p1
      } else {
        plot_list_multiple[[length(plot_list_multiple) + 1]] <- p1
      }
    }

    n_plots <- length(plot_list)
    n_col <- 5
    n_row <- ceiling(n_plots / n_col)
    combined_plot <- wrap_plots(plot_list, ncol = n_col, nrow = n_row) +
      plot_annotation(
        title = "Rho Comparison: Control vs Disease Across All Cohorts",
        subtitle = paste("Total plots:", n_plots, "| Methods:",
                         paste(outdir_last, collapse = ", ")),
        theme = theme(
          plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )

    combined_output_dir <- ensure_dir(
      file.path(cfg$results_dir, cfg$paths$rho, "combined")
    )
    ggsave(
      file.path(combined_output_dir, paste0(outdir_last, "_all_plots_combined.png")),
      combined_plot, width = 20, height = 24, dpi = 300, limitsize = FALSE
    )

    if (length(plot_list_single) > 0) {
      n_plots_single <- length(plot_list_single)
      n_col_single <- ceiling(sqrt(n_plots_single))
      n_row_single <- ceiling(n_plots_single / n_col_single)
      combined_single <- wrap_plots(
        plot_list_single, ncol = n_col_single, nrow = n_row_single
      ) +
        plot_annotation(
          title = "Single Cohort Diseases: Rho Comparison (Control vs Disease)",
          subtitle = paste("Total:", n_plots_single,
                           "diseases | Each with only 1 cohort"),
          theme = theme(
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 12, hjust = 0.5)
          )
        )
      ggsave(
        file.path(combined_output_dir,
                  paste0(outdir_last, "_single_cohort_combined.png")),
        combined_single, width = 16, height = 16, dpi = 300, limitsize = FALSE
      )
    }

    if (length(plot_list_multiple) > 0) {
      n_plots_multiple <- length(plot_list_multiple)
      n_col_multiple <- 5
      n_row_multiple <- ceiling(n_plots_multiple / n_col_multiple)
      plot_list_multiple_title <- keep_first_title_each_row(
        plot_list_multiple, n_col_multiple
      )
      combined_multiple <- wrap_plots(
        plot_list_multiple_title, ncol = n_col_multiple, nrow = n_row_multiple
      ) +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(size = 16, face = "bold")) +
        plot_layout(heights = c(1, 0.3, 1))

      ggsave(
        file.path(combined_output_dir,
                  paste0(outdir_last, "multiple_cohort_combined.png")),
        combined_multiple, width = 12, height = 10, dpi = 300, limitsize = FALSE
      )
      ggsave(file.path(cfg$results_dir, cfg$paths$final, "fig_net_5.pdf"),
             combined_multiple, width = 12, height = 10)
    }
  }

  invisible(TRUE)
}

run_rho_by_edge_type <- function(cfg) {
  disease_map <- read_disease_map(cfg)
  suffixes <- analysis_suffixes(cfg)

  for (outdir_last in suffixes) {
    plot_list_multiple <- list()
    disease_cohorts <- get_disease_cohorts(network_dir(cfg))

    combined_output_dir <- ensure_dir(
      file.path(cfg$results_dir, cfg$paths$rho, "combined")
    )

    for (disease in names(disease_cohorts)) {
      if (disease %in% cfg$skip_diseases$rho_edge) next

      all_rho_data <- data.frame()
      cohorts <- disease_cohorts[[disease]]
      output_dir <- ensure_dir(
        file.path(cfg$results_dir, cfg$paths$rho, disease, outdir_last)
      )

      for (cohort in cohorts) {
        rho_file <- file.path(
          network_dir(cfg), cohort, outdir_last, "data", "rho_results_edge_type.csv"
        )
        if (!file.exists(rho_file)) next

        rho_data <- read.csv(rho_file, stringsAsFactors = FALSE)
        control_row <- rho_data %>%
          filter(grepl("Control", Group, ignore.case = TRUE))
        if (nrow(control_row) == 0) next

        disease_rows <- rho_data %>%
          filter(!grepl("Control", Group, ignore.case = TRUE))
        if (nrow(disease_rows) == 0) next

        if (disease == "IBD") {
          for (i in seq_len(nrow(disease_rows))) {
            disease_row <- disease_rows[i, ]
            for (edge_type in c("Rich_Rich", "Rich_Rare", "Rare_Rare")) {
              all_rho_data <- bind_rows(
                all_rho_data,
                data.frame(
                  Cohort = cohort,
                  Disease = "IBD",
                  Subtype = disease_row$Group,
                  Outdir = outdir_last,
                  EdgeType = edge_type,
                  Group = c("Control", disease_row$Group),
                  Rho = c(control_row[[edge_type]], disease_row[[edge_type]]),
                  Pair_ID = paste(cohort, "Control", disease_row$Group,
                                  edge_type, sep = "_"),
                  stringsAsFactors = FALSE
                )
              )
            }
          }
        } else {
          disease_row <- pick_disease_row(disease_rows, disease, cohort)
          for (edge_type in c("Rich_Rich", "Rich_Rare", "Rare_Rare")) {
            all_rho_data <- bind_rows(
              all_rho_data,
              data.frame(
                Cohort = cohort,
                Disease = disease,
                Subtype = NA,
                Outdir = outdir_last,
                EdgeType = edge_type,
                Group = c("Control", disease_row$Group),
                Rho = c(control_row[[edge_type]], disease_row[[edge_type]]),
                Pair_ID = paste(cohort, "Control", disease_row$Group,
                                edge_type, sep = "_"),
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }

      all_rho_data <- replace_disease_name2(
        all_rho_data, disease_map, cols = c("Disease", "Group")
      )
      if (nrow(all_rho_data) == 0) {
        message("No edge-type rho data found for ", disease)
        next
      }

      rho_df <- all_rho_data %>%
        mutate(EdgeType = rename_rich_to_abundant(EdgeType))
      rho_df$Group <- factor(
        rho_df$Group,
        levels = c("Control", unique(rho_df$Group[rho_df$Group != "Control"]))
      )

      edge_labels <- c(
        "Abundant_Abundant" = "Abundant-Abundant",
        "Abundant_Rare" = "Abundant-Rare",
        "Rare_Rare" = "Rare-Rare"
      )
      rho_df$EdgeType_Label <- edge_labels[rho_df$EdgeType]
      rho_df$EdgeType_Label <- factor(
        rho_df$EdgeType_Label,
        levels = c("Rare-Rare", "Abundant-Rare", "Abundant-Abundant")
      )

      p1 <- ggplot(rho_df, aes(x = Group, y = Rho, fill = Group)) +
        geom_boxplot(alpha = 0.7, width = 0.6,
                     outlier.shape = 21, outlier.size = 1.5) +
        geom_line(
          data = rho_df %>% group_by(Pair_ID) %>% filter(n() == 2),
          aes(group = Pair_ID),
          color = "gray50", alpha = 0.6, linewidth = 0.5
        ) +
        geom_point(
          data = rho_df %>% group_by(Pair_ID) %>% filter(n() == 2),
          aes(color = Group), size = 2, alpha = 0.8
        ) +
        facet_wrap(~EdgeType_Label, ncol = 3, scales = "free_y") +
        labs(x = "", y = "Net interaction") +
        theme_bw() +
        theme(
          legend.position = "none",
          strip.background = element_rect(fill = "lightgray", color = NA),
          strip.text = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
        )

      ggsave(file.path(output_dir, paste0("rho_boxplot_facet_", disease, ".png")),
             p1, width = 12, height = 6)

      if (length(unique(rho_df$Cohort)) > 1) {
        plot_list_multiple[[length(plot_list_multiple) + 1]] <- p1
      }

      for (edge_type in unique(rho_df$EdgeType_Label)) {
        sub_df <- rho_df %>% filter(EdgeType_Label == edge_type)
        control_rho <- sub_df %>% filter(Group == "Control") %>% pull(Rho)
        disease_rho <- sub_df %>% filter(Group != "Control") %>% pull(Rho)
        if (length(control_rho) > 0 && length(disease_rho) > 0) {
          wilcox_result <- wilcox.test(control_rho, disease_rho, paired = TRUE)
          message(sprintf("%s: paired Wilcoxon p = %.4f",
                          edge_type, wilcox_result$p.value))
        }
      }
    }

    if (length(plot_list_multiple) > 0) {
      n_col_multiple <- 3
      n_row_multiple <- ceiling(length(plot_list_multiple) / n_col_multiple)
      plot_list_multiple_title <- keep_first_title_each_row(
        plot_list_multiple, n_col_multiple
      )
      combined_multiple <- wrap_plots(
        plot_list_multiple_title, ncol = n_col_multiple, nrow = n_row_multiple
      ) +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(size = 16, face = "bold")) +
        plot_layout(heights = c(1, 0.3, 1))

      ggsave(
        file.path(combined_output_dir,
                  paste0(outdir_last, "_multiple_cohort_faceted.png")),
        combined_multiple, width = 18, height = 18, dpi = 300, limitsize = FALSE
      )
      ggsave(file.path(cfg$results_dir, cfg$paths$final, "fig_net_6.pdf"),
             combined_multiple, width = 18, height = 18)
    }
  }

  invisible(TRUE)
}
